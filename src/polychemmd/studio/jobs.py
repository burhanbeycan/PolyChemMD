"""SQLite job ownership and immutable successful attempts for local computations."""

from __future__ import annotations

from contextlib import closing
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
import threading
import time
import uuid

KINDS = ("geometry", "packing", "md", "radonpy", "experiment", "trajectory")
LEASE_SECONDS = 120
MAX_ATTEMPTS = 3


def stamp():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def encode(value):
    return json.dumps(
        value,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )


def file_digest(file):
    digest = hashlib.sha256()
    with file.open("rb") as handle:
        for piece in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(piece)
    return digest.hexdigest()


class JobStore:
    def __init__(self, root):
        self.root = Path(root).expanduser().resolve()
        self.root.mkdir(parents=True, exist_ok=True)
        self.db = self.root / "jobs.sqlite"
        with closing(self.connect()) as db:
            db.execute("PRAGMA journal_mode=WAL")
            db.executescript("""
            CREATE TABLE IF NOT EXISTS jobs (
              id TEXT PRIMARY KEY, kind TEXT NOT NULL, config TEXT NOT NULL,
              status TEXT NOT NULL, created_at TEXT NOT NULL, updated_at TEXT NOT NULL,
              attempts INTEGER NOT NULL DEFAULT 0, owner TEXT, expires REAL,
              result TEXT, error TEXT, files TEXT NOT NULL DEFAULT '[]', output TEXT,
              request_key TEXT UNIQUE
            );
            CREATE INDEX IF NOT EXISTS by_status ON jobs(status,created_at);
            """)

    def connect(self):
        db = sqlite3.connect(self.db, timeout=20, isolation_level=None)
        db.row_factory = sqlite3.Row
        return db

    def submit(self, kind, config, request_key=None):
        if kind not in KINDS or not isinstance(config, dict):
            raise ValueError(
                "A supported job kind and a configuration object are required."
            )
        data = encode(config)
        if len(data.encode("utf-8")) > 100_000:
            raise ValueError("A recipe must be smaller than 100 kB.")
        if request_key is not None and (
            not isinstance(request_key, str) or not 1 <= len(request_key) <= 128
        ):
            raise ValueError("Invalid idempotency key.")
        moment = stamp()
        job_id = uuid.uuid4().hex
        with closing(self.connect()) as db:
            db.execute("BEGIN IMMEDIATE")
            if request_key:
                old = db.execute(
                    "SELECT id,kind,config FROM jobs WHERE request_key=?",
                    (request_key,),
                ).fetchone()
                if old:
                    if old["kind"] != kind or old["config"] != data:
                        raise ValueError(
                            "The request key already belongs to a different recipe."
                        )
                    db.commit()
                    return self.get(old["id"])
            db.execute(
                "INSERT INTO jobs(id,kind,config,status,created_at,updated_at,request_key) VALUES(?,?,?,'queued',?,?,?)",
                (job_id, kind, data, moment, moment, request_key),
            )
            db.commit()
        return self.get(job_id)

    def get(self, job_id):
        with closing(self.connect()) as db:
            row = db.execute("SELECT * FROM jobs WHERE id=?", (job_id,)).fetchone()
        if row is None:
            raise KeyError("Unknown job")
        result = {
            k: row[k]
            for k in (
                "id",
                "kind",
                "status",
                "created_at",
                "updated_at",
                "attempts",
                "error",
            )
        }
        for k in ("config", "result", "files"):
            result[k] = json.loads(row[k]) if row[k] is not None else None
        return result

    def list(self, limit=100):
        with closing(self.connect()) as db:
            records = db.execute(
                "SELECT id,kind,status,created_at,updated_at,attempts,error FROM jobs ORDER BY rowid DESC LIMIT ?",
                (max(1, min(limit, 500)),),
            ).fetchall()
        return [dict(row) for row in records]

    def claim(self, job_id=None):
        with closing(self.connect()) as db:
            db.execute("BEGIN IMMEDIATE")
            now = time.time()
            db.execute(
                "UPDATE jobs SET status='failed',error='Worker lease expired after maximum attempts',updated_at=?,owner=NULL,expires=NULL WHERE status='running' AND expires<? AND attempts>=?",
                (stamp(), now, MAX_ATTEMPTS),
            )
            where = (
                "(status='queued' OR (status='running' AND expires<? AND attempts<?))"
            )
            params = [now, MAX_ATTEMPTS]
            if job_id is not None:
                where += " AND id=?"
                params.append(job_id)
            job = db.execute(
                "SELECT * FROM jobs WHERE " + where + " ORDER BY rowid LIMIT 1", params
            ).fetchone()
            if job is None:
                db.commit()
                return None
            owner = uuid.uuid4().hex
            db.execute(
                "UPDATE jobs SET status='running',attempts=attempts+1,owner=?,expires=?,updated_at=? WHERE id=?",
                (owner, now + LEASE_SECONDS, stamp(), job["id"]),
            )
            db.commit()
            return job["id"], owner, job["kind"], json.loads(job["config"])

    def renew(self, job_id, owner):
        with closing(self.connect()) as db:
            return bool(
                db.execute(
                    "UPDATE jobs SET expires=? WHERE id=? AND status='running' AND owner=?",
                    (time.time() + LEASE_SECONDS, job_id, owner),
                ).rowcount
            )

    def output_dir(self, job_id):
        with closing(self.connect()) as db:
            row = db.execute(
                "SELECT output FROM jobs WHERE id=? AND status='succeeded'", (job_id,)
            ).fetchone()
        if not row or not row["output"]:
            raise ValueError("A completed source job is required.")
        return self.root / row["output"]

    def finish(self, job_id, owner, *, result=None, error=None, output=None, files=()):
        with closing(self.connect()) as db:
            return bool(
                db.execute(
                    "UPDATE jobs SET status=?,result=?,error=?,output=?,files=?,updated_at=?,owner=NULL,expires=NULL WHERE id=? AND owner=? AND status='running'",
                    (
                        "failed" if error is not None else "succeeded",
                        encode(result),
                        error,
                        output,
                        encode(files),
                        stamp(),
                        job_id,
                        owner,
                    ),
                ).rowcount
            )


def execute_next(store, job_id=None):
    from .operations import execute, versions

    claim = store.claim(job_id)
    if claim is None:
        return False
    job_id, owner, kind, config = claim
    directory = store.root / "runs" / job_id / owner
    stopped = threading.Event()

    def heartbeat():
        while not stopped.wait(LEASE_SECONDS / 4):
            if not store.renew(job_id, owner):
                return

    monitor = threading.Thread(target=heartbeat, daemon=True)
    monitor.start()
    try:
        directory.mkdir(parents=True)
        (directory / "recipe.json").write_text(
            encode({"kind": kind, "config": config}), encoding="utf-8"
        )
        started = stamp()
        result = execute(kind, config, directory, store)
        (directory / "result.json").write_text(encode(result), encoding="utf-8")
        artifacts = sorted(
            p for p in directory.rglob("*") if p.is_file() and not p.is_symlink()
        )
        manifest = {
            "job": job_id,
            "kind": kind,
            "recipe_sha256": file_digest(directory / "recipe.json"),
            "started_utc": started,
            "finished_utc": stamp(),
            "software": versions(),
            "files": {
                p.relative_to(directory).as_posix(): file_digest(p) for p in artifacts
            },
        }
        (directory / "manifest.json").write_text(encode(manifest), encoding="utf-8")
        files = sorted(
            [p.relative_to(directory).as_posix() for p in artifacts] + ["manifest.json"]
        )
        store.finish(
            job_id,
            owner,
            result=result,
            output=str(directory.relative_to(store.root)),
            files=files,
        )
    except Exception as err:
        store.finish(job_id, owner, error=f"{type(err).__name__}: {err}"[:4000])
    finally:
        stopped.set()
        monitor.join()
    return True


def work(store, stop=None, poll=1.0):
    stop = stop or threading.Event()
    while not stop.is_set():
        if not execute_next(store):
            stop.wait(poll)
