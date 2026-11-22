"""Matched campaigns with shared starting seeds and an explicit paired estimand.

For the two packing comparisons, hydration is a geometric property of a
starting configuration. For fixed-substrate MD, it is averaged over sampled
frames including step zero. Neither branch is a long-time adsorption assay.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import uuid
from contextlib import closing
from random import Random
from statistics import mean, stdev

from .jobs import encode, stamp
from .planner import PAIRS

_ARM_PARAMETER = {
    "termination": "pattern",
    "hydrolysis": "acetate_pattern",
    "model": "md_profile",
}
_METRICS = ("OH", "F")
_CUTOFF_A = 3.5
# The version enters the stored recipe. A change to the estimand or experimental
# design can never silently substitute for a previous campaign under the same ID.
_VERSION = 1


def _integer(value, name, lo, hi):
    if type(value) is not int or not lo <= value <= hi:
        raise ValueError(f"{name} must be an integer in [{lo},{hi}]")
    return value


def _canonical(config):
    if not isinstance(config, dict) or set(config) - {
        "objective",
        "replicates",
        "base_seed",
        "md_steps",
        "experiment_id",
    }:
        raise ValueError("Unsupported campaign options")
    obj = config.get("objective", "termination")
    if not isinstance(obj, str) or obj not in _ARM_PARAMETER or obj not in PAIRS:
        raise ValueError("Unknown paired campaign objective")
    label = config.get("experiment_id")
    if label is not None and (
        not isinstance(label, str)
        or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,79}", label)
    ):
        raise ValueError(
            "experiment_id must be 1–80 ASCII letters, numbers, dots, underscores or dashes"
        )
    recipe = {
        "version": _VERSION,
        "objective": obj,
        "experiment_id": label,
        "replicates": _integer(config.get("replicates", 3), "replicates", 1, 12),
        "base_seed": _integer(config.get("base_seed", 17), "base_seed", 1, 900_000_000),
    }
    if obj == "model":
        recipe["md_steps"] = _integer(config.get("md_steps", 20), "md_steps", 1, 2000)
    elif "md_steps" in config:
        raise ValueError("md_steps applies only to the model comparison")
    return recipe


def _tables(db):
    db.execute("""CREATE TABLE IF NOT EXISTS studio_campaigns (
               id TEXT PRIMARY KEY, recipe TEXT NOT NULL, created_at TEXT NOT NULL)""")
    db.execute("""CREATE TABLE IF NOT EXISTS studio_campaign_jobs (
               campaign_id TEXT NOT NULL, replicate INTEGER NOT NULL, seed INTEGER NOT NULL,
               arm_position INTEGER NOT NULL, arm TEXT NOT NULL, job_id TEXT NOT NULL UNIQUE,
               PRIMARY KEY(campaign_id,replicate,arm_position))""")


def _seeds(first, n):
    rng = Random(first)
    seeds = [first]
    while len(seeds) < n:
        candidate = rng.randint(1, 900_000_000)
        if candidate not in seeds:
            seeds.append(candidate)
    return seeds


def submit_campaign(store, config):
    """Atomically enqueue both arms of each replicate; exact retries reuse all IDs."""
    recipe = _canonical(config)
    obj = recipe["objective"]
    arms = PAIRS[obj]
    if len(arms) != 2:
        raise ValueError("Paired design requires exactly two registered arms")
    kind = "md" if obj == "model" else "packing"
    # Freeze defaults for packing in the submission, so later software versions
    # do not change the intended dose, cell, salt, or water content on replay.
    from .operations import DEFAULTS

    arm_recipes = [
        {**({} if kind == "md" else DEFAULTS), **arm["config"]} for arm in arms
    ]
    recipe["arms"] = [
        {"arm": arm[_ARM_PARAMETER[obj]], "kind": kind, "config": arm}
        for arm in arm_recipes
    ]
    data = encode(recipe)
    key = (
        ("label:" + recipe["experiment_id"])
        if recipe["experiment_id"]
        else "recipe:" + data
    )
    campaign_id = "cmp_" + hashlib.sha256(key.encode()).hexdigest()[:32]
    moment = stamp()
    with closing(store.connect()) as db, db:
        db.execute("BEGIN IMMEDIATE")
        _tables(db)
        existing = db.execute(
            "SELECT recipe FROM studio_campaigns WHERE id=?", (campaign_id,)
        ).fetchone()
        if existing is not None:
            if existing["recipe"] != data:
                raise ValueError(
                    "This experiment_id belongs to a different campaign design"
                )
        else:
            db.execute(
                "INSERT INTO studio_campaigns VALUES(?,?,?)",
                (campaign_id, data, moment),
            )
            for rep, seed in enumerate(
                _seeds(recipe["base_seed"], recipe["replicates"]), 1
            ):
                for index, arm in enumerate(recipe["arms"]):
                    spec = {**arm["config"], "seed": seed}
                    if kind == "md":
                        spec["md_steps"] = recipe["md_steps"]
                    job = uuid.uuid4().hex
                    # Separate relationship table reserves campaign membership. An
                    # ordinary job submitted with any request_key cannot join it.
                    db.execute(
                        """INSERT INTO jobs(id,kind,config,status,created_at,updated_at)
                               VALUES (?,?,?,'queued',?,?)""",
                        (job, kind, encode(spec), moment, moment),
                    )
                    db.execute(
                        "INSERT INTO studio_campaign_jobs VALUES (?,?,?,?,?,?)",
                        (campaign_id, rep, seed, index, arm["arm"], job),
                    )
    return campaign_summary(store, campaign_id)


def _packing_contacts(path):
    """Counts pairs, normalized per exposed site, with lateral minimum images."""
    import numpy as np

    doc = json.loads(path.read_text(encoding="utf-8"))
    atoms = doc["atoms"]
    box = np.asarray(doc["box_A"], dtype=float)
    if box.shape != (3,) or not np.isfinite(box).all() or not (box > 0).all():
        raise ValueError("Packed coordinate artifact has an invalid cell")
    sheets = {a["sheet_index"] for a in atoms if a["role"] == "slab"}
    if sheets not in ({0}, {0, 1}):
        raise ValueError("Campaign requires a supported single-surface or slit model")

    def exposed(a):
        return a["role"] == "slab" and (
            (a["sheet_index"] == 0 and a["layer"].startswith("top_"))
            or (a["sheet_index"] == 1 and a["layer"].startswith("bottom_"))
        )

    sites = [a for a in atoms if exposed(a) and a["element"] in ("O", "F")]
    waters = np.asarray(
        [
            [a[k] for k in ("x", "y", "z")]
            for a in atoms
            if a["role"] == "water" and a["element"] == "O"
        ],
        dtype=float,
    ).reshape((-1, 3))
    result = {}
    for label, el in (("OH", "O"), ("F", "F")):
        selected = [[a["x"], a["y"], a["z"]] for a in sites if a["element"] == el]
        if not selected:
            continue
        contacts = 0
        for site in selected:
            dr = waters - np.asarray(site)
            dr[:, :2] -= box[:2] * np.rint(dr[:, :2] / box[:2])
            contacts += int(
                np.count_nonzero(np.einsum("ij,ij->i", dr, dr) <= _CUTOFF_A**2)
            )
        result[label] = contacts / len(selected)
    return result


def _job_measurements(store, row):
    if row["kind"] == "md":
        data = json.loads(row["result"] or "{}").get("water_contacts", {})
    else:
        data = _packing_contacts(store.output_dir(row["job_id"]) / "packed.json")
    # Missing site types have no denominator; omit them. Do not treat them as zero.
    if not isinstance(data, dict):
        raise ValueError("Unavailable contact measurements")
    for value in data.values():
        if type(value) not in (float, int) or not math.isfinite(value) or value < 0:
            raise ValueError("Invalid water contact count")
    return {key: data[key] for key in _METRICS if key in data}


def campaign_summary(store, campaign_id):
    """Read an atomic status snapshot; aggregate second-minus-first complete pairs."""
    if not isinstance(campaign_id, str) or not re.fullmatch(
        r"cmp_[0-9a-f]{32}", campaign_id
    ):
        raise KeyError("Unknown campaign")
    with closing(store.connect()) as db, db:
        db.execute("BEGIN")
        schema = db.execute(
            "SELECT 1 FROM sqlite_master WHERE name='studio_campaigns'"
        ).fetchone()
        if not schema:
            raise KeyError("Unknown campaign")
        record = db.execute(
            "SELECT * FROM studio_campaigns WHERE id=?", (campaign_id,)
        ).fetchone()
        if record is None:
            raise KeyError("Unknown campaign")
        rows = db.execute(
            """SELECT x.*,j.kind,j.status,j.error,j.result FROM studio_campaign_jobs x
              JOIN jobs j ON x.job_id=j.id WHERE x.campaign_id=? ORDER BY replicate,arm_position""",
            (campaign_id,),
        ).fetchall()
    recipe = json.loads(record["recipe"])
    results = []
    differences = {k: [] for k in _METRICS}
    complete = 0
    for first, second in zip(rows[::2], rows[1::2], strict=True):
        if first["replicate"] != second["replicate"] or first["seed"] != second["seed"]:
            raise ValueError("Invalid campaign pairing")
        jobs = [
            {"id": r["job_id"], **{k: r[k] for k in ("kind", "arm", "status", "error")}}
            for r in (first, second)
        ]
        pair = {
            "replicate": first["replicate"],
            "seed": first["seed"],
            "jobs": jobs,
            "values": None,
            "difference": None,
        }
        if all(j["status"] == "succeeded" for j in jobs):
            try:
                a, b = (_job_measurements(store, r) for r in (first, second))
                if not any(k in a and k in b for k in _METRICS):
                    raise ValueError("No site class has a matched hydration observable")
                pair["values"] = [a, b]
                pair["difference"] = {
                    k: b[k] - a[k] for k in _METRICS if k in a and k in b
                }
                for k, d in pair["difference"].items():
                    differences[k].append(d)
                complete += 1
            except (OSError, ValueError, KeyError, IndexError, TypeError) as exc:
                pair["analysis_error"] = str(exc)
        results.append(pair)
    states = {r["status"] for r in rows}
    failed = bool(
        "failed" in states or any("analysis_error" in pair for pair in results)
    )
    if states == {"queued"}:
        state = "queued"
    elif states & {"queued", "running"}:
        state = "running"
    elif failed:
        state = "partial_failed" if complete else "failed"
    else:
        state = "succeeded"
    metrics = {
        k: {
            "n": len(d),
            "mean_difference": mean(d) if d else None,
            "sd_difference": stdev(d) if len(d) > 1 else None,
        }
        for k, d in differences.items()
    }
    return {
        "id": campaign_id,
        "created_at": record["created_at"],
        "experiment_id": recipe["experiment_id"],
        "objective": recipe["objective"],
        "status": state,
        "replicates": recipe["replicates"],
        "base_seed": recipe["base_seed"],
        "md_steps": recipe.get("md_steps"),
        "complete_pairs": complete,
        "arms": [arm["arm"] for arm in recipe["arms"]],
        "pairs": results,
        "metrics": metrics,
        "units": "water O neighbors per exposed site",
        "cutoff_A": _CUTOFF_A,
        "difference_convention": "second minus first",
        "scope": (
            "Fixed substrate. Mean over saved MD frames including step zero; "
            "replicate SD is between paired initializations, not an equilibrium interval."
            if recipe["objective"] == "model"
            else "Initial rigid packing; geometric water-O neighbors per fluid-facing OH/F site. "
            "Distances do not establish an equilibrated hydration or adsorption profile."
        ),
    }
