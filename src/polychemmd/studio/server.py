"""Loopback-only HTTP application. Long-running jobs are owned by the SQLite worker."""

from contextlib import asynccontextmanager
from pathlib import Path
import threading
from urllib.parse import urlsplit

from fastapi import FastAPI, HTTPException, Request
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles

from .jobs import JobStore, work
from .operations import status
from .planner import suggest

LOCAL = {"localhost", "127.0.0.1", "::1"}


def app_factory(root, *, start_worker=True):
    store = JobStore(root)

    @asynccontextmanager
    async def lifespan(_app):
        stop = threading.Event()
        runner = None
        if start_worker:
            runner = threading.Thread(
                target=work, args=(store, stop), name="studio-worker", daemon=True
            )
            runner.start()
        try:
            yield
        finally:
            stop.set()
            if runner:
                runner.join(timeout=2)

    app = FastAPI(
        title="PolyChemMD • MXene–Polimer",
        lifespan=lifespan,
        docs_url=None,
        redoc_url=None,
    )
    app.state.store = store

    @app.middleware("http")
    async def loopback(request, call_next):
        host = request.headers.get("host", "")
        try:
            allowed = urlsplit("//" + host).hostname in LOCAL
        except ValueError:
            allowed = False
        if not allowed:
            return JSONResponse(
                {"detail": "Use the local Studio address."}, status_code=403
            )
        if request.method == "POST":
            origin = request.headers.get("origin")
            if origin and origin != f"{request.url.scheme}://{host}":
                return JSONResponse(
                    {"detail": "Unknown request origin."}, status_code=403
                )
            if (
                request.headers.get("content-type", "").split(";")[0].strip()
                != "application/json"
            ):
                return JSONResponse({"detail": "Send a JSON request."}, status_code=415)
        return await call_next(request)

    async def body(request):
        raw = bytearray()
        async for chunk in request.stream():
            raw.extend(chunk)
            if len(raw) > 120000:
                raise HTTPException(413, "Request too large.")
        import json

        try:
            item = json.loads(raw)
            if not isinstance(item, dict):
                raise ValueError()
            return item
        except (ValueError, UnicodeError) as e:
            raise HTTPException(422, "Expected a JSON object.") from e

    @app.get("/api/status")
    def service_status():
        return status()

    @app.post("/api/campaigns", status_code=201)
    async def add_campaign(request: Request):
        from .campaign import submit_campaign

        data = await body(request)
        try:
            return submit_campaign(store, data)
        except (ValueError, TypeError) as e:
            raise HTTPException(422, str(e)) from e

    @app.get("/api/campaigns/{campaign_id}")
    def inspect_campaign(campaign_id: str):
        from .campaign import campaign_summary

        try:
            return campaign_summary(store, campaign_id)
        except KeyError as e:
            raise HTTPException(404, "Unknown campaign.") from e
        except ValueError as e:
            raise HTTPException(422, str(e)) from e

    @app.get("/api/profiles")
    def registered_profiles():
        from .register_topology import profiles

        return {"profiles": profiles()}

    @app.get("/api/atlas")
    def atlas_data():
        from .experiment import atlas

        return atlas()

    @app.get("/api/jobs")
    def jobs(limit: int = 100):
        return {"jobs": store.list(limit)}

    @app.post("/api/jobs", status_code=201)
    async def create_job(request: Request):
        data = await body(request)
        if set(data) - {"kind", "config"}:
            raise HTTPException(422, "Only kind and config are accepted.")
        try:
            return store.submit(
                data.get("kind"),
                data.get("config"),
                request.headers.get("idempotency-key"),
            )
        except (ValueError, TypeError, OverflowError) as e:
            raise HTTPException(422, str(e)) from e

    @app.post("/api/plan")
    async def plan(request: Request):
        data = await body(request)
        try:
            from starlette.concurrency import run_in_threadpool

            return await run_in_threadpool(
                suggest,
                data.get("mode", "rules"),
                data.get("model"),
                data.get("objective", "termination"),
            )
        except (OSError, ValueError, KeyError, TypeError) as e:
            raise HTTPException(422, str(e)) from e

    def get_job(job_id):
        try:
            return store.get(job_id)
        except KeyError as e:
            raise HTTPException(404, "Unknown job.") from e

    @app.get("/api/jobs/{job_id}")
    def job(job_id: str):
        return get_job(job_id)

    @app.get("/api/jobs/{job_id}/structure")
    def structure(job_id: str):
        j = get_job(job_id)
        for f in ("packed.json", "structure.json"):
            if f in j["files"]:
                return FileResponse(
                    store.output_dir(job_id) / f, media_type="application/json"
                )
        raise HTTPException(404, "This job has no molecular coordinates.")

    @app.get("/api/jobs/{job_id}/files/{name:path}")
    def artifact(job_id: str, name: str):
        j = get_job(job_id)
        if name not in j["files"]:
            raise HTTPException(404, "Unknown job artifact.")
        return FileResponse(store.output_dir(job_id) / name, filename=name)

    app.mount(
        "/",
        StaticFiles(directory=Path(__file__).with_name("web"), html=True),
        name="studio",
    )
    return app
