import tempfile
import unittest
from contextlib import closing
from importlib.util import find_spec
from polychemmd.studio.jobs import JobStore, execute_next
from polychemmd.studio.planner import suggest


class TestStudioJobs(unittest.TestCase):
    def setUp(self):
        t = tempfile.TemporaryDirectory()
        self.addCleanup(t.cleanup)
        self.store = JobStore(t.name)

    def test_lease_excludes_old_worker(self):
        j = self.store.submit(
            "experiment", {"assay": "xrd", "reference": "xing-pei"}, "example"
        )
        self.assertEqual(
            j["id"], self.store.submit("experiment", j["config"], "example")["id"]
        )
        with self.assertRaises(ValueError):
            self.store.submit("md", {}, "example")
        _, worker, _, _ = self.store.claim()
        self.assertIsNone(self.store.claim())
        with closing(self.store.connect()) as db:
            db.execute("UPDATE jobs SET expires=0")
        self.assertTrue(execute_next(self.store))
        self.assertFalse(self.store.finish(j["id"], worker, error="stale"))
        out = self.store.get(j["id"])
        self.assertEqual(out["status"], "succeeded")
        self.assertTrue((self.store.output_dir(j["id"]) / "manifest.json").is_file())
        self.assertNotIn("result", self.store.list()[0])

    def test_invalid_and_missing_inputs(self):
        for kind, config in (("unknown", {}), ("experiment", {"x": float("nan")})):
            with self.assertRaises(ValueError):
                self.store.submit(kind, config)
        job = self.store.submit("experiment", {"assay": "xrd", "rows": []})
        self.assertTrue(execute_next(self.store))
        self.assertEqual(self.store.get(job["id"])["status"], "failed")
        self.assertEqual(len(suggest(objective="hydrolysis")["jobs"]), 2)

    @unittest.skipUnless(
        find_spec("fastapi") and find_spec("httpx"), "Optional HTTP stack absent"
    )
    def test_http_workflow(self):
        from fastapi.testclient import TestClient
        from polychemmd.studio.server import app_factory

        with TestClient(
            app_factory(self.store.root, start_worker=False),
            base_url="http://127.0.0.1",
        ) as app:
            self.assertEqual(app.get("/").status_code, 200)
            recipe = {
                "kind": "experiment",
                "config": {"reference": "xing-pei", "assay": "xrd"},
            }
            self.assertEqual(
                app.post(
                    "/api/jobs",
                    json=recipe,
                    headers={"origin": "https://elsewhere.example"},
                ).status_code,
                403,
            )
            j = app.post("/api/jobs", json=recipe).json()
            self.assertTrue(execute_next(self.store))
            self.assertEqual(
                app.get("/api/jobs/" + j["id"]).json()["status"], "succeeded"
            )
            self.assertEqual(
                app.get("/api/jobs/" + j["id"] + "/files/experiment.csv").status_code,
                200,
            )
            self.assertEqual(
                app.get("/api/jobs/" + j["id"] + "/files/unknown.txt").status_code, 404
            )
            self.assertEqual(
                app.get("/api/jobs/" + j["id"] + "/structure").status_code, 404
            )
            self.assertGreater(len(app.get("/api/atlas").json()["records"]), 10)
