import json
import sqlite3
import tempfile
import unittest
from concurrent.futures import ThreadPoolExecutor
from contextlib import closing
from importlib.util import find_spec

from polychemmd.studio.campaign import campaign_summary, submit_campaign
from polychemmd.studio.jobs import JobStore


class PairedCampaignTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.store = JobStore(self.tmp.name)

    def test_concurrent_idempotency_and_correct_arm_config(self):
        cfg = {
            "objective": "hydrolysis",
            "base_seed": 1234,
            "replicates": 3,
            "experiment_id": "test-A",
        }
        with ThreadPoolExecutor(2) as pool:
            a, b = pool.map(lambda _: submit_campaign(self.store, cfg), range(2))
        self.assertEqual(a["id"], b["id"])
        self.assertEqual(len(self.store.list()), 6)
        self.assertEqual(len({p["seed"] for p in a["pairs"]}), 3)
        self.assertEqual(a["pairs"][0]["seed"], 1234)
        self.assertEqual(a["complete_pairs"], 0)
        self.assertEqual(a["status"], "queued")
        self.assertIsNone(a["metrics"]["OH"]["mean_difference"])
        for pair in a["pairs"]:
            jobs = [self.store.get(j["id"]) for j in pair["jobs"]]
            for j in jobs:
                self.assertEqual(j["kind"], "packing")
                self.assertEqual(j["config"]["seed"], pair["seed"])
                self.assertEqual(j["config"]["hydrolysis"], 0.75)
            self.assertEqual(
                [j["config"]["acetate_pattern"] for j in jobs], ["random", "clustered"]
            )
        with self.assertRaises(ValueError):
            submit_campaign(self.store, {**cfg, "replicates": 1})
        # Personal normal-job idempotency keys cannot be mistaken for campaign members.
        self.store.submit("geometry", {}, a["id"])
        self.assertEqual(len(campaign_summary(self.store, a["id"])["pairs"]), 3)

    def test_atomic_failure_and_bad_parameters(self):
        for bad in (
            {"objective": "unregistered"},
            {"replicates": False},
            {"replicates": 13},
            {"objective": "model", "md_steps": 0},
            {"md_steps": 10},
            {"experiment_id": "../x"},
        ):
            with self.subTest(bad=bad), self.assertRaises(ValueError):
                submit_campaign(self.store, bad)
        with closing(self.store.connect()) as db:
            db.execute("""CREATE TRIGGER stop_second BEFORE INSERT ON jobs WHEN
                    (SELECT count(*) FROM jobs)>0 BEGIN SELECT RAISE(ABORT,'trigger rollback'); END""")
        with self.assertRaises(sqlite3.IntegrityError):
            submit_campaign(self.store, {"replicates": 2})
        self.assertEqual(self.store.list(), [])
        with self.assertRaises(KeyError):
            campaign_summary(self.store, "cmp_" + "c" * 32)

    def test_md_differences_uncertainty_and_partial_failure(self):
        new = submit_campaign(
            self.store, {"objective": "model", "replicates": 3, "md_steps": 2}
        )
        self.assertEqual(new["arms"], ["iff-opls", "iff-charmm"])
        for pair in new["pairs"]:
            for i, j in enumerate(pair["jobs"]):
                job_id, owner, _, cfg = self.store.claim(j["id"])
                self.assertEqual(cfg["md_steps"], 2)
                if pair["replicate"] == 3 and i == 1:
                    self.store.finish(job_id, owner, error="Test fixture failure")
                else:
                    self.store.finish(
                        job_id,
                        owner,
                        result={
                            "water_contacts": {
                                "OH": 1 + i * pair["replicate"],
                                "F": 2 - i,
                            }
                        },
                    )
        s = campaign_summary(self.store, new["id"])
        self.assertEqual(s["status"], "partial_failed")
        self.assertEqual(s["complete_pairs"], 2)
        self.assertEqual(
            s["metrics"]["F"], {"n": 2, "mean_difference": -1, "sd_difference": 0.0}
        )
        self.assertEqual(s["metrics"]["OH"]["mean_difference"], 1.5)
        self.assertAlmostEqual(s["metrics"]["OH"]["sd_difference"], 2**-0.5)
        self.assertIsNone(s["pairs"][-1]["difference"])

    @unittest.skipUnless(find_spec("numpy"), "requires NumPy")
    def test_only_fluid_facing_sites_and_lateral_periodicity(self):
        new = submit_campaign(self.store, {"objective": "termination", "replicates": 1})

        def atom(element, role, x, z, **extra):
            return {"element": element, "role": role, "x": x, "y": 5.0, "z": z, **extra}

        for index, j in enumerate(new["pairs"][0]["jobs"]):
            job_id, owner, _, _ = self.store.claim(j["id"])
            area = self.store.root / job_id
            area.mkdir()
            # Water wraps across x, remains 8 Å from the bottom-face control.
            water = atom("O", "water", 19.5 if index else 10.0, 12.0)
            atoms = [
                atom("O", "slab", 0.5, 10.0, layer="top_OH", sheet_index=0),
                atom("O", "slab", 0.5, 4.0, layer="bottom_OH", sheet_index=0),
                atom("F", "slab", 12.0, 10.0, layer="top_F", sheet_index=0),
                water,
            ]
            (area / "packed.json").write_text(
                json.dumps({"box_A": [20, 20, 30], "atoms": atoms})
            )
            self.store.finish(job_id, owner, result={}, output=job_id)
        s = campaign_summary(self.store, new["id"])
        self.assertEqual(s["scope"].split(";")[0], "Initial rigid packing")
        self.assertEqual(s["pairs"][0]["difference"], {"OH": 1.0, "F": -1.0})
        self.assertIsNone(s["metrics"]["F"]["sd_difference"])
