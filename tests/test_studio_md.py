import hashlib
import json
from pathlib import Path
import tempfile
import unittest

from polychemmd.studio.iff_sources import HASHES, install, load, status
from polychemmd.studio.md import _checkpoint, recipe, restart_coefficients, run_md


class IffCommissioningTests(unittest.TestCase):
    def test_rejects_user_coordinates_and_nonphysical_controls(self):
        for cfg in (
            {"salt_M": 0.1},
            {"md_steps": False},
            {"md_steps": 0},
            {"md_steps": 2001},
            {"seed": 0},
            {"md_profile": "arbitrary"},
            {"md_sequence": "unassigned"},
            {"restart_job": "job-without-resolver"},
        ):
            with self.subTest(cfg=cfg), self.assertRaises(ValueError):
                run_md(cfg, Path("unused"))

    def test_data_integrity_before_cache_mutation(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td)
            payload = b"private test bytes"
            digest = hashlib.sha256(payload).hexdigest()
            install(path, "test.data", payload, digest)
            with self.assertRaises(ValueError):
                install(path, "test.data", b"wrong", digest)
            self.assertEqual((path / "test.data").read_bytes(), payload)
            self.assertEqual(status(path), dict.fromkeys(HASHES, False))
            with self.assertRaises(RuntimeError):
                load("iff-water", cache=path)
            (path / "iff-water.data").write_bytes(payload)
            with self.assertRaises(ValueError):
                load("iff-water", cache=path)

    def test_parameter_profiles_and_periodic_stack(self):
        opls = recipe("iff-opls", 1101, 11)
        charmm = recipe("iff-charmm", 1101, 11)
        self.assertIn("pair_modify mix geometric", opls)
        self.assertIn("special_bonds lj/coul 0 0 0.5", opls)
        self.assertIn("include fixed_iff_pairs.in", charmm)
        self.assertIn("special_bonds lj/coul 0 0 1", charmm)
        for script in (opls, charmm):
            self.assertIn("read_data interface.data\nkspace_style pppm", script)
            self.assertIn("fix slab solid setforce 0 0 0", script)
            self.assertIn("run 11\nwrite_restart state.restart", script)

    def test_restart_inherits_coordinates_velocities_and_thermostat(self):
        script = recipe("iff-charmm", 1101, 11, continuation=True)
        self.assertIn("read_restart parent.restart", script)
        self.assertIn("include restore-coeffs.in", script)
        self.assertIn("fix propagate mobile nvt", script)
        self.assertNotIn("velocity ", script)
        self.assertNotIn("minimize ", script)
        self.assertNotIn("reset_timestep ", script)
        self.assertIn(
            b"pair_coeff 1 1 1.2 3",
            restart_coefficients(b"My table\n\nPair Coeffs\n\n1 1.2 3\n"),
        )

    def test_modified_checkpoint_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            path = Path(d)
            (path / "engine").mkdir()
            original = b"original restart"
            script = b"source commands"
            digest = lambda x: hashlib.sha256(x).hexdigest()
            (path / "engine/input.lammps").write_bytes(script)
            (path / "engine/state.restart").write_bytes(original)
            (path / "md_summary.json").write_text(
                json.dumps({"input_script_sha256": digest(script)})
            )
            (path / "engine/run_record.json").write_text(
                json.dumps(
                    {
                        "status": "completed",
                        "engine_identity": {"version": "fixture"},
                        "artifact_sha256": {"state.restart": digest(original)},
                    }
                )
            )
            self.assertEqual(_checkpoint(path)[2], original)
            (path / "engine/state.restart").write_bytes(b"changed")
            with self.assertRaises(ValueError):
                _checkpoint(path)
