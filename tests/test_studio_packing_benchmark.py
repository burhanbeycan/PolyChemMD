"""Scientific regression checks for confined rigid poses, parity and paired benchmark records."""

import importlib.util
import json
from pathlib import Path
from unittest.mock import patch
import tempfile
import unittest

from polychemmd.studio.packing import _audit, _check_rigid_fit, _slit_pose, pack_system
from polychemmd.studio.packing_benchmark import (
    _list,
    _paired_comparisons,
    _summary,
    run_benchmark,
)

HAS_NUMERICS = all(importlib.util.find_spec(x) for x in ("numpy", "scipy"))
HAS_STUDIO = HAS_NUMERICS and all(importlib.util.find_spec(x) for x in ("ase", "rdkit"))


def atom(v):
    return dict(element="C", x=v[0], y=v[1], z=v[2], role="polymer")


@unittest.skipUnless(HAS_NUMERICS, "Needs optional numerical dependencies")
class RigidProposalTests(unittest.TestCase):
    def test_feasible_translation_preserves_handedness_and_slit(self):
        import numpy as np
        from scipy.spatial.distance import pdist

        pts = np.array([[-4.0, 0.2, 0], [4, 0.1, 0], [1.1, 2.0, 0.1], [-1, -1, 0.5]])
        centered = pts - pts.mean(axis=0)
        _, _, vh = np.linalg.svd(centered, full_matrices=False)
        axes = vh.T.copy()
        axes[:, -1] *= np.sign(np.linalg.det(axes))
        rng = np.random.default_rng(97)
        low = np.array([0, 0, 4])
        high = np.array([20, 20, 5])
        placed = _slit_pose(centered, axes, rng, low, high, 0)
        self.assertIsNotNone(placed)
        self.assertTrue(np.all(placed >= low) and np.all(placed <= high))
        self.assertTrue(np.allclose(pdist(pts), pdist(placed)))
        signed = lambda a: np.linalg.det(a[1:] - a[0])
        self.assertAlmostEqual(signed(pts), signed(placed), places=9)

    def test_fit_bound_rejects_impossible_diagonal_and_width(self):
        import numpy as np
        from itertools import product

        slab = ("slab", [atom((0, 0, 0))], [], 1)
        narrow = np.array([4.0, 4.0, 4.0])
        bounds = np.array([2.0, 3.0])
        for pts, message in (
            ([(0, 0, 0), (8, 0, 0)], "diameter"),
            (list(product((-1, 1), repeat=3)), "width"),
        ):
            with (
                self.subTest(message=message),
                self.assertRaisesRegex(ValueError, message),
            ):
                _check_rigid_fit(
                    [slab, ("polymer", [atom(x) for x in pts], [], 1)], narrow, bounds
                )

    def test_independent_audit_rejects_mirrored_tetrahedron(self):
        import numpy as np

        pts = np.array([[0, 0, 0], [1, 0, 0], [0, 2, 0], [0, 0, 3]], dtype=float)
        atoms = [atom((1, 1, 1))] + [atom(x) for x in pts]
        mirror = np.vstack(([1, 1, 1], pts * [-1, 1, 1] + 7))
        specs = [("slab", atoms[:1], [], 1), ("polymer", atoms[1:], [], 1)]
        with self.assertRaisesRegex(ValueError, "handedness"):
            _audit(
                atoms,
                mirror,
                [("slab", 0, 1), ("polymer", 1, 4)],
                specs,
                np.array([20.0, 20.0, 20.0]),
                np.array([1.0, 19.0]),
                1.8,
            )

    def test_confined_restarts_have_distinct_seeds_and_account_for_rejected_work(self):
        from polychemmd.studio.packing import _InsertionExhausted, _run_builtin

        calls = []

        def simulate(*args, **kwargs):
            calls.append(args[4])
            if len(calls) == 1:
                raise _InsertionExhausted("water", 2, 3, 777)
            return "coordinates", {"placement_attempts": 12}

        with patch("polychemmd.studio.packing._builtin", side_effect=simulate):
            _, info = _run_builtin(
                None, None, None, 1.8, 17, slit=True, matched_constraints=True
            )
        self.assertNotEqual(calls[0], calls[1])
        self.assertEqual(info["rejected_starts"][0]["seed"], 17)
        self.assertEqual(info["successful_placement_seed"], calls[1])
        self.assertEqual(info["total_placement_attempts"], 789)

    def test_seed_override_controls_only_placement(self):
        shape = dict(
            box_A=[14, 14, 15],
            accessible_z_A=[2, 15],
            geometry="single",
            templates={"slab": [atom((1, 1, 1))], "polymer": []},
            polymer_bonds=[],
        )
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            (out / "structure.json").write_text(json.dumps(shape))
            config = dict(
                packing_engine="slit",
                water_count=10,
                salt_M=0,
                chains=0,
                packing_seed=101,
            )
            first = pack_system({**config, "seed": 99}, out)
            second = pack_system({**config, "seed": 7}, out)
            self.assertEqual(first["packed_xyz_sha256"], second["packed_xyz_sha256"])
            third = pack_system({**config, "packing_seed": 102}, out)
            self.assertNotEqual(first["packed_xyz_sha256"], third["packed_xyz_sha256"])
            self.assertEqual(
                first["geometry_input_sha256"], third["geometry_input_sha256"]
            )


class PackingBenchmarkTests(unittest.TestCase):
    def test_validation_and_failed_pairs(self):
        for values in ([True], [], [3, 3], [2.5], list(range(1, 10))):
            with self.subTest(values=values), self.assertRaises(ValueError):
                _list({"seeds": values}, "seeds", (17,))
        rows = [
            dict(
                case="c",
                engine=engine,
                seed=seed,
                status=status,
                minimum_mobile_to_other_A=1.9,
                wall_seconds=clock,
                force_status="not_requested",
            )
            for engine, seed, status, clock in [
                ("builtin", 1, "success", 2),
                ("slit", 1, "success", 1),
                ("builtin", 2, "failed", 3),
                ("slit", 2, "success", 2),
                ("builtin", 3, "unavailable", 0),
                ("slit", 3, "timeout", 2),
            ]
        ]
        comp = _paired_comparisons(rows, ["c"], ["builtin", "slit"])[0]
        self.assertEqual(comp["joint_success_count"], 1)
        self.assertEqual(comp["engine_only_success_count"], 1)
        self.assertEqual(comp["neither_success_count"], 1)
        self.assertEqual(comp["median_paired_seconds_engine_minus_reference"], -1)
        self.assertEqual(_summary(rows, ["c"], ["builtin"])[0]["unavailable_count"], 1)

    @unittest.skipUnless(HAS_STUDIO, "Full geometry and packing optional dependencies")
    def test_live_paired_subprocess_writes_measurements(self):
        with tempfile.TemporaryDirectory() as tmp:
            out = Path(tmp)
            result = run_benchmark(
                {"cases": ["aqueous"], "engines": ["builtin", "slit"], "seeds": [7]},
                out,
            )
            a, b = result["rows"]
            self.assertEqual([a["status"], b["status"]], ["success", "success"])
            self.assertEqual(a["geometry_sha256"], b["geometry_sha256"])
            self.assertEqual(a["atom_count"], b["atom_count"])
            self.assertTrue(
                a["proper_rigid_transforms"] and b["proper_rigid_transforms"]
            )
            self.assertEqual(result["paired_comparisons"][0]["joint_success_count"], 1)
            self.assertTrue((out / "benchmark_trials.csv").is_file())
