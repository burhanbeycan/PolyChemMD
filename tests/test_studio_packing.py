"""Geometry, units and native-engine boundaries of the optional studio packers."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

from polychemmd.studio.packing import _audit, pack_system, packmol_path


HAS_NUMERICS = all(
    importlib.util.find_spec(pkg) is not None for pkg in ("scipy", "numpy")
)
ATOM = {
    "element": "Ti",
    "x": 0.2,
    "y": 0.2,
    "z": 2.0,
    "role": "slab",
    "molecule_id": "sheet:1",
}


@unittest.skipUnless(HAS_NUMERICS, "Studio packing needs NumPy and SciPy")
class TestStudioPacking(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.out = Path(self.temp.name)
        self.shape = {
            "units": "angstrom",
            "box_A": [20.0, 20.0, 22.0],
            "accessible_z_A": [2.0, 22.0],
            "geometry": "single",
            "templates": {"slab": [ATOM], "polymer": []},
            "polymer_bonds": [],
        }
        self.config = {"chains": 0, "water_count": 16, "salt_M": 0.25, "seed": 17}
        (self.out / "structure.json").write_text(json.dumps(self.shape))

    def check_backend(self, backend):
        result = pack_system({**self.config, "packing_engine": backend}, self.out)
        self.assertEqual(result["water_count"], 16)
        self.assertEqual(result["ion_pairs"], 1)
        self.assertAlmostEqual(
            result["salt_M_realized"], 1 / (6.02214076e23 * 8000e-27)
        )
        self.assertTrue(result["audit"]["fixed_surface_invariant"])
        self.assertGreaterEqual(
            result["audit"]["minimum_mobile_to_other_A"], 1.8 - 0.002
        )
        self.assertEqual(result["water_count_basis"], "explicit")
        geom = json.loads((self.out / "packed.json").read_text())
        self.assertEqual(len(geom["atoms"]), 1 + 3 * 16 + 2)
        self.assertEqual(len(geom["bonds"]), 2 * 16)
        self.assertEqual(geom["atoms"][0]["x"], ATOM["x"])
        self.assertEqual(len(geom["atoms"]), len({a["id"] for a in geom["atoms"]}))
        return result

    def test_builtin_reproduces_coordinates_and_mass(self):
        a = self.check_backend("builtin")
        b = self.check_backend("builtin")
        self.assertEqual(a["packed_xyz_sha256"], b["packed_xyz_sha256"])
        self.assertAlmostEqual(
            a["water_density_realized_g_cm3"],
            16 * 18.01528 / (6.02214076e23 * 8000e-24),
        )

    @unittest.skipUnless(
        importlib.util.find_spec("molpack"), "MolPack Python bindings unavailable"
    )
    def test_molpack_runs(self):
        a = self.check_backend("molpack")
        self.assertEqual(a["engine_info"]["package"], "molcrafts-molpack")

    @unittest.skipUnless(packmol_path(), "Canonical Packmol executable unavailable")
    def test_packmol_runs(self):
        self.check_backend("packmol")
        self.assertTrue((self.out / "packing_inputs" / "packmol.log").is_file())

    def test_invalid_counts_and_engines(self):
        for change in (
            {"water_count": 1.1},
            {"packing_engine": "other"},
            {"salt_M": float("nan")},
            {"chains": True},
            {"min_distance_A": float("inf")},
        ):
            with self.subTest(change=change), self.assertRaises(ValueError):
                pack_system({**self.config, **change}, self.out)

    def test_known_periodic_cross_boundary_collision_fails(self):
        import numpy as np

        atoms = [
            ATOM,
            dict(
                ATOM, element="Na", x=19.8, z=2.2, molecule_id="sodium:1", role="sodium"
            ),
        ]
        coords = np.array([[a[k] for k in ("x", "y", "z")] for a in atoms])
        with self.assertRaisesRegex(ValueError, "Minimum XY-image"):
            _audit(
                atoms,
                coords,
                [("slab", 0, 1), ("sodium", 1, 1)],
                [("slab", atoms[:1], [], 1), ("sodium", atoms[1:], [], 1)],
                np.array([20.0, 20.0, 22.0]),
                np.array([2.0, 22.0]),
                1.8,
            )

    def test_dry_empty_cell_emits_no_invented_distance(self):
        r = pack_system({**self.config, "water_count": 0, "salt_M": 0}, self.out)
        self.assertIsNone(r["audit"]["minimum_mobile_to_other_A"])
        self.assertEqual(r["engine_info"]["package"], "none")


if __name__ == "__main__":
    unittest.main()
