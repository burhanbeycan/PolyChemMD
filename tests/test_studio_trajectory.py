"""Analytic controls in Å. The synthetic coordinates define no material properties."""

import csv
import gzip
import json
from pathlib import Path
import unittest
from importlib.util import find_spec
from tempfile import TemporaryDirectory

from polychemmd.studio.trajectory import analyze_trajectory, frames, run_trajectory

HAS_NUMERICS = find_spec("numpy") is not None and find_spec("scipy") is not None
GROUPS = {
    "reference_plane": [1, 2],
    "comparison_plane": [7],
    "surface_sites": [3],
    "water_oxygen": [4, 5],
    "polymer_heavy": [6],
    "plane_restraint": "fixed",
}
ATOMS = [
    (1, 0.5, 0.5, 0.95),
    (2, 0.5, 0.5, 0.05),
    (3, 0.98, 0.20, 0.98),
    (4, 0.01, 0.20, 0.01),
    (5, 0.5, 0.5, 0.51),
    (6, 0.0, 0.20, 0.01),
    (7, 0.5, 0.5, 0.4),
]


def dump_frame(step=0, atoms=ATOMS, *, tilt=False, zbc="pp", lz=10, scaled=True):
    boxes = (
        f"xy xz yz pp pp {zbc}\n0 15 4\n0 10 1\n0 {lz} 0\n"
        if tilt
        else (f"pp pp {zbc}\n0 10\n0 10\n0 {lz}\n")
    )
    return (
        f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{len(atoms)}\n"
        f"ITEM: BOX BOUNDS {boxes}ITEM: ATOMS id "
        + ("xs ys zs\n" if scaled else "x y z\n")
        + "".join(" ".join(str(v) for v in row) + "\n" for row in atoms)
    )


@unittest.skipUnless(HAS_NUMERICS, "NumPy and SciPy are optional Studio dependencies")
class TestStudioTrajectory(unittest.TestCase):
    def setUp(self):
        self.temp = TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.dump = self.root / "engine" / "trajectory.dump"
        self.dump.parent.mkdir()
        self.output = self.root / "post"
        self.metadata = self.root / "trajectory-groups.json"
        self.metadata.write_text(json.dumps(GROUPS))

    def run_case(self, text, groups=GROUPS, **kw):
        self.dump.write_text(text)
        return analyze_trajectory(self.dump, groups, self.output, **kw)

    def csv(self, name):
        with (self.output / name).open() as handle:
            return list(csv.DictReader(handle))

    def test_orthorhombic_contacts_order_periodic_plane_and_mass_conservation(self):
        records = dump_frame() + dump_frame(10, list(reversed(ATOMS)))
        res = self.run_case(records, bin_width_A=1, contact_cutoff_A=0.6)
        self.assertEqual(res["n_frames"], 2)
        self.assertAlmostEqual(res["mean_plane_separation_nm"], 0.4)
        self.assertEqual(res["mean_water_O_contact_fraction"], 0.5)
        self.assertEqual(res["mean_polymer_heavy_contact_fraction"], 1)
        densities = self.csv("density.csv")
        self.assertEqual(len(densities), 10)
        self.assertAlmostEqual(
            sum(
                float(r["water_O_number_density_A-3"]) * float(r["mean_bin_volume_A3"])
                for r in densities
            ),
            2,
        )
        self.assertTrue(all(float(r["frame_sd_A-3"]) == 0 for r in densities))
        self.assertAlmostEqual(
            float(self.csv("contacts.csv")[0]["complementary_separation_A"]), 6
        )

    def test_restricted_triclinic_periodic_y_contact(self):
        atoms = list(ATOMS)
        atoms[2] = (3, 0.2, 0.99, 0.5)
        atoms[3] = (4, 0.2, 0.01, 0.5)
        self.run_case(dump_frame(atoms=atoms, tilt=True), contact_cutoff_A=0.23)
        selected = self.csv("contacts.csv")[0]
        # sqrt((.02 * 4)^2 + (.02 * 10)^2) = 0.2154 Å
        self.assertEqual(int(selected["water_O_in_contact"]), 1)
        self.run_case(dump_frame(atoms=atoms, tilt=True), contact_cutoff_A=0.21)
        self.assertEqual(int(self.csv("contacts.csv")[0]["water_O_in_contact"]), 0)

    def test_cartesian_columns_nonsingular_box(self):
        a = [(i, x * 10, y * 10, z * 10) for i, x, y, z in ATOMS]
        res = self.run_case(dump_frame(atoms=a, scaled=False), contact_cutoff_A=0.6)
        self.assertEqual(res["mean_water_O_contact_fraction"], 0.5)

    def test_per_frame_volume_on_changing_box(self):
        self.run_case(dump_frame() + dump_frame(10, lz=20), bin_width_A=1)
        densities = self.csv("density.csv")
        nonzero = [r for r in densities if float(r["water_O_number_density_A-3"]) > 0]
        self.assertEqual(len(nonzero), 2)
        self.assertAlmostEqual(float(nonzero[0]["water_O_number_density_A-3"]), 0.0075)

    def test_empty_polymer_and_absent_comparison(self):
        self.dump.write_text(dump_frame())
        self.metadata.write_text(
            json.dumps(
                {
                    **GROUPS,
                    "reference_plane": [1],
                    "comparison_plane": [],
                    "polymer_heavy": [],
                }
            )
        )
        gzip_path = self.dump.with_suffix(".dump.gz")
        with gzip.open(gzip_path, "wb") as f:
            f.write(self.dump.read_bytes())
        self.dump.unlink()
        result = run_trajectory(
            {"source_job": "e01"}, self.output, source_dir=self.root
        )
        self.assertIsNone(result["mean_polymer_heavy_contact_fraction"])
        self.assertIsNone(result["plane_distance_observable"])
        self.assertEqual(len(result["source_dump_sha256"]), 64)
        self.assertEqual(result["source_job"], "e01")

    def test_nonperiodic_cross_boundary_contact_excluded(self):
        res = self.run_case(dump_frame(zbc="ff"), contact_cutoff_A=0.6)
        self.assertEqual(res["mean_water_O_contact_fraction"], 0)
        self.assertEqual(res["density_reference"], "nonperiodic lower cell boundary")

    def test_invalid_data_rejected(self):
        repeated = list(ATOMS)
        repeated[1] = repeated[0]
        invalid = [
            dump_frame(atoms=repeated),
            dump_frame() + dump_frame(),
            dump_frame().replace("0.95", "nan"),
            dump_frame()[:-20],
            dump_frame() + dump_frame(1, list(ATOMS[:-1])),
        ]
        for text in invalid:
            with self.subTest(text=text[:20]), self.assertRaises(ValueError):
                self.run_case(text)
        with self.assertRaises(ValueError):
            self.run_case(dump_frame(), groups={**GROUPS, "surface_sites": [4]})
        with self.assertRaises(ValueError):
            self.run_case(dump_frame(), contact_cutoff_A=True)

    @unittest.skipUnless(
        find_spec("MDAnalysis"), "MDAnalysis is an optional cross-check"
    )
    def test_minimum_image_against_mdanalysis(self):
        import numpy as np
        from MDAnalysis.lib.distances import distance_array
        from MDAnalysis.lib.mdamath import triclinic_box
        from polychemmd.studio.trajectory import _contact_mask

        rng = np.random.default_rng(37)
        atoms = [(i + 1, *row) for i, row in enumerate(rng.random((40, 3)))]
        self.dump.write_text(dump_frame(atoms=atoms, tilt=True))
        f = next(frames(self.dump))
        xyz = f.fractional @ f.cell
        truth = (
            distance_array(xyz[10:], xyz[:10], box=triclinic_box(*f.cell)).min(axis=1)
            <= 2.4
        )
        np.testing.assert_array_equal(
            _contact_mask(f, list(range(10)), list(range(10, 40)), 2.4), truth
        )
