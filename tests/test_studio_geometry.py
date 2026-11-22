"""Chemical identities, surface periodicity, and limitations of the preparatory model."""

import hashlib
import math
import tempfile
import unittest
from collections import Counter
from importlib.util import find_spec
from pathlib import Path

from polychemmd.studio.geometry import build_geometry


class GeometryInputTests(unittest.TestCase):
    def test_bad_inputs_do_not_reach_optional_dependencies(self):
        invalid = [
            {"nx": True},
            {"gap_A": float("nan")},
            {"dp": 0},
            {"pattern": "other"},
            {"polymer": "PEO", "hydrolysis": 0.8},
            {"polymer": "bare", "chains": 1},
            {"polymer": "PEI", "side_dp": 3},
        ]
        for config in invalid:
            with self.subTest(config=config), tempfile.TemporaryDirectory() as folder:
                with self.assertRaises(ValueError):
                    build_geometry(config, Path(folder))


@unittest.skipUnless(find_spec("ase"), "requires the Studio ASE extra")
class CrystalGeometryTests(unittest.TestCase):
    def test_orthorhombic_bare_hydrated_slit_and_stoichiometry(self):
        cfg = {
            "geometry": "slit",
            "nx": 2,
            "ny": 3,
            "gap_A": 17,
            "polymer": "bare",
            "oh_top": 0.25,
            "oh_bottom": 1,
            "pattern": "patch",
        }
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            s = build_geometry(cfg, out)
            self.assertEqual(
                s["xyz_sha256"],
                hashlib.sha256((out / "structure.xyz").read_bytes()).hexdigest(),
            )
        counts = Counter(a["element"] for a in s["atoms"])
        nface = 2 * 2 * 3
        self.assertEqual(
            counts, {"Ti": 3 * nface * 2, "C": 2 * nface * 2, "O": 30, "H": 30, "F": 18}
        )
        self.assertAlmostEqual(s["accessible_z_A"][1] - s["accessible_z_A"][0], 17)
        self.assertGreater(s["center_spacing_A"], 17)
        self.assertEqual(s["pbc"], [True, True, False])
        lx, ly, lz = s["box_A"]
        self.assertAlmostEqual(s["accessible_volume_A3"], lx * ly * 17)
        for a in s["atoms"]:
            self.assertTrue(0 <= a["x"] < lx and 0 <= a["y"] < ly and 0 <= a["z"] < lz)
        # Check the physical O-H bond even for a wrapped atom at a box edge.
        lookup = {
            (a["sheet_index"], a["layer"], a["site"], a["element"]): a
            for a in s["atoms"]
        }
        for a in s["atoms"]:
            if a["element"] != "H":
                continue
            o = lookup[a["sheet_index"], a["layer"], a["site"], "O"]
            dx = a["x"] - o["x"]
            dy = a["y"] - o["y"]
            dz = a["z"] - o["z"]
            self.assertAlmostEqual(
                math.hypot(dx - lx * round(dx / lx), dy - ly * round(dy / ly), dz),
                0.94,
                places=6,
            )


@unittest.skipUnless(
    find_spec("ase") and find_spec("rdkit"), "requires the Studio chemistry extras"
)
class PolymerGeometryTests(unittest.TestCase):
    def test_chemical_formula_and_mapped_end_groups(self):
        from rdkit import Chem

        cases = [
            ("PEO", "C16H34O9"),
            ("PVA", "C17H36O8"),
            ("PEtOx", "C41H76N8O8"),
            ("PVP", "C49H76N8O8"),
            ("PEI", "C17H44N8"),
            ("PCE", "C35H56O20"),
        ]
        for polymer, formula in cases:
            with self.subTest(polymer=polymer), tempfile.TemporaryDirectory() as folder:
                out = Path(folder)
                s = build_geometry({"polymer": polymer}, out)
                p = s["polymer"]
                self.assertEqual(p["molecular_formula"], formula)
                self.assertEqual(p["formal_charge"], 0)
                if polymer == "PCE":
                    self.assertEqual(p["graft_residues_1based"], [1, 5])
                m = Chem.SDMolSupplier(str(out / "polymer.sdf"), removeHs=False)[0]
                atoms = s["templates"]["polymer"]
                self.assertEqual(m.GetNumAtoms(), len(atoms))
                self.assertEqual(len(s["polymer_bonds"]), m.GetNumBonds())
                for row, point in zip(atoms, m.GetConformer().GetPositions()):
                    for axis, value in zip(("x", "y", "z"), point):
                        self.assertAlmostEqual(row[axis], value, delta=5.1e-5)
                self.assertEqual(
                    [a.GetAtomMapNum() for a in m.GetAtoms()],
                    [a["atom_map"] for a in atoms],
                )
                self.assertEqual(
                    p["geometry_checks"]["stereochemistry"][
                        "specified_double_bonds_checked"
                    ],
                    0,
                )

    def test_pva_sequence_and_geometry_deterministic(self):
        cfg = {
            "polymer": "PVA",
            "dp": 8,
            "hydrolysis": 0.75,
            "acetate_pattern": "clustered",
            "seed": 11,
        }
        with tempfile.TemporaryDirectory() as folder:
            out = Path(folder)
            a = build_geometry(cfg, out)
            b = build_geometry(cfg, out)
        residues = a["polymer"]["acetate_residues_1based"]
        self.assertEqual(len(residues), 2)
        self.assertEqual(residues[1] - residues[0], 1)
        self.assertEqual(a["polymer"]["molecular_formula"], "C21H40O10")
        self.assertEqual(a["xyz_sha256"], b["xyz_sha256"])
        self.assertNotEqual(
            a["templates"]["polymer"][0]["atom_id"],
            a["templates"]["polymer"][1]["atom_id"],
        )
