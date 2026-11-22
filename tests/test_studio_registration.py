import copy
import hashlib
import json
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from polychemmd.studio.register_topology import register, profiles, run_registered
from polychemmd.studio.topology_file import validate

# Four Ar atoms (one immobilised). This is a numerical LJ fixture and is not a MXene model.
COORDS = [(5.0, 5.0, 6.0), (9.0, 5.0, 6.0), (9.0, 9.0, 6.0), (5.0, 9.0, 6.0)]


def make(root):
    packed = {
        "pbc": [True, True, False],
        "units": "angstrom",
        "box_A": [16, 16, 16],
        "bonds": [],
        "atoms": [
            dict(
                id=i,
                element="Ar",
                x=x,
                y=y,
                z=z,
                role="slab" if i == 0 else "probe",
                molecule_id=i + 1,
            )
            for i, (x, y, z) in enumerate(COORDS)
        ],
    }
    packed_path = root / "packed.json"
    packed_path.write_text(json.dumps(packed))
    params = dict(
        schema="studio-full-slab-v1",
        label="Argon LJ verification fixture",
        parameter_source="Model LJ coefficients for an analytical regression",
        applicability="Static four-atom LJ test, no MXene chemistry",
        water_model="none",
        water_constraints="no_water",
        type_elements={"1": "Ar"},
        bond_orders={},
        mixing="arithmetic",
        lj14=0.5,
        coul14=0.5,
        cutoff_A=6,
        dt_fs=0.25,
        packed_sha256=hashlib.sha256(packed_path.read_bytes()).hexdigest(),
    )
    profile_path = root / "profile.json"
    profile_path.write_text(json.dumps(params))
    data = (
        "LJ test fixture\n\n4 atoms\n1 atom types\n0 16 xlo xhi\n0 16 ylo yhi\n0 16 zlo zhi\n\nMasses\n\n1 39.948\n\nPair Coeffs\n\n1 0.238 3.405\n\nAtoms # full\n\n"
        + "\n".join(f"{i} {i} 1 0 {x} {y} {z}" for i, (x, y, z) in enumerate(COORDS, 1))
        + "\n"
    )
    data_path = root / "typed.data"
    data_path.write_text(data)
    return packed_path, data_path, profile_path


class RegisteredTopologyTests(unittest.TestCase):
    def test_register_reuse_and_reject_drift(self):
        with (
            TemporaryDirectory() as d,
            patch.dict(os.environ, {"POLYCHEMMD_TOPOLOGY_REGISTRY": d + "/registry"}),
        ):
            root = Path(d)
            p, t, m = make(root)
            model = register(p, t, m)
            self.assertEqual(register(p, t, m)["id"], model["id"])
            self.assertEqual(profiles()[0]["packed_sha256"], model["packed_sha256"])
            p.write_text(p.read_text() + " ")
            with self.assertRaisesRegex(ValueError, "checksum"):
                run_registered({"topology_id": model["id"]}, root / "run", root)

    def test_reject_incompatible_topology(self):
        with TemporaryDirectory() as d:
            root = Path(d)
            p, t, m = make(root)
            template = json.loads(p.read_text())
            data = t.read_text()
            manifest = json.loads(m.read_text())
            for reason, mod in [
                ("Element", lambda v: v["atoms"][0].update(element="Ti")),
                ("Coordinates", lambda v: v["atoms"][0].update(x=8)),
                ("Bond", lambda v: v.update(bonds=[[0, 1, 1]])),
                ("periodic", lambda v: v.update(pbc=[True, True, True])),
            ]:
                altered = copy.deepcopy(template)
                mod(altered)
                with (
                    self.subTest(reason=reason),
                    self.assertRaisesRegex(ValueError, reason),
                ):
                    validate(altered, data, manifest)
            with self.assertRaisesRegex(ValueError, "neutral"):
                validate(template, data.replace("1 1 1 0 5.0", "1 1 1 1 5.0"), manifest)
            with self.assertRaises(ValueError):
                validate(template, data + "\nrun 100\n", manifest)

    @unittest.skipUnless(
        os.getenv("POLYCHEMMD_TEST_REGISTERED_LAMMPS") == "1",
        "optional actual LAMMPS process",
    )
    def test_real_static_lj_and_bounded_md(self):
        with (
            TemporaryDirectory() as d,
            patch.dict(os.environ, {"POLYCHEMMD_TOPOLOGY_REGISTRY": d + "/registry"}),
        ):
            root = Path(d)
            model = register(*make(root))
            cfg = {"topology_id": model["id"]}
            results = [
                run_registered(dict(cfg, steps=n), root / f"run{n}", root)
                for n in (0, 10)
            ]
            lj = lambda r: 4 * 0.238 * ((3.405 / r) ** 12 - (3.405 / r) ** 6)
            reference = 4 * lj(4) + 2 * lj(math.sqrt(32))
            self.assertAlmostEqual(results[0]["initial_kcal_mol"], reference, places=8)
            dump = (
                (root / "run0" / "registered-trajectory.dump").read_text().splitlines()
            )
            fx = float(dump[9].split()[5])
            fx_expected = sum(
                24 * 0.238 * (2 * (3.405 / r) ** 12 - (3.405 / r) ** 6) * dx / r**2
                for dx, r in [(-4, 4), (-4, math.sqrt(32))]
            )
            self.assertAlmostEqual(fx, fx_expected, places=8)
            self.assertAlmostEqual(results[1]["initial_kcal_mol"], reference, places=8)
            self.assertEqual(results[1]["steps"], 10)
            self.assertTrue(math.isfinite(results[1]["final_kcal_mol"]))

    @unittest.skipUnless(
        os.getenv("POLYCHEMMD_TEST_REGISTERED_LAMMPS") == "1",
        "optional actual LAMMPS process",
    )
    def test_real_charged_flexible_water_fixture(self):
        with (
            TemporaryDirectory() as d,
            patch.dict(os.environ, {"POLYCHEMMD_TOPOLOGY_REGISTRY": d + "/registry"}),
        ):
            root = Path(d)
            p, t, m = make(root)
            packed = json.loads(p.read_text())
            manifest = json.loads(m.read_text())
            packed["atoms"][1:] = [
                dict(
                    id=i, element=e, x=x, y=y, z=6, role="water", molecule_id="water:1"
                )
                for i, e, x, y in [
                    (1, "O", 9, 5),
                    (2, "H", 9.96, 5),
                    (3, "H", 8.76, 5.93),
                ]
            ]
            packed["bonds"] = [[1, 2, 1], [1, 3, 1]]
            p.write_text(json.dumps(packed))
            manifest.update(
                water_model="Triatomic numerical verification fixture",
                water_constraints="flexible",
                type_elements={"1": "Ar", "2": "O", "3": "H"},
                bond_orders={"1": 1},
                packed_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),
            )
            m.write_text(json.dumps(manifest))
            data = (
                "Toy charged slab fixture\n\n4 atoms\n2 bonds\n1 angles\n3 atom types\n1 bond types\n1 angle types\n"
                "0 16 xlo xhi\n0 16 ylo yhi\n0 16 zlo zhi\n\nMasses\n\n1 39.948\n2 15.999\n3 1.008\n"
                "\nPair Coeffs\n\n1 .238 3.405\n2 .15 3.1\n3 0 0\n\nBond Coeffs\n\n1 380 .96\n"
                "\nAngle Coeffs\n\n1 42 104.5\n\nAtoms # full\n\n1 1 1 0 5 5 6\n2 2 2 -.8 9 5 6\n"
                "3 2 3 .4 9.96 5 6\n4 2 3 .4 8.76 5.93 6\n\nBonds\n\n1 1 2 3\n2 1 2 4\n\nAngles\n\n1 1 3 2 4\n"
            )
            t.write_text(data)
            model = register(p, t, m)
            r = run_registered(
                dict(topology_id=model["id"], steps=8), root / "run", root
            )
            self.assertTrue(math.isfinite(r["final_kcal_mol"]))
            self.assertIn(
                "kspace_modify slab 3.0", (root / "run" / "registered.in").read_text()
            )
            with self.assertRaises(ValueError):
                validate(packed, data.replace("1 1 3 2 4", "1 1 2 3 4"), manifest)

    def test_contact_groups_use_only_accessible_termination_faces(self):
        from polychemmd.studio.register_topology import _groups

        def slab(z, sheet, layer, elem):
            return dict(z=z, sheet_index=sheet, layer=layer, element=elem, role="slab")

        atoms = [
            slab(4, 9, "middle_Ti", "Ti"),
            slab(2, 9, "bottom_OH", "O"),
            slab(6, 9, "top_F", "F"),
            slab(16, 1, "middle_Ti", "Ti"),
            slab(14, 1, "bottom_OH", "O"),
            slab(14, 1, "bottom_OH", "H"),
            slab(18, 1, "top_F", "F"),
        ]
        slit = _groups(dict(atoms=atoms, geometry="slit"))
        self.assertEqual(slit["reference_plane"], [1])
        self.assertEqual(slit["comparison_plane"], [4])
        self.assertEqual(slit["surface_sites"], [3, 5])
        self.assertEqual(slit["plane_restraint"], "fixed")
        single = _groups(dict(atoms=atoms[:3], geometry="single"))
        self.assertEqual(single["surface_sites"], [3])
        self.assertEqual(single["comparison_plane"], [])
