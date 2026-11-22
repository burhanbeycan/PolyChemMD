"""Audit the crystal composition, mobile interactions and actual generated topology."""

import importlib.util
import json
import math
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from polychemmd.studio import iff_sources
from polychemmd.studio.topology_build import (
    _audit_sheet,
    _native_key,
    _write_data,
    run_topology,
)
from polychemmd.studio.topology_file import parse_data
from polychemmd.studio.register_topology import run_registered


class GeneratedTopologyTests(unittest.TestCase):
    def test_each_sheet_must_contain_complete_decorations_and_correct_stoichiometry(
        self,
    ):
        atoms = [
            dict(role="slab", sheet_index=0, layer=k, element=e, site=s)
            for k, e, s in [
                ("middle_Ti", "Ti", None),
                ("top_Ti", "Ti", None),
                ("bottom_Ti", "Ti", None),
                ("top_C", "C", None),
                ("bottom_C", "C", None),
                ("top_F", "F", 0),
                ("bottom_OH", "O", 0),
                ("bottom_OH", "H", 0),
            ]
        ]
        _audit_sheet(atoms, [_native_key(a) for a in atoms])
        for missing in (2, 7):
            edited = atoms[:missing] + atoms[missing + 1 :]
            with self.assertRaises(ValueError):
                _audit_sheet(edited, [_native_key(a) for a in edited])

    def test_pair_matrix_keeps_all_mobile_cross_interactions(self):
        # Small numeric fixture with arbitrary LJ values, solely for pair mixing.
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "test.data"
            types = {
                ("IFF", 8): dict(
                    element="O", mass=16.0, epsilon=0.25, sigma=4.0, charge=-0.8
                ),
                ("IFF", 10): dict(
                    element="O", mass=16.0, epsilon=0.36, sigma=9.0, charge=-0.8
                ),
                ("OPLS", 1): dict(
                    element="H", mass=1.0, epsilon=0.0, sigma=0.0, charge=0.4
                ),
            }
            labels = list(types)
            atoms = [
                dict(
                    element=types[k]["element"],
                    molecule_id=m,
                    x=float(i + 1),
                    y=1.0,
                    z=1.0,
                )
                for i, (k, m) in enumerate(zip(labels, ("slab", "water", "chain")))
            ]
            packed = dict(box_A=[20.0, 20.0, 20.0], atoms=atoms)
            sections = ["Bonds", "Angles", "Dihedrals", "Impropers"]
            _write_data(
                path,
                packed,
                types,
                labels,
                {s: {} for s in sections},
                {s: [] for s in sections},
                {},
            )
            rows, _, _ = parse_data(path.read_text(), cvff_impropers=True)
            pairs = {
                (int(a), int(b)): (float(e), float(s))
                for a, b, e, s in rows["PairIJ Coeffs"]
            }
            self.assertEqual(pairs[1, 1], (0.0, 4.0))
            self.assertEqual(pairs[1, 2], (0.3, 6.0))
            self.assertEqual(pairs[2, 2], (0.36, 9.0))
            self.assertEqual(pairs[1, 3], (0.0, 0.0))

    def test_reject_out_of_range_before_reading_inputs(self):
        for cutoff in (True, float("nan"), 0, 5.9, 12.1):
            with (
                self.subTest(cutoff=cutoff),
                self.assertRaisesRegex(ValueError, "cutoff"),
            ):
                run_topology({"cutoff_A": cutoff}, "unused", source_dir="unused")

    def test_nonzero_requested_salt_is_rejected_even_when_rounded_to_no_ions(self):
        from hashlib import sha256

        with TemporaryDirectory() as td:
            base = Path(td)
            data = b'{"polymer":{"name":"bare"}}'
            (base / "structure.json").write_bytes(data)
            (base / "packing_report.json").write_text(
                json.dumps(
                    {
                        "geometry_input_sha256": sha256(data).hexdigest(),
                        "salt_M_requested": 0.001,
                        "salt_M_realized": 0,
                        "ion_pairs": 0,
                    }
                )
            )
            (base / "packed.json").write_text(
                json.dumps(
                    {
                        "pbc": [True, True, False],
                        "units": "angstrom",
                        "box_A": [20, 20, 20],
                        "atoms": [],
                    }
                )
            )
            with self.assertRaisesRegex(ValueError, "salt_M=0"):
                run_topology({}, base / "output", source_dir=base)

    @unittest.skipUnless(
        os.getenv("POLYCHEMMD_TEST_REGISTERED_LAMMPS") == "1"
        and (iff_sources.cache_directory() / "iff-water.data").is_file()
        and all(importlib.util.find_spec(m) for m in ["ase", "rdkit", "moltemplate"]),
        "real LAMMPS/Moltemplate and verified source required",
    )
    def test_generated_pva_reaches_lammps_and_preserves_graph(self):
        from polychemmd.studio.geometry import build_geometry
        from polychemmd.studio.packing import pack_system

        with (
            TemporaryDirectory() as td,
            patch.dict(os.environ, {"POLYCHEMMD_TOPOLOGY_REGISTRY": td + "/registry"}),
        ):
            d = Path(td)
            cell = d / "pack"
            cfg = dict(
                polymer="PVA",
                dp=4,
                hydrolysis=0.5,
                acetate_pattern="clustered",
                seed=21,
                nx=6,
                ny=4,
                chains=2,
                salt_M=0.0,
                water_count=12,
            )
            shape = build_geometry(cfg, cell)
            pack_system(cfg, cell)
            result = run_topology({}, d / "typed", source_dir=cell)
            packed = json.loads((cell / "packed.json").read_text())
            rows, counts, _ = parse_data(
                (d / "typed" / "typed.data").read_text(), cvff_impropers=True
            )
            self.assertEqual(counts["atoms"], len(packed["atoms"]))
            self.assertEqual(result["term_counts"]["Bonds"], len(packed["bonds"]))
            self.assertLess(abs(result["charge_e"]), 1e-8)
            self.assertEqual(result["term_counts"]["Impropers"], 4)
            self.assertEqual(result["chain_count"], 2)
            self.assertEqual(
                result["parameter_provenance"]["polymer"]["acetate_residues_1based"],
                shape["polymer"]["acetate_residues_1based"],
            )
            self.assertTrue((d / "typed" / "atom-types.csv").is_file())
            md = run_registered(
                dict(topology_id=result["topology_id"], minimize_steps=5, steps=2),
                d / "md",
                cell,
            )
            self.assertTrue(math.isfinite(md["final_mobile_fmax_kcal_mol_A"]))
            self.assertLess(md["post_minimize_kcal_mol"], md["initial_kcal_mol"])
            self.assertEqual(md["steps"], 2)
            self.assertTrue((d / "md" / "registered-trajectory.dump").is_file())


if __name__ == "__main__":
    unittest.main()
