import importlib.util
from importlib.resources import files
import math
import os
from pathlib import Path
import re
import tempfile
import unittest

from polychemmd.studio.polymer_opls import compile_polymer, source_atom_types


def studio_fixture(sequence):
    """Reconstruct element/connectivity from the independently typed Atlas file."""
    source = (
        files("polychemmd.studio")
        .joinpath("md_assets", f"pva_{sequence}.lt")
        .read_text()
    )
    elements = {
        135: "C",
        136: "C",
        140: "H",
        154: "O",
        155: "H",
        158: "C",
        282: "H",
        465: "C",
        466: "O",
        467: "O",
        469: "H",
        491: "C",
    }
    atoms, types, bonds = [], [], []
    for line in source.splitlines():
        if hit := re.match(
            r"\s*\$atom:a(\d+) \$mol:\. @atom:(\d+) (\S+) (\S+) (\S+) (\S+)", line
        ):
            k = int(hit[2])
            atoms.append(
                dict(
                    element=elements[k],
                    x=float(hit[4]),
                    y=float(hit[5]),
                    z=float(hit[6]),
                )
            )
            types.append(k)
        elif hit := re.match(r"\s*\$bond:b\d+ \$atom:a(\d+) \$atom:a(\d+)", line):
            i, j = int(hit[1]) - 1, int(hit[2]) - 1
            order = 2 if {types[i], types[j]} == {465, 466} else 1
            bonds.append([i, j, order])
    return atoms, bonds, types


def small_chain(smiles):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    AllChem.Compute2DCoords(mol)
    atoms = [
        dict(element=a.GetSymbol(), x=float(x), y=float(y), z=float(z))
        for a, (x, y, z) in zip(mol.GetAtoms(), mol.GetConformer().GetPositions())
    ]
    bonds = [
        (b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble())
        for b in mol.GetBonds()
    ]
    return atoms, bonds


class TestLocalOPLS(unittest.TestCase):
    def test_all_atom_labels_reproduce_both_legacy_sequence_assignments(self):
        for name in ("dispersed", "clustered"):
            with self.subTest(name=name):
                atoms, bonds, reference = studio_fixture(name)
                self.assertEqual(len(atoms), 104)
                self.assertEqual(source_atom_types(atoms, bonds, "PVA"), reference)

    def test_missing_hydrogen_and_duplicate_edge_are_rejected(self):
        atoms, bonds, _ = studio_fixture("clustered")
        with self.assertRaisesRegex(ValueError, "valence"):
            source_atom_types(atoms, bonds[:-1], "PVA")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            source_atom_types(atoms, bonds + [bonds[0]], "PVA")

    @unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit optional")
    def test_dialkyl_ether_and_diol_end_types(self):
        atoms, bonds = small_chain("OCCOCCO")
        types = source_atom_types(atoms, bonds, "PEO")
        self.assertEqual(types[:7], [169, 173, 182, 180, 182, 173, 169])
        self.assertEqual(types.count(185), 4)
        self.assertEqual(types.count(170), 2)
        with self.assertRaises(ValueError):
            source_atom_types(atoms, bonds, "PVA")

    @unittest.skipUnless(importlib.util.find_spec("rdkit"), "RDKit optional")
    def test_missing_tertiary_amide_inter_repeat_torsion_blocks_PEtOx(self):
        atoms, bonds = small_chain("CN(C(=O)CC)CCN(C(=O)CC)CC")
        with tempfile.TemporaryDirectory() as td:
            with self.assertRaisesRegex(ValueError, "NM–CT–CT–NM"):
                compile_polymer(atoms, bonds, "PEtOx", td)
            self.assertEqual(list(Path(td).iterdir()), [])

    @unittest.skipUnless(
        importlib.util.find_spec("moltemplate") and os.name == "posix",
        "optional POSIX Moltemplate compiler",
    )
    def test_moltemplate_compiles_every_pva_interaction_with_exact_atom_identity(self):
        atoms, bonds, ids = studio_fixture("dispersed")
        with tempfile.TemporaryDirectory() as td:
            result = compile_polymer(atoms, bonds, "PVA", td)
        self.assertEqual(
            {k: len(v) for k, v in result["terms"].items()},
            {"Bonds": 103, "Angles": 189, "Dihedrals": 276, "Impropers": 3},
        )
        self.assertEqual(
            [result["types"][t]["source_type"] for t in result["atom_types"]], ids
        )
        self.assertAlmostEqual(
            math.fsum(result["types"][t]["charge"] for t in result["atom_types"]),
            0.0,
            places=12,
        )
        self.assertEqual(result["source"]["version"], "2.23.0")
