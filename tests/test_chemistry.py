import json
from pathlib import Path
import subprocess
import sys
import unittest

from polychemmd.chemistry import AtomLineage, RepeatUnitSpec, UnsupportedChemistryError

try:
    import rdkit
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


class OptionalChemistryTests(unittest.TestCase):
    def test_imports_remain_lazy_when_rdkit_is_unavailable(self):
        source = Path(__file__).resolve().parents[1] / "src"
        script = f"""
import sys
sys.path.insert(0, {str(source)!r})
class BlockRDKit:
    def find_spec(self, fullname, path=None, target=None):
        if fullname == 'rdkit' or fullname.startswith('rdkit.'):
            raise ImportError('intentionally unavailable')
sys.meta_path.insert(0, BlockRDKit())
import polychemmd
from polychemmd.chemistry import RepeatUnitSpec, ChemistryDependencyError
from polychemmd.builders import build_linear
assert 'rdkit' not in sys.modules
try:
    RepeatUnitSpec('PE', '[*:1]CC[*:2]', 1, 2)
except ChemistryDependencyError:
    pass
else:
    raise AssertionError('missing optional dependency was silently ignored')
"""
        completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
        self.assertEqual(completed.returncode, 0, completed.stderr)


@unittest.skipUnless(HAS_RDKIT, "optional RDKit dependency not installed")
class RepeatUnitContractTests(unittest.TestCase):
    def spec(self, smiles="[*:1]CC[*:2]", **kwargs):
        return RepeatUnitSpec("unit", smiles, kwargs.get("head_map", 1), kwargs.get("tail_map", 2))

    def test_lossless_roundtrip_and_fresh_graph(self):
        data = {"label": "mapped PE", "smiles": "[*:11][CH2:7]C[*:22]", "head_map": 11, "tail_map": 22}
        spec = RepeatUnitSpec.from_mapping(data)
        self.assertEqual(spec.to_dict(), json.loads(json.dumps(spec.to_dict())))
        self.assertEqual(RepeatUnitSpec.from_mapping(spec.to_dict()), spec)
        mol = spec.to_mol()
        mol.GetAtomWithIdx(1).SetAtomMapNum(99)
        self.assertEqual(spec.to_mol().GetAtomWithIdx(1).GetAtomMapNum(), 7)
        self.assertEqual(spec.attachment_indices(), {"head_dummy": 0, "head_anchor": 1, "tail_dummy": 3, "tail_anchor": 2})

    def test_strict_object_keys_roles_and_identifiers(self):
        good = {"label": "PE", "smiles": "[*:1]CC[*:2]", "head_map": 1, "tail_map": 2}
        for key in good:
            bad = dict(good)
            del bad[key]
            with self.subTest(missing=key), self.assertRaisesRegex(ValueError, "missing keys"):
                RepeatUnitSpec.from_mapping(bad)
        for key, value in [("label", " "), ("smiles", ""), ("head_map", True), ("tail_map", 2.0), ("head_map", 0), ("tail_map", 1)]:
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                RepeatUnitSpec.from_mapping({**good, key: value})
        with self.assertRaisesRegex(ValueError, "unknown keys"):
            RepeatUnitSpec.from_mapping({**good, "monomer_reaction": "addition"})
        with self.assertRaisesRegex(ValueError, "object"):
            RepeatUnitSpec.from_mapping([])

    def test_free_monomer_is_not_transformed_implicitly(self):
        with self.assertRaisesRegex(ValueError, "free monomers"):
            self.spec("C=C")

    def test_single_connected_repeat_and_two_terminal_dummies(self):
        for smiles, pattern in [
            ("[*:1]CC[*:2].C", "connected"),
            ("[*:1]C([*:2])[*:3]", "exactly two"),
            ("C[*:1](C)CC[*:2]", "terminal"),
            ("[*:1][*:2]", "non-dummy"),
        ]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, pattern):
                self.spec(smiles)

    def test_all_nonzero_maps_unique_and_dummy_roles_match(self):
        for smiles, pattern in [
            ("[*:1]CC[*:1]", "unique"),
            ("[*:1][CH2:9][CH2:9][*:2]", "unique"),
            ("[*:1][CH2:2]C[*:2]", "unique"),
            ("*CC[*:2]", "specified"),
            ("[*:7]CC[*:2]", "specified"),
        ]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, pattern):
                self.spec(smiles)

    def test_invalid_valence_and_unsupported_attachment_bonds(self):
        for smiles in ["[*:1][CH4]C[*:2]", "this is not smiles"]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, "invalid|valence"):
                self.spec(smiles)
        for smiles in ["[*:1]=CC[*:2]", "[*:1]#CC[*:2]"]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(UnsupportedChemistryError, "single-bond"):
                self.spec(smiles)

    def test_names_and_cxsmiles_extensions_are_not_silently_discarded(self):
        for smiles in ["[*:1]CC[*:2] ignored name", "[*:1]CC[*:2] |atomProp:1.foo.bar|"]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, "invalid"):
                self.spec(smiles)

    def test_dummy_charge_isotope_and_radicals_rejected(self):
        for smiles in ["[*+:1]CC[*:2]", "[1*:1]CC[*:2]"]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(ValueError, "neutral terminal"):
                self.spec(smiles)
        with self.assertRaisesRegex(UnsupportedChemistryError, "radical"):
            self.spec("[*:1][CH]C[*:2]")

    def test_attachment_stereochemistry_is_explicitly_unsupported(self):
        for smiles in ["[*:1][C@H](F)CC[*:2]", "[*:1]CC[C@@H](F)[*:2]"]:
            with self.subTest(smiles=smiles), self.assertRaisesRegex(UnsupportedChemistryError, "attachment anchor"):
                self.spec(smiles)
        with self.assertRaisesRegex(UnsupportedChemistryError, "double-bond stereochemistry"):
            self.spec("[*:1]/C=C/[*:2]")
        with self.assertRaisesRegex(UnsupportedChemistryError, "non-tetrahedral"):
            self.spec("[*:1]C[Pt@SP1](Cl)(Br)C[*:2]")

    def test_internal_stereo_is_accepted_without_a_tacticity_claim(self):
        spec = self.spec("[*:1]CC[C@H](F)CCC[*:2]")
        self.assertIn("@", spec.smiles)
        self.assertNotIn("tacticity", spec.to_dict())

    def test_unique_path_roles_allow_sidegroups_and_sidegroup_rings(self):
        spec = self.spec("[*:1]CC(C)C[*:2]")
        roles = spec.structural_roles()
        self.assertEqual(roles["method"], "unique_acyclic_attachment_path")
        self.assertEqual(roles["backbone_source_atom_indices"], [1, 2, 4])
        self.assertEqual(roles["atom_roles"][3], "sidechain")
        self.assertEqual(roles["atom_roles"][0], "attachment")
        ring_sidegroup = self.spec("[*:1]CC(c1ccccc1)C[*:2]").structural_roles()
        self.assertEqual(ring_sidegroup["method"], "unique_acyclic_attachment_path")
        self.assertEqual(list(ring_sidegroup["atom_roles"].values()).count("sidechain"), 6)

    def test_cyclic_backbone_paths_are_marked_ambiguous(self):
        roles = self.spec("[*:1]c1ccc([*:2])cc1").structural_roles()
        self.assertEqual(roles["method"], "ambiguous_cyclic_path")
        self.assertEqual(roles["backbone_source_atom_indices"], [])
        self.assertEqual(list(roles["atom_roles"].values()).count("unassigned"), 6)


class AtomLineageTests(unittest.TestCase):
    def test_json_roundtrip_preserves_copy_source_and_current_mapping(self):
        item = AtomLineage("repeat:2:atom:5", 7, 8, "repeat", 2, 5, 77, None, None, "backbone")
        encoded = json.loads(json.dumps(item.to_dict()))
        self.assertEqual(AtomLineage.from_mapping(encoded), item)

    def test_invalid_or_incomplete_lineage_rejected(self):
        good = AtomLineage("repeat:2:atom:5", 7, 8, "repeat", 2, 5, 77, None, None, "backbone").to_dict()
        for key, value in [("atom_map", True), ("atom_index", -1), ("repeat_index", None), ("source_atom_map", 0), ("structural_role", "guessed"), ("terminus", "head")]:
            with self.subTest(key=key), self.assertRaises(ValueError):
                AtomLineage.from_mapping({**good, key: value})
        with self.assertRaisesRegex(ValueError, "missing keys"):
            AtomLineage.from_mapping({"atom_id": "x"})
