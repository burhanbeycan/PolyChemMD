from collections import Counter
import json
import unittest

from polychemmd.chemistry import AtomLineage, RepeatUnitSpec, UnsupportedChemistryError
from polychemmd.builders import build_linear

try:
    from rdkit import Chem
    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False


@unittest.skipUnless(HAS_RDKIT, "optional RDKit dependency not installed")
class LinearBuilderTests(unittest.TestCase):
    def spec(self, smiles="[*:1]CC[*:2]"):
        return RepeatUnitSpec("test unit", smiles, 1, 2)

    def canonical(self, smiles):
        return Chem.MolToSmiles(Chem.MolFromSmiles(smiles), isomericSmiles=True)

    def test_hand_checked_pe_formula_and_identity_for_one_two_several_units(self):
        # n ethylene repeat graphs contain C_(2n)H_(4n); two H caps add H2.
        for dp in (1, 2, 5, 12):
            with self.subTest(dp=dp):
                result = build_linear(self.spec(), dp)
                self.assertEqual(result.molecular_formula, f"C{2 * dp}H{4 * dp + 2}")
                self.assertEqual(result.canonical_smiles, "C" * (2 * dp))
                self.assertEqual(result.mol.GetNumAtoms(), 6 * dp + 2)
                self.assertEqual(result.formal_charge, 0)
                self.assertEqual(result.to_dict()["element_counts"], {"C": 2 * dp, "H": 4 * dp + 2})

    def test_manual_methyl_and_mixed_cap_atom_balance(self):
        for dp in (1, 2, 5):
            for caps, extra_c, extra_h in [(('methyl', 'methyl'), 2, 6), (('methyl', 'H'), 1, 4), (('H', 'methyl'), 1, 4)]:
                with self.subTest(dp=dp, caps=caps):
                    result = build_linear(self.spec(), dp, *caps)
                    self.assertEqual(result.molecular_formula, f"C{2 * dp + extra_c}H{4 * dp + extra_h}")
                    self.assertEqual(result.canonical_smiles, "C" * (2 * dp + extra_c))
                    self.assertEqual(result.mol.GetNumAtoms(), 6 * dp + extra_c + extra_h)

    def test_head_tail_direction_and_heteroatom_linkage(self):
        # *-CH2-CH2-O-* gives H-CH2-CH2-O-CH2-CH2-O-H at DP=2.
        result = build_linear(self.spec("[*:1]CCO[*:2]"), 2)
        self.assertEqual(result.molecular_formula, "C4H10O2")
        self.assertEqual(result.canonical_smiles, self.canonical("CCOCCO"))
        # A one-carbon backbone is also a valid two-ended graph, not a branch.
        self.assertEqual(build_linear(self.spec("[*:1]C[*:2]"), 2).canonical_smiles, "CC")

    def test_sidebranches_remain_sidebranches(self):
        result = build_linear(self.spec("[*:1]CC(C)C[*:2]"), 2)
        self.assertEqual(result.molecular_formula, "C8H18")
        self.assertEqual(result.canonical_smiles, self.canonical("CC(C)CCC(C)C"))
        carbons = [a for a in result.to_dict()["atoms"] if a["element"] == "C"]
        self.assertEqual(Counter(a["structural_role"] for a in carbons), {"backbone": 6, "sidechain": 2})

    def test_formal_charge_multiplies_but_is_not_a_partial_charge_model(self):
        result = build_linear(self.spec("[*:1]C[NH2+]C[*:2]"), 3)
        self.assertEqual(result.formal_charge, 3)
        self.assertEqual(result.molecular_formula, "C6H20N3+3")
        self.assertNotIn("partial_charges", result.to_dict())

    def test_lineage_is_unique_stable_and_roundtrips(self):
        spec = self.spec("[*:1][CH2:7]C[*:2]")
        first, longer = build_linear(spec, 2), build_linear(spec, 3)
        self.assertEqual(len(first.lineage), len({a.atom_id for a in first.lineage}))
        self.assertEqual([a.atom_map for a in first.lineage], list(range(1, first.mol.GetNumAtoms() + 1)))
        roundtrip = [AtomLineage.from_mapping(json.loads(json.dumps(a.to_dict()))) for a in first.lineage]
        self.assertEqual(roundtrip, list(first.lineage))
        old_repeat = {a.atom_id: (a.repeat_index, a.source_atom_index, a.source_atom_map) for a in first.lineage if a.origin == "repeat"}
        new_repeat = {a.atom_id: (a.repeat_index, a.source_atom_index, a.source_atom_map) for a in longer.lineage if a.origin == "repeat"}
        for atom_id, source in old_repeat.items():
            self.assertEqual(new_repeat[atom_id], source)
        self.assertEqual([a.source_atom_map for a in first.lineage if a.origin == "repeat"], [7, None, 7, None])
        ids = {a.atom_id for a in first.lineage}
        for atom in first.lineage:
            if atom.origin == "added_hydrogen":
                self.assertIn(atom.parent_atom_id, ids)

    def test_explicit_source_hydrogen_isotopes_and_source_indices_preserved(self):
        result = build_linear(self.spec("[*:1][13CH2]C([2H:9])([H:8])[*:2]"), 2)
        isotope_atoms = [(a.GetAtomicNum(), a.GetIsotope()) for a in result.mol.GetAtoms() if a.GetIsotope()]
        self.assertEqual(Counter(isotope_atoms), {(6, 13): 2, (1, 2): 2})
        self.assertEqual(result.molecular_formula, "C2[13C]2H8D2")
        source_hydrogens = [a for a in result.to_dict()["atoms"] if a["origin"] == "repeat" and a["element"] == "H"]
        self.assertEqual(len(source_hydrogens), 4)
        self.assertEqual(Counter(a["source_atom_map"] for a in source_hydrogens), {8: 2, 9: 2})

    def test_structural_roles_do_not_depend_on_explicit_source_hydrogen_notation(self):
        implicit = build_linear(self.spec(), 1)
        explicit = build_linear(self.spec("[H]C([H])([*:1])C[*:2]"), 1)
        self.assertEqual(implicit.canonical_smiles, explicit.canonical_smiles)
        for result in (implicit, explicit):
            mol = result.mol
            for atom, item in zip(mol.GetAtoms(), result.lineage):
                if atom.GetAtomicNum() == 1 and item.origin != "cap":
                    parent = atom.GetNeighbors()[0]
                    self.assertEqual(item.structural_role, result.lineage[parent.GetIdx()].structural_role)
            counts = Counter((a["element"], a["structural_role"]) for a in result.to_dict()["atoms"])
            self.assertEqual(counts, {("C", "backbone"): 2, ("H", "backbone"): 4, ("H", "terminus"): 2})

    def test_canonical_identity_is_unmapped_and_separate_from_mapped_graph(self):
        result = build_linear(self.spec(), 2)
        self.assertNotIn(":", result.canonical_smiles)
        parameters = Chem.SmilesParserParams()
        parameters.removeHs = False
        mapped = Chem.MolFromSmiles(result.mapped_smiles, parameters)
        maps = [a.GetAtomMapNum() for a in mapped.GetAtoms()]
        self.assertEqual(sorted(maps), list(range(1, result.mol.GetNumAtoms() + 1)))
        self.assertEqual(mapped.GetNumAtoms(), result.mol.GetNumAtoms())
        self.assertFalse(any(a.GetAtomicNum() == 0 for a in mapped.GetAtoms()))
        self.assertEqual(len(Chem.GetMolFrags(mapped)), 1)

    def test_deterministic_result_and_mol_defensive_copy(self):
        first = build_linear(self.spec(), 4)
        self.assertEqual(first.to_dict(), build_linear(self.spec(), 4).to_dict())
        copy = first.mol
        copy.GetAtomWithIdx(0).SetFormalCharge(4)
        self.assertEqual(first.to_dict()["formal_charge"], 0)
        self.assertEqual(first.mol.GetAtomWithIdx(0).GetFormalCharge(), 0)
        self.assertEqual(first.to_dict(), json.loads(json.dumps(first.to_dict())))

    def test_remote_tetrahedral_stereochemistry_matches_hand_written_product(self):
        spec = self.spec("[*:1]CC[C@H](F)CCC[*:2]")
        single = build_linear(spec, 1)
        self.assertEqual(single.canonical_smiles, self.canonical("CC[C@H](F)CCC"))
        dimer = build_linear(spec, 2)
        self.assertEqual(dimer.canonical_smiles, self.canonical("CC[C@H](F)CCCCC[C@H](F)CCC"))
        self.assertEqual(sum(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in dimer.mol.GetAtoms()), 2)
        self.assertNotEqual(dimer.canonical_smiles, self.canonical("CC[C@@H](F)CCCCC[C@@H](F)CCC"))

    def test_internal_double_bond_geometry_matches_hand_written_product(self):
        result = build_linear(self.spec("[*:1]CC/C=C/CC[*:2]"), 2)
        self.assertEqual(result.canonical_smiles, self.canonical("CC/C=C/CCCC/C=C/CC"))
        self.assertEqual(sum(b.GetStereo() == Chem.BondStereo.STEREOE for b in result.mol.GetBonds()), 2)

    def test_stereo_that_becomes_undefined_under_caps_is_explicitly_rejected(self):
        # Identical capped arms at DP=1 remove the stated repeat-center chirality.
        with self.assertRaisesRegex(UnsupportedChemistryError, "stereochemistry is lost"):
            build_linear(self.spec("[*:1]C[C@H](C)C[*:2]"), 1)

    def test_strict_dp_and_configured_safety_limits(self):
        for value in (True, False, 0, -1, 1.0, "2", None, float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaisesRegex(ValueError, "positive integer"):
                build_linear(self.spec(), value)
        for kwargs in ({"max_atoms": True}, {"max_degree_of_polymerization": 0}):
            with self.subTest(kwargs=kwargs), self.assertRaisesRegex(ValueError, "positive integer"):
                build_linear(self.spec(), 1, **kwargs)
        with self.assertRaisesRegex(ValueError, "degree_of_polymerization exceeds"):
            build_linear(self.spec(), 3, max_degree_of_polymerization=2)
        with self.assertRaisesRegex(ValueError, "max_atoms"):
            build_linear(self.spec(), 2, max_atoms=13)
        self.assertEqual(build_linear(self.spec(), 2, max_atoms=14).mol.GetNumAtoms(), 14)

    def test_caps_and_input_type_fail_explicitly(self):
        for value in ("OH", "CH3", "h", None, ["H"]):
            with self.subTest(value=value), self.assertRaisesRegex(UnsupportedChemistryError, "H or methyl"):
                build_linear(self.spec(), 1, head_cap=value)
        with self.assertRaisesRegex(ValueError, "RepeatUnitSpec"):
            build_linear({"smiles": "C=C"}, 1)

    def test_output_claims_only_graph_assembly(self):
        result = build_linear(self.spec(), 2).to_dict()
        self.assertEqual(result["operation"], "repeat_unit_assembly")
        self.assertEqual(result["evidence_level"], "graph_validated")
        self.assertFalse(result["validation"]["geometry_generated"])
        self.assertFalse(result["validation"]["force_field_assigned"])
        for key in ("density", "coordinates", "reaction", "tacticity", "energy"):
            self.assertNotIn(key, result)
