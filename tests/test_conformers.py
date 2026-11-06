import hashlib
import importlib.util
import json
import math
import unittest
from unittest import mock

from polychemmd.builders.linear import build_linear
from polychemmd.builders.conformers import (
    ConformerGenerationError, _check_embedded_geometry, _verify_stereochemistry,
    generate_conformer,
)
from polychemmd.chemistry import RepeatUnitSpec
from polychemmd.geometry import MolecularGeometry


HAS_RDKIT = importlib.util.find_spec("rdkit") is not None


@unittest.skipUnless(HAS_RDKIT, "optional RDKit dependency not installed")
class ConformerTests(unittest.TestCase):
    def chain(self, smiles="[*:1]CC[*:2]", dp=3):
        return build_linear(RepeatUnitSpec("conformer test", smiles, 1, 2), dp)

    def embedded(self, chain, seed=17):
        from rdkit.Chem import rdDistGeom
        source = chain.mol
        embedded = chain.mol
        parameters = rdDistGeom.ETKDGv3()
        parameters.randomSeed = seed
        parameters.numThreads = 1
        self.assertEqual(rdDistGeom.EmbedMolecule(embedded, parameters), 0)
        return source, embedded

    def test_seeded_real_geometry_preserves_graph_and_lineage(self):
        chain = self.chain(dp=5)
        graph = chain.to_dict()
        geometry = generate_conformer(chain, seed=43)
        self.assertEqual(len(geometry.atoms), 32)
        self.assertEqual(len(geometry.bonds), 31)
        self.assertEqual([a.atom_id for a in geometry.atoms], [a["atom_id"] for a in graph["atoms"]])
        self.assertEqual([json.loads(a.lineage_json) for a in geometry.atoms], graph["atoms"])
        self.assertEqual([a.atom_map for a in geometry.atoms], list(range(1, 33)))
        self.assertEqual(chain.to_dict(), graph)
        self.assertEqual(chain.mol.GetNumConformers(), 0)
        expected = hashlib.sha256(json.dumps(graph, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest()
        self.assertEqual(geometry.source_graph_sha256, expected)
        restored = MolecularGeometry.from_mapping(json.loads(json.dumps(geometry.to_dict())))
        self.assertEqual(restored, geometry)

    def test_hand_checked_ethane_geometry_and_mass(self):
        geometry = generate_conformer(self.chain(dp=1), seed=17)
        carbons = [i for i, atom in enumerate(geometry.atoms) if atom.element == "C"]
        self.assertEqual(len(carbons), 2)
        distance = math.dist(*(geometry.coordinates_angstrom[i] for i in carbons))
        # Independent broad known molecular dimensions, not a fitted reference.
        self.assertGreater(distance, 1.35)
        self.assertLess(distance, 1.7)
        for bond in geometry.bonds:
            atom_by_id = {a.atom_id: i for i, a in enumerate(geometry.atoms)}
            indices = [atom_by_id[identity] for identity in bond.atom_ids]
            if any(geometry.atoms[i].element == "H" for i in indices):
                self.assertGreater(math.dist(*(geometry.coordinates_angstrom[i] for i in indices)), .9)
                self.assertLess(math.dist(*(geometry.coordinates_angstrom[i] for i in indices)), 1.25)
        self.assertAlmostEqual(sum(a.mass_g_mol for a in geometry.atoms), 2 * 12.011 + 6 * 1.008, places=6)

    def test_identical_seed_is_reproducible_and_distinct_seed_changes_geometry(self):
        chain = self.chain(dp=5)
        first = generate_conformer(chain, seed=0)
        second = generate_conformer(chain, seed=0)
        other = generate_conformer(chain, seed=193)
        self.assertEqual(first.to_dict(), second.to_dict())
        self.assertNotEqual(first.coordinates_angstrom, other.coordinates_angstrom)
        self.assertEqual(first.source_graph_sha256, other.source_graph_sha256)

    def test_isotope_and_charge_are_preserved_without_partial_charge_claim(self):
        geometry = generate_conformer(self.chain("[*:1][13CH2]C([2H])[*:2]", dp=2), seed=19)
        self.assertEqual(sum(a.isotope == 13 for a in geometry.atoms), 2)
        self.assertEqual(sum(a.isotope == 2 for a in geometry.atoms), 2)
        self.assertAlmostEqual(next(a.mass_g_mol for a in geometry.atoms if a.isotope == 13), 13.00335484, places=6)
        ionic = generate_conformer(self.chain("[*:1]C[NH2+]C[*:2]", dp=2), seed=19)
        self.assertEqual(sum(a.formal_charge for a in ionic.atoms), 2)
        self.assertFalse(json.loads(ionic.provenance_json)["force_field_assigned"])
        self.assertNotIn("partial_charges", ionic.to_dict())

    def test_specified_tetrahedral_stereo_checked_from_actual_3d(self):
        geometry = generate_conformer(self.chain("[*:1]CC[C@H](F)CCC[*:2]", dp=2), seed=17)
        checks = json.loads(geometry.provenance_json)["validation"]["stereochemistry"]
        self.assertEqual(checks["specified_tetrahedral_centers_checked"], 2)
        self.assertEqual(checks["specified_double_bonds_checked"], 0)

    def test_mirror_coordinates_are_rejected_despite_copied_atom_tags(self):
        chain = self.chain("[*:1]CC[C@H](F)CCC[*:2]", dp=2)
        source, embedded = self.embedded(chain)
        conformer = embedded.GetConformer()
        for i in range(embedded.GetNumAtoms()):
            point = conformer.GetAtomPosition(i)
            conformer.SetAtomPosition(i, (-point.x, point.y, point.z))
        # Chemical atom tags still agree, yet the actual coordinates are mirrored.
        self.assertEqual([a.GetChiralTag() for a in source.GetAtoms()], [a.GetChiralTag() for a in embedded.GetAtoms()])
        with self.assertRaisesRegex(ConformerGenerationError, "tetrahedral stereo"):
            _verify_stereochemistry(source, embedded)

    def test_specified_e_and_z_double_bond_geometry(self):
        for smiles in ("[*:1]CC/C=C/CC[*:2]", "[*:1]CC/C=C\\CC[*:2]"):
            with self.subTest(smiles=smiles):
                geometry = generate_conformer(self.chain(smiles, dp=2), seed=17)
                checks = json.loads(geometry.provenance_json)["validation"]["stereochemistry"]
                self.assertEqual(checks["specified_double_bonds_checked"], 2)

    def test_inverted_e_z_coordinates_are_rejected_despite_copied_bond_tags(self):
        from rdkit import Chem
        source, embedded = self.embedded(self.chain("[*:1]CC/C=C/CC[*:2]", dp=1))
        bond = next(b for b in embedded.GetBonds() if b.GetStereo() == Chem.BondStereo.STEREOE)
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        affected, pending = {right}, [right]
        while pending:
            index = pending.pop()
            for neighbor in embedded.GetAtomWithIdx(index).GetNeighbors():
                other = neighbor.GetIdx()
                if other != left and other not in affected:
                    affected.add(other)
                    pending.append(other)
        conformer = embedded.GetConformer()
        origin = tuple(conformer.GetAtomPosition(right))
        axis = tuple(a-b for a, b in zip(conformer.GetAtomPosition(right), conformer.GetAtomPosition(left)))
        norm = math.sqrt(sum(v*v for v in axis))
        axis = tuple(v/norm for v in axis)
        for i in affected:
            vector = tuple(a-b for a, b in zip(conformer.GetAtomPosition(i), origin))
            projection = sum(a*b for a, b in zip(vector, axis))
            conformer.SetAtomPosition(i, tuple(o + 2*projection*u-v for o,u,v in zip(origin,axis,vector)))
        self.assertEqual(embedded.GetBondWithIdx(bond.GetIdx()).GetStereo(), Chem.BondStereo.STEREOE)
        with self.assertRaisesRegex(ConformerGenerationError, "double-bond stereo"):
            _verify_stereochemistry(source, embedded)

    def test_unspecified_stereo_is_not_assigned_to_chemical_source(self):
        chain = self.chain("[*:1]CC(C)CCC[*:2]", dp=1)
        before = chain.to_dict()
        geometry = generate_conformer(chain, seed=17)
        stereo = json.loads(geometry.provenance_json)["validation"]["stereochemistry"]
        self.assertEqual(stereo["specified_tetrahedral_centers_checked"], 0)
        self.assertEqual(stereo["unspecified_stereochemistry"], "not_assigned_to_source_graph")
        self.assertEqual(chain.to_dict(), before)

    def test_collapsed_nonfinite_and_2d_coordinates_are_rejected(self):
        for failure in ("collapse", "nonfinite", "2d"):
            with self.subTest(failure=failure):
                source, embedded = self.embedded(self.chain(dp=1))
                conformer = embedded.GetConformer()
                if failure == "collapse":
                    conformer.SetAtomPosition(0, conformer.GetAtomPosition(1))
                elif failure == "nonfinite":
                    conformer.SetAtomPosition(0, (float("nan"), 0, 0))
                else:
                    conformer.Set3D(False)
                with self.assertRaises(ConformerGenerationError):
                    _check_embedded_geometry(source, embedded)

    def test_bond_distance_sanity_rejects_a_detached_atom(self):
        source, embedded = self.embedded(self.chain(dp=1))
        embedded.GetConformer().SetAtomPosition(0, (100, 100, 100))
        with self.assertRaisesRegex(ConformerGenerationError, "covalent-distance"):
            _check_embedded_geometry(source, embedded)

    def test_embedding_failure_is_explicit_not_a_placeholder_geometry(self):
        chain = self.chain()
        with mock.patch("rdkit.Chem.rdDistGeom.EmbedMolecule", return_value=-1), self.assertRaisesRegex(ConformerGenerationError, "did not produce"):
            generate_conformer(chain, seed=17)
        with mock.patch("rdkit.Chem.rdDistGeom.EmbedMolecule", side_effect=RuntimeError("controlled backend failure")), self.assertRaisesRegex(ConformerGenerationError, "coordinate generation failed"):
            generate_conformer(chain, seed=17)

    def test_backend_uses_recorded_bounds_seed_and_single_thread(self):
        from rdkit.Chem import rdDistGeom
        actual = rdDistGeom.EmbedMolecule
        observed = {}
        def inspect(molecule, parameters):
            observed.update(seed=parameters.randomSeed, iterations=parameters.maxIterations,
                            threads=parameters.numThreads, timeout=parameters.timeout,
                            chirality=parameters.enforceChirality)
            return actual(molecule, parameters)
        with mock.patch("rdkit.Chem.rdDistGeom.EmbedMolecule", side_effect=inspect):
            generate_conformer(self.chain(dp=1), seed=79, max_iterations=27)
        self.assertEqual(observed, {"seed":79, "iterations":27, "threads":1, "timeout":30, "chirality":True})

    def test_limits_and_numeric_input_are_strict(self):
        chain = self.chain(dp=1)
        for seed in (True, False, -1, 1.0, "17", None, 2**31, float("nan")):
            with self.subTest(seed=seed), self.assertRaisesRegex(ValueError, "seed must be an integer"):
                generate_conformer(chain, seed=seed)
        for kwargs in ({"max_iterations":0}, {"max_iterations":True}, {"max_iterations":2**31},
                       {"max_atoms":False}, {"max_atoms":0}, {"max_atoms":1.5}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                generate_conformer(chain, seed=17, **kwargs)
        with self.assertRaisesRegex(ValueError, "exceeds max_atoms"):
            generate_conformer(chain, seed=17, max_atoms=7)
        with self.assertRaisesRegex(ValueError, "LinearChainResult"):
            generate_conformer({"smiles":"CC"}, seed=17)

    def test_evidence_states_geometry_only(self):
        provenance = json.loads(generate_conformer(self.chain(dp=1), seed=17).provenance_json)
        self.assertEqual(provenance["method"], "RDKit_ETKDGv3")
        self.assertEqual(provenance["parameters"]["ETversion"], 2)
        for key in ("force_field_assigned", "energy_minimization_performed", "conformational_search_performed", "simulation_performed"):
            self.assertFalse(provenance[key])
        self.assertEqual(provenance["material_property_validation"], "not_performed")


if __name__ == "__main__":
    unittest.main()
