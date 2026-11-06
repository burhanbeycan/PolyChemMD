from collections import Counter
from dataclasses import replace
from itertools import combinations
import json
import math
import unittest
from unittest.mock import patch

from polychemmd.forcefields.alkanes import AssignedAlkane, UnsupportedAlkaneError, assign_alkane
from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry


def graph_fixture(carbon_count, carbon_edges=None):
    """Explicit hand-defined graph; coordinates are not a computed geometry."""
    if carbon_edges is None:
        carbon_edges = [(i, i+1) for i in range(carbon_count-1)]
    degree = Counter(endpoint for edge in carbon_edges for endpoint in edge)
    atoms = [GeometryAtom(f'carbon:{i}',i+1,'C',0,0,12.011,'{"origin":"analytic_graph_fixture"}')
             for i in range(carbon_count)]
    points = [(1.5*i,0.,0.) for i in range(carbon_count)]
    bonds = [GeometryBond((f'carbon:{first}',f'carbon:{second}'),1)
             for first,second in carbon_edges]
    for i in range(carbon_count):
        for h in range(4-degree[i]):
            name = f'H:{i}:{h}'
            atoms.append(GeometryAtom(name,len(atoms)+1,'H',0,0,1.008,'{"origin":"analytic_graph_fixture"}'))
            points.append((1.5*i+0.1*(h+1),1.+(h%2),.3*h))
            bonds.append(GeometryBond((f'carbon:{i}',name),1))
    return MolecularGeometry(tuple(atoms),tuple(bonds),tuple(points),'0'*64,
                             '{"source_kind":"analytic_topology_fixture","geometry_validation":"not_performed"}')


def replace_atom(geometry, index, **changes):
    atoms = list(geometry.atoms)
    atoms[index] = replace(atoms[index],**changes)
    return replace(geometry,atoms=tuple(atoms))


class AlkaneAssignmentTests(unittest.TestCase):
    def test_ethane_hand_counted_complete_topology_and_pair_scaling(self):
        record=assign_alkane(graph_fixture(2)).to_dict()
        self.assertEqual(Counter(a['type_id'] for a in record['atom_assignments']),
                         {'opls_135':2,'opls_140':6})
        self.assertEqual([len(record['interactions'][kind]) for kind in ('bonds','angles','dihedrals','impropers')],
                         [7,12,9,0])
        self.assertEqual(Counter(d['type_id'] for d in record['interactions']['dihedrals']),
                         {'torsion_HCCH':9})
        self.assertEqual(len(record['nonbonded_pairs']),28)
        self.assertEqual(Counter(pair['topological_distance'] for pair in record['nonbonded_pairs']),
                         {1:7,2:12,3:9})
        for pair in record['nonbonded_pairs']:
            expected=.5 if pair['topological_distance']==3 else 0.
            self.assertEqual(pair['lj_scale'],expected)
            self.assertEqual(pair['coulomb_scale'],expected)

    def test_butane_retains_hydrogen_ended_torsions_and_all_angles(self):
        record=assign_alkane(graph_fixture(4)).to_dict()
        self.assertEqual(Counter(a['type_id'] for a in record['atom_assignments']),
                         {'opls_135':2,'opls_136':2,'opls_140':10})
        self.assertEqual([len(record['interactions'][kind]) for kind in ('bonds','angles','dihedrals')],
                         [13,24,27])
        self.assertEqual(Counter(d['type_id'] for d in record['interactions']['dihedrals']),
                         {'torsion_CCCC':1,'torsion_CCCH':10,'torsion_HCCH':16})
        self.assertEqual(Counter(a['type_id'] for a in record['interactions']['angles']),
                         {'angle_CCC':2,'angle_CCH':14,'angle_HCH':8})
        self.assertEqual(record['validation']['pair_scaling_counts'],
                         {'excluded':37,'scaled_14':27,'unscaled':27})
        # Two terminal carbon atoms are three graph edges apart, irrespective
        # of their deliberately nonphysical fixture-coordinate separation.
        selected=next(p for p in record['nonbonded_pairs'] if p['atom_ids']==['carbon:0','carbon:3'])
        self.assertEqual(selected['topological_distance'],3)
        self.assertEqual(selected['lj_scale'],.5)

    def test_complete_coverage_and_charge_balance_have_scoped_evidence(self):
        for n in (2,3,4,10):
            with self.subTest(carbon_count=n):
                record=assign_alkane(graph_fixture(n)).to_dict()
                self.assertEqual(len(record['atom_assignments']),3*n+2)
                self.assertEqual(len(record['interactions']['bonds']),3*n+1)
                self.assertEqual(len(record['interactions']['angles']),6*n)
                self.assertEqual(len(record['interactions']['dihedrals']),9*(n-1))
                self.assertAlmostEqual(math.fsum(a['charge_e'] for a in record['atom_assignments']),0.,places=14)
                self.assertTrue(record['coverage']['assigned_type_coverage_complete'])
                self.assertTrue(record['validation']['parameters_assigned'])
                self.assertFalse(record['validation']['production_ready'])
                self.assertEqual(record['validation']['physical_accuracy'],'not_assessed')
                self.assertEqual(record['parameter_set']['source']['kind'],'literature')

    def test_input_geometry_lineage_and_order_are_preserved(self):
        geometry=graph_fixture(4)
        before=geometry.to_dict()
        assigned=assign_alkane(geometry)
        self.assertEqual(geometry.to_dict(),before)
        self.assertEqual(assigned.to_dict()['geometry'],before)
        self.assertEqual(assigned.geometry,geometry)
        out=assigned.to_dict()
        out['geometry']['atoms'][0]['lineage']['origin']='modified'
        out['atom_assignments'][0]['charge_e']=99
        self.assertEqual(assigned.to_dict()['geometry'],before)
        self.assertEqual(assigned.to_dict()['atom_assignments'][0]['charge_e'],-.18)
        self.assertEqual(assign_alkane(geometry).to_dict(),assigned.to_dict())

    def test_topology_interactions_are_unique_complete_graph_paths(self):
        geometry=graph_fixture(5)
        record=assign_alkane(geometry).to_dict()
        edges={frozenset(b.atom_ids) for b in geometry.bonds}
        for kind in ('bonds','angles','dihedrals'):
            rows=record['interactions'][kind]
            normalised=[min(tuple(row['atom_ids']),tuple(reversed(row['atom_ids']))) for row in rows]
            self.assertEqual(len(normalised),len(set(normalised)))
            self.assertEqual(len({row['id'] for row in rows}),len(rows))
            for row in rows:
                self.assertEqual(len(row['atom_ids']),len(set(row['atom_ids'])))
                self.assertTrue(all(frozenset(pair) in edges for pair in zip(row['atom_ids'],row['atom_ids'][1:])))
        expected={frozenset(pair) for pair in combinations([a.atom_id for a in geometry.atoms],2)}
        actual={frozenset(pair['atom_ids']) for pair in record['nonbonded_pairs']}
        self.assertEqual(actual,expected)

    def test_atom_ordering_changes_no_chemical_assignment(self):
        original=graph_fixture(4)
        reordered=replace(original,atoms=tuple(reversed(original.atoms)),
                          coordinates_angstrom=tuple(reversed(original.coordinates_angstrom)))
        first,second=(assign_alkane(g).to_dict() for g in (original,reordered))
        self.assertEqual({a['atom_id']:(a['type_id'],a['charge_e']) for a in first['atom_assignments']},
                         {a['atom_id']:(a['type_id'],a['charge_e']) for a in second['atom_assignments']})
        for kind in ('bonds','angles','dihedrals'):
            def canonical(record):
                return {(min(tuple(row['atom_ids']),tuple(reversed(row['atom_ids']))),row['type_id'])
                        for row in record['interactions'][kind]}
            self.assertEqual(canonical(first),canonical(second))

    def test_branches_cycles_and_methane_are_not_silently_typed(self):
        fixtures=[graph_fixture(4,[(0,1),(0,2),(0,3)]),
                  graph_fixture(3,[(0,1),(1,2),(2,0)]),graph_fixture(1)]
        for geometry in fixtures:
            with self.subTest(atoms=len(geometry.atoms)),self.assertRaises(UnsupportedAlkaneError):
                assign_alkane(geometry)

    def test_heteroatoms_charges_isotopes_and_mass_changes_are_rejected(self):
        original=graph_fixture(2)
        for changes in ({'element':'O','mass_g_mol':15.999}, {'formal_charge':1},
                        {'isotope':13,'mass_g_mol':13.00335484},{'mass_g_mol':12.01}):
            with self.subTest(changes=changes),self.assertRaises(UnsupportedAlkaneError):
                assign_alkane(replace_atom(original,0,**changes))

    def test_unsaturated_aromatic_and_ambiguous_bond_tags_are_rejected(self):
        geometry=graph_fixture(2)
        for changes in ({'order':2},{'order':1.5,'aromatic':True},{'stereo':'STEREOANY'}):
            bad=replace(geometry,bonds=(replace(geometry.bonds[0],**changes),)+geometry.bonds[1:])
            with self.subTest(changes=changes),self.assertRaises(UnsupportedAlkaneError):
                assign_alkane(bad)

    def test_missing_explicit_hydrogen_and_hydrogen_bridges_are_rejected(self):
        geometry=graph_fixture(2)
        # Remove a terminal H and its bond; graph stays connected but carbon
        # no longer has its required explicit tetravalent environment.
        victim=geometry.atoms[-1].atom_id
        missing=replace(geometry,atoms=geometry.atoms[:-1],coordinates_angstrom=geometry.coordinates_angstrom[:-1],
                        bonds=tuple(b for b in geometry.bonds if victim not in b.atom_ids))
        with self.assertRaises(UnsupportedAlkaneError):
            assign_alkane(missing)
        # An H bridge is a legal generic connected topology but no alkane H.
        bridge=replace(geometry,bonds=geometry.bonds+(GeometryBond(('H:0:0','carbon:1'),1),))
        with self.assertRaises(UnsupportedAlkaneError):
            assign_alkane(bridge)

    def test_invalid_inputs_and_explicit_pair_budget(self):
        geometry=graph_fixture(2)
        for max_atoms in (True,0,-1,8.,'8'):
            with self.subTest(max_atoms=max_atoms),self.assertRaises(ValueError):
                assign_alkane(geometry,max_atoms=max_atoms)
        with self.assertRaisesRegex(ValueError,'exceeds'):
            assign_alkane(geometry,max_atoms=7)
        with self.assertRaisesRegex(ValueError,'MolecularGeometry'):
            assign_alkane(geometry.to_dict())

    def test_public_serialized_constructor_round_trip_revalidates_frozen_model(self):
        original=assign_alkane(graph_fixture(4))
        restored=AssignedAlkane(json.dumps(original.to_dict(),indent=2))
        self.assertEqual(restored.to_dict(),original.to_dict())
        self.assertEqual(restored.geometry,original.geometry)
        # Assignment depends on chemical topology, not a particular geometry.
        # Coordinate-only edits remain explicit supplied coordinate data.
        moved=original.to_dict()
        moved['geometry']['coordinates_angstrom'][0][0]+=.25
        restored=AssignedAlkane(json.dumps(moved))
        self.assertEqual(restored.to_dict()['geometry'],moved['geometry'])
        self.assertFalse(restored.to_dict()['validation']['production_ready'])

    def test_changed_partial_charges_cannot_keep_the_frozen_model_label(self):
        record=assign_alkane(graph_fixture(4)).to_dict()
        for atom in record['atom_assignments']:
            atom['charge_e']=0.
        with self.assertRaisesRegex(ValueError,'frozen n-alkane model'):
            AssignedAlkane(json.dumps(record))

    def test_missing_proper_torsion_cannot_keep_the_frozen_model_label(self):
        record=assign_alkane(graph_fixture(4)).to_dict()
        record['interactions']['dihedrals'].pop()
        with self.assertRaisesRegex(ValueError,'frozen n-alkane model'):
            AssignedAlkane(json.dumps(record))

    def test_changed_source_parameter_cannot_keep_the_frozen_model_label(self):
        record=assign_alkane(graph_fixture(4)).to_dict()
        record['parameter_set']['parameters']['bonds'][0]['values']['k_kj_mol_angstrom2']+=1.
        with self.assertRaisesRegex(ValueError,'frozen n-alkane model'):
            AssignedAlkane(json.dumps(record))

    def test_changed_one_four_scale_cannot_keep_the_frozen_model_label(self):
        record=assign_alkane(graph_fixture(4)).to_dict()
        pair=next(pair for pair in record['nonbonded_pairs'] if pair['topological_distance']==3)
        pair['lj_scale']=1.
        with self.assertRaisesRegex(ValueError,'frozen n-alkane model'):
            AssignedAlkane(json.dumps(record))

    def test_serialized_constructor_rejects_bad_json_duplicate_and_type_tampering(self):
        for raw in ('not json','[]','{"geometry":{},"geometry":{}}','{"geometry":NaN}'):
            with self.subTest(raw=raw),self.assertRaises(ValueError):
                AssignedAlkane(raw)
        with self.assertRaises(ValueError):
            AssignedAlkane({})
        record=assign_alkane(graph_fixture(2)).to_dict()
        # Python mapping equality would consider these values equal. The
        # canonical frozen-model check preserves their different JSON types.
        record['validation']['formal_charge_e']=False
        with self.assertRaisesRegex(ValueError,'frozen n-alkane model'):
            AssignedAlkane(json.dumps(record))
        with self.assertRaises(ValueError):
            AssignedAlkane(json.dumps(assign_alkane(graph_fixture(2)).to_dict()),max_atoms=True)
        with self.assertRaises(ValueError):
            AssignedAlkane(json.dumps(assign_alkane(graph_fixture(2)).to_dict()),max_atoms=7)


if __name__=='__main__':
    unittest.main()
