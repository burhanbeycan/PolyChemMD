from collections import Counter
from copy import deepcopy
import json
import math
from pathlib import Path
import unittest

from polychemmd.builders.packing import pack_cell
from polychemmd.forcefields.periodic import AssignedPeriodicAlkaneCell, assign_periodic_alkanes
from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry


ROOT = Path(__file__).resolve().parents[1]


def analytic_ethane():
    # Explicit topology/geometry fixture, not a measured or optimized molecule.
    names = ['C0','C1']+[f'H{i}' for i in range(6)]
    atoms = tuple(GeometryAtom(name,i+1,'C' if i < 2 else 'H',0,0,
                              12.011 if i < 2 else 1.008,
                              '{"source_kind":"analytic_geometry_fixture"}')
                  for i,name in enumerate(names))
    points = ((0.,0.,0.),(1.53,0.,0.),(-.363,1.03,0.),
              (-.363,-.515,.892),(-.363,-.515,-.892),
              (1.893,.515,.892),(1.893,.515,-.892),(1.893,-1.03,0.))
    pairs = [('C0','C1')]+[('C0',f'H{i}') for i in range(3)]+[('C1',f'H{i}') for i in range(3,6)]
    return MolecularGeometry(atoms,tuple(GeometryBond(pair,1.) for pair in pairs),points,
                             '0'*64,'{"kind":"analytic_geometry_fixture"}')


def packed_ethanes():
    return pack_cell([analytic_ethane(),analytic_ethane()],box_lengths_angstrom=(24.,25.,26.),
                     seed=353,min_distance_angstrom=1.5).to_dict()


def translate_molecule(cell, molecule_index, delta):
    moved=deepcopy(cell); lengths=moved['box']['lengths_angstrom']
    placement=moved['molecules'][molecule_index]['placement']
    placement['translated_centroid_angstrom']=[v+d for v,d in zip(placement['translated_centroid_angstrom'],delta)]
    for atom in moved['atoms']:
        if atom['molecule_index'] != molecule_index:
            continue
        point=[v+d for v,d in zip(atom['coordinates_angstrom'],delta)]
        images=[math.floor(v/l) for v,l in zip(point,lengths)]
        atom['coordinates_angstrom']=point
        atom['wrapped_coordinates_angstrom']=[v-n*l for v,n,l in zip(point,images,lengths)]
        atom['image_flags']=images
    return moved


class PeriodicAssignmentTests(unittest.TestCase):
    def test_two_ethanes_have_complete_local_topology_and_cross_molecule_pairs(self):
        record=assign_periodic_alkanes(packed_ethanes()).to_dict()
        self.assertEqual(len(record['atoms']),16)
        self.assertEqual(Counter(a['type_id'] for a in record['atoms']),{'opls_135':4,'opls_140':12})
        self.assertEqual({k:len(v) for k,v in record['interactions'].items()},
                         {'bonds':14,'angles':24,'dihedrals':18,'impropers':0})
        self.assertEqual(len(record['nonbonded_pairs']),120)
        pairs=record['nonbonded_pairs']
        self.assertEqual(sum(not p['same_molecule'] for p in pairs),64)
        self.assertEqual(sum(p['topological_distance']==3 for p in pairs),18)
        for pair in pairs:
            if not pair['same_molecule']:
                self.assertIsNone(pair['topological_distance'])
                self.assertEqual((pair['lj_scale'],pair['coulomb_scale']),(1.,1.))
            elif pair['topological_distance'] <= 2:
                self.assertEqual((pair['lj_scale'],pair['coulomb_scale']),(0.,0.))
            else:
                self.assertEqual((pair['lj_scale'],pair['coulomb_scale']),(.5,.5))

    def test_source_and_global_identity_remain_distinct_and_owned(self):
        source=packed_ethanes(); before=deepcopy(source)
        assigned=assign_periodic_alkanes(source); record=assigned.to_dict()
        self.assertEqual(source,before)
        self.assertEqual(record['source_cell'],source)
        for old,new in zip(source['atoms'],record['atoms']):
            self.assertEqual({k:v for k,v in new.items() if k not in ('type_id','charge_e')},old)
        self.assertEqual(len({a['cell_atom_id'] for a in record['atoms']}),16)
        self.assertEqual(len({a['source_atom_id'] for a in record['atoms']}),8)
        source['atoms'][0]['lineage']['changed']=True
        record['atoms'][0]['charge_e']=99
        self.assertEqual(assigned.to_dict()['source_cell'],before)
        self.assertNotEqual(assigned.to_dict()['atoms'][0]['charge_e'],99)

    def test_no_bonded_interaction_crosses_molecules(self):
        record=assign_periodic_alkanes(packed_ethanes()).to_dict()
        lookup={a['cell_atom_id']:a['molecule_index'] for a in record['atoms']}
        for kind in ('bonds','angles','dihedrals'):
            self.assertEqual(len({r['id'] for r in record['interactions'][kind]}),len(record['interactions'][kind]))
            for row in record['interactions'][kind]:
                self.assertEqual(len({lookup[i] for i in row['atom_ids']}),1)
                self.assertEqual(len(set(row['atom_ids'])),len(row['atom_ids']))

    def test_actual_four_decane_starting_cell_with_cli_extra_metadata(self):
        source=json.loads((ROOT/'evidence'/'geometry_v03'/'cell.json').read_text())['payload']
        record=assign_periodic_alkanes(source).to_dict()
        self.assertIn('preparation',record['source_cell'])
        self.assertEqual(len(record['atoms']),128)
        self.assertEqual([len(record['interactions'][k]) for k in ('bonds','angles','dihedrals')],[124,240,324])
        self.assertEqual(sum(not p['same_molecule'] for p in record['nonbonded_pairs']),6144)
        self.assertAlmostEqual(record['validation']['total_partial_charge_e'],0.,places=14)
        self.assertEqual(record['validation']['total_molar_mass_g_mol'],source['total_molar_mass_g_mol'])
        self.assertFalse(record['validation']['production_ready'])
        self.assertEqual(record['validation']['cutoff_and_long_range_model'],'not_selected_by_assignment')

    def test_boundary_crossing_and_whole_lattice_translation_are_preserved(self):
        source=packed_ethanes()
        centroid=source['molecules'][0]['placement']['translated_centroid_angstrom']
        moved=translate_molecule(source,0,[-centroid[0]+.1,0,0])
        self.assertTrue(any(a['image_flags'][0] != 0 for a in moved['atoms'] if a['molecule_index']==0))
        first=assign_periodic_alkanes(moved).to_dict()
        shifted=translate_molecule(moved,0,[24.,-25.,52.])
        second=assign_periodic_alkanes(shifted).to_dict()
        self.assertEqual(first['nonbonded_pairs'],second['nonbonded_pairs'])
        self.assertEqual(first['interactions'],second['interactions'])
        for old,new in zip(first['atoms'],second['atoms']):
            if old['molecule_index']==0:
                for a,b in zip(old['wrapped_coordinates_angstrom'],new['wrapped_coordinates_angstrom']):
                    self.assertAlmostEqual(a,b,places=12)

    def test_changed_wrapped_image_or_rigid_coordinates_are_rejected(self):
        source=packed_ethanes()
        for field in ('coordinates_angstrom','wrapped_coordinates_angstrom','image_flags'):
            changed=deepcopy(source)
            changed['atoms'][0][field][0]+=.1 if field!='image_flags' else 1
            with self.subTest(field=field),self.assertRaises(ValueError):
                assign_periodic_alkanes(changed)
        changed=deepcopy(source)
        for atom in changed['atoms']:
            if atom['molecule_index']==0:
                atom['coordinates_angstrom'][0]+=.1
                atom['wrapped_coordinates_angstrom'][0]+=.1
        with self.assertRaisesRegex(ValueError,'rigid coordinates'):
            assign_periodic_alkanes(changed)

    def test_improper_rotation_or_false_source_centroid_is_rejected(self):
        for field in ('reflection','centroid'):
            source=packed_ethanes()
            if field=='reflection':
                source['molecules'][0]['placement']['rotation_matrix'][0]=[-x for x in source['molecules'][0]['placement']['rotation_matrix'][0]]
            else:
                source['molecules'][0]['placement']['source_centroid_angstrom'][0]+=.2
            with self.subTest(field=field),self.assertRaises(ValueError):
                assign_periodic_alkanes(source)

    def test_unrepresentable_large_translation_cannot_collapse_a_molecule(self):
        source=packed_ethanes()
        moved=translate_molecule(source,0,[24.*2**54,0.,0.])
        with self.assertRaises(ValueError):
            assign_periodic_alkanes(moved)

    def test_atom_id_map_molecule_lineage_mass_and_charge_tampering_fail(self):
        variants=[('cell_atom_id','other'),('atom_map',2),('source_atom_map',99),
                  ('molecule_index',1),('lineage',{'fake':True}),('mass_g_mol',12.0),
                  ('formal_charge',1),('isotope',13),('source_atom_id','missing')]
        for key,value in variants:
            source=packed_ethanes();source['atoms'][0][key]=value
            with self.subTest(key=key),self.assertRaises(ValueError):
                assign_periodic_alkanes(source)
        source=packed_ethanes();source['atoms'][0]['charge_e']=0
        with self.assertRaises(ValueError):
            assign_periodic_alkanes(source)

    def test_wrong_missing_and_cross_molecule_source_bonds_fail(self):
        for change in ('order','map','missing','cross'):
            source=packed_ethanes()
            if change=='order':source['bonds'][0]['order']=2
            elif change=='map':source['bonds'][0]['atom_maps'][0]=999
            elif change=='missing':
                source['bonds'].pop();source['bond_count']-=1
            else:
                source['bonds'][0]['cell_atom_ids'][1]=source['atoms'][-1]['cell_atom_id']
                source['bonds'][0]['atom_maps'][1]=source['atoms'][-1]['atom_map']
            with self.subTest(change=change),self.assertRaises(ValueError):
                assign_periodic_alkanes(source)

    def test_schema_box_counts_and_nonfinite_records_fail(self):
        source=packed_ethanes()
        mutations=[lambda s:s.update(coordinate_units='nm'),
                   lambda s:s['box'].update(periodic=[1,1,1]),
                   lambda s:s['box'].update(origin_angstrom=[1,0,0]),
                   lambda s:s['box'].update(volume_angstrom3=1),
                   lambda s:s.update(atom_count=True),
                   lambda s:s['atoms'][0]['image_flags'].__setitem__(0,True),
                   lambda s:s['atoms'][0]['coordinates_angstrom'].__setitem__(0,float('nan')),
                   lambda s:s.update(total_molar_mass_g_mol=1),
                   lambda s:s.update(construction_density_g_cm3=1)]
        for mutate in mutations:
            changed=deepcopy(source);mutate(changed)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError):
                assign_periodic_alkanes(changed)

    def test_source_validation_flags_are_not_promoted_to_new_claims(self):
        source=packed_ethanes()
        source['validation']={'physical_accuracy':'perfect','equilibrium':True}
        source['material_validation']='experimentally_validated'
        record=assign_periodic_alkanes(source).to_dict()
        self.assertEqual(record['source_cell']['validation'],source['validation'])
        self.assertEqual(record['validation']['equilibrium_status'],'not_assessed')
        self.assertEqual(record['validation']['material_property_validation'],'not_performed')
        self.assertFalse(record['validation']['simulation_performed'])

    def test_public_constructor_rejects_changed_model_and_preserves_exact_round_trip(self):
        record=assign_periodic_alkanes(packed_ethanes()).to_dict()
        self.assertEqual(AssignedPeriodicAlkaneCell(json.dumps(record)).to_dict(),record)
        for kind in ('charge','coefficient','pair','topology','claim'):
            changed=deepcopy(record)
            if kind=='charge':changed['atoms'][0]['charge_e']=0
            elif kind=='coefficient':changed['parameter_set']['parameters']['bonds'][0]['values']['k_kj_mol_angstrom2']+=1
            elif kind=='pair':changed['nonbonded_pairs'][-1]['lj_scale']=0
            elif kind=='topology':changed['interactions']['dihedrals'].pop()
            else:changed['validation']['production_ready']=True
            with self.subTest(kind=kind),self.assertRaises(ValueError):
                AssignedPeriodicAlkaneCell(json.dumps(changed))

    def test_constructor_json_and_budget_gates(self):
        for raw in ('[]','{"x":NaN}','{"source_cell":{},"source_cell":{}}'):
            with self.subTest(raw=raw),self.assertRaises(ValueError):
                AssignedPeriodicAlkaneCell(raw)
        source=packed_ethanes()
        for limit in (True,0,-1,15,16.):
            with self.subTest(limit=limit),self.assertRaises(ValueError):
                assign_periodic_alkanes(source,max_atoms=limit)
        self.assertEqual(len(assign_periodic_alkanes(source,max_atoms=16).to_dict()['atoms']),16)


if __name__=='__main__':
    unittest.main()
