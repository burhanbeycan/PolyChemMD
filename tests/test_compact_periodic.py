from collections import Counter, deque
from copy import deepcopy
import hashlib
from itertools import combinations
import json
import math
from pathlib import Path
import unittest
from unittest.mock import patch

from polychemmd.analysis.density import density_from_mass_volume
from polychemmd.builders.packing import pack_cell
from polychemmd.forcefields import periodic
from polychemmd.forcefields.periodic import (
    AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell,
    assign_periodic_alkanes, assign_compact_periodic_alkanes,
)
from test_periodic_assignment import analytic_ethane, packed_ethanes, translate_molecule


ROOT = Path(__file__).resolve().parents[1]


def replicated_ethanes(count):
    """Deterministic sparse geometry fixture, not an equilibrium material."""
    source = pack_cell([analytic_ethane()], box_lengths_angstrom=(12.,12.,12.),
                       seed=532,min_distance_angstrom=1.5).to_dict()
    result = deepcopy(source)
    lengths = [12.*count,12.,12.]
    result['box']['lengths_angstrom'] = lengths
    result['box']['volume_angstrom3'] = math.prod(lengths)
    result['atoms'],result['bonds'],result['molecules'] = [],[],[]
    for index in range(count):
        molecule = deepcopy(source['molecules'][0]);molecule['molecule_index'] = index
        molecule['placement']['translated_centroid_angstrom'][0] += 12.*index
        result['molecules'].append(molecule)
        for old in source['atoms']:
            atom=deepcopy(old)
            atom['molecule_index']=index
            atom['cell_atom_id']=f"molecule:{index}/{old['source_atom_id']}"
            atom['atom_map']=8*index+old['atom_map']
            atom['coordinates_angstrom'][0]+=12.*index
            atom['image_flags']=[math.floor(v/l) for v,l in zip(atom['coordinates_angstrom'],lengths)]
            atom['wrapped_coordinates_angstrom']=[v-n*l for v,n,l in
                zip(atom['coordinates_angstrom'],atom['image_flags'],lengths)]
            result['atoms'].append(atom)
        for old in source['bonds']:
            bond=deepcopy(old);bond['molecule_index']=index
            bond['cell_atom_ids']=[f"molecule:{index}/{name.split('/',1)[1]}" for name in old['cell_atom_ids']]
            bond['atom_maps']=[8*index+value for value in old['atom_maps']]
            result['bonds'].append(bond)
    result.update(atom_count=8*count,molecule_count=count,bond_count=7*count,
                  total_molar_mass_g_mol=source['total_molar_mass_g_mol']*count,
                  construction_density_g_cm3=density_from_mass_volume(
                      [atom['mass_g_mol'] for atom in result['atoms']],math.prod(lengths)))
    result['packing']={'method':'deterministic_sparse_analytic_fixture_replication',
                       'independent_cells':1,'equilibrium_claim':False}
    return result


def graph_pair_audit(source):
    """Independent all-pair BFS oracle; deliberately used only for small cases."""
    atoms=source['atoms'];adjacent={atom['cell_atom_id']:set() for atom in atoms}
    for bond in source['bonds']:
        first,second=bond['cell_atom_ids'];adjacent[first].add(second);adjacent[second].add(first)
    distances={}
    for atom in atoms:
        identity=atom['cell_atom_id'];current={identity:0};queue=deque([identity])
        while queue:
            first=queue.popleft()
            for second in adjacent[first]:
                if second not in current:current[second]=current[first]+1;queue.append(second)
        distances[identity]=current
    pairs=[]
    for first,second in combinations(atoms,2):
        ids=[first['cell_atom_id'],second['cell_atom_id']]
        distance=distances[ids[0]].get(ids[1])
        scale=0. if distance in (1,2) else (.5 if distance==3 else 1.)
        pairs.append({'atom_ids':ids,'same_molecule':distance is not None,
                      'topological_distance':distance,'lj_scale':scale,'coulomb_scale':scale})
    return pairs


class CompactPeriodicAssignmentTests(unittest.TestCase):
    def test_small_cell_matches_independent_all_pair_graph_oracle(self):
        source=packed_ethanes();legacy=assign_periodic_alkanes(source).to_dict()
        compact=assign_compact_periodic_alkanes(source).to_dict();expected=graph_pair_audit(source)
        self.assertEqual(legacy['nonbonded_pairs'],expected)
        exceptions={tuple(pair['atom_ids']):pair for pair in compact['special_pairs']}
        self.assertEqual(compact['special_pairs'],[row for row in expected if row['lj_scale']!=1])
        for pair in expected:
            actual=exceptions.get(tuple(pair['atom_ids']),compact['nonbonded_rules']['default'])
            self.assertEqual((actual['lj_scale'],actual['coulomb_scale']),
                             (pair['lj_scale'],pair['coulomb_scale']))
        self.assertNotIn('nonbonded_pairs',compact)
        self.assertEqual(compact['validation']['pair_scaling_counts'],legacy['validation']['pair_scaling_counts'])
        self.assertEqual(compact['validation']['total_unordered_pair_count'],120)
        for field in ('atoms','interactions','parameter_set','source_cell','source_snapshot','model_id','box'):
            self.assertEqual(compact[field],legacy[field])

    def test_four_decane_nonlocal_pairs_have_implicit_unit_scale(self):
        source=json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        compact=assign_compact_periodic_alkanes(source).to_dict()
        expected=graph_pair_audit(source)
        self.assertEqual(compact['special_pairs'],[row for row in expected if row['lj_scale']!=1])
        self.assertEqual(len(compact['special_pairs']),688)
        self.assertEqual(compact['validation']['pair_scaling_counts'],
            {'intermolecular':6144,'excluded':364,'scaled_14':324,'unscaled_intramolecular':1296})
        self.assertTrue(all(row['same_molecule'] and row['topological_distance']<=3
                            for row in compact['special_pairs']))

    def test_global_all_pair_enumeration_is_absent_beyond_legacy_capacity(self):
        source=replicated_ethanes(80)
        original=periodic.combinations
        def guarded(sequence,size):
            self.assertLessEqual(len(sequence),8,'global quadratic enumeration reappeared')
            return original(sequence,size)
        with patch.object(periodic,'combinations',side_effect=guarded):
            compact=assign_compact_periodic_alkanes(source).to_dict()
        self.assertEqual(compact['validation']['atom_count'],640)
        self.assertEqual(len(compact['special_pairs']),80*28)
        self.assertEqual(compact['validation']['intramolecular_pair_work'],80*28)
        self.assertEqual(compact['validation']['total_unordered_pair_count'],640*639//2)
        self.assertEqual(len(compact['molecules']),80)
        self.assertTrue(all(set(row)=={'molecule_index','atom_count','source_assignment_sha256'}
                            for row in compact['molecules']))
        with self.assertRaisesRegex(ValueError,'max_atoms'):assign_periodic_alkanes(source)

    def test_sparse_record_growth_tracks_molecule_count_not_global_pairs(self):
        small=assign_compact_periodic_alkanes(replicated_ethanes(8)).to_dict()
        large=assign_compact_periodic_alkanes(replicated_ethanes(80)).to_dict()
        self.assertEqual(len(large['special_pairs']),10*len(small['special_pairs']))
        self.assertLess(len(json.dumps(large)),12*len(json.dumps(small)))
        self.assertGreater(large['validation']['total_unordered_pair_count'],
                           90*small['validation']['total_unordered_pair_count'])

    def test_molecular_assignment_digest_is_of_recomputed_complete_local_model(self):
        source=packed_ethanes();legacy=assign_periodic_alkanes(source).to_dict()
        compact=assign_compact_periodic_alkanes(source).to_dict()
        for old,new in zip(legacy['molecules'],compact['molecules']):
            digest=hashlib.sha256(json.dumps(old['source_assignment'],sort_keys=True,
                separators=(',',':'),allow_nan=False).encode()).hexdigest()
            self.assertEqual(new['source_assignment_sha256'],digest)
            self.assertEqual(new['atom_count'],len(old['source_assignment']['geometry']['atoms']))

    def test_round_trip_and_input_output_ownership(self):
        source=packed_ethanes();before=deepcopy(source)
        assigned=assign_compact_periodic_alkanes(source);record=assigned.to_dict()
        self.assertEqual(CompactAssignedPeriodicAlkaneCell(json.dumps(record)).to_dict(),record)
        source['atoms'][0]['coordinates_angstrom'][0]+=1
        record['nonbonded_rules']['default']['lj_scale']=0
        self.assertEqual(assigned.to_dict()['source_cell'],before)
        self.assertEqual(assigned.to_dict()['nonbonded_rules']['default']['lj_scale'],1.)

    def test_modified_rules_pairs_topology_model_or_validation_are_rejected(self):
        record=assign_compact_periodic_alkanes(packed_ethanes()).to_dict()
        changes=[lambda d:d['nonbonded_rules']['default'].update(lj_scale=.5),
                 lambda d:d['nonbonded_rules']['one_four'].update(coulomb_scale=1.),
                 lambda d:d['special_pairs'][0].update(lj_scale=1.),
                 lambda d:d['special_pairs'].pop(),
                 lambda d:d['special_pairs'].append(deepcopy(d['special_pairs'][0])),
                 lambda d:d['special_pairs'][0].update(topological_distance=99),
                 lambda d:d['atoms'][0].update(charge_e=0.),
                 lambda d:d['interactions']['bonds'].pop(),
                 lambda d:d['parameter_set']['parameters']['bonds'][0]['values'].update(k_kj_mol_angstrom2=99.),
                 lambda d:d['validation']['pair_scaling_counts'].update(intermolecular=0),
                 lambda d:d['validation'].update(production_ready=True),
                 lambda d:d['molecules'][0].update(source_assignment_sha256='0'*64),
                 lambda d:d.update(schema_version='1'),lambda d:d.update(model_id='other'),
                 lambda d:d.update(nonbonded_pairs=[])]
        for change in changes:
            data=deepcopy(record);change(data)
            with self.subTest(change=change),self.assertRaises(ValueError):
                CompactAssignedPeriodicAlkaneCell(json.dumps(data))

    def test_sparse_assignment_retains_boundary_crossing_images(self):
        source=packed_ethanes();center=source['molecules'][0]['placement']['translated_centroid_angstrom']
        moved=translate_molecule(source,0,[-center[0]+.1,0,0])
        original=assign_compact_periodic_alkanes(source).to_dict()
        translated=assign_compact_periodic_alkanes(moved).to_dict()
        self.assertTrue(any(a['image_flags'][0]!=0 for a in translated['atoms']))
        self.assertEqual(translated['special_pairs'],original['special_pairs'])
        self.assertEqual(translated['atoms'][0]['image_flags'],moved['atoms'][0]['image_flags'])

    def test_source_geometry_maps_lineage_and_mass_are_still_revalidated(self):
        mutations=[lambda d:d['atoms'][0]['coordinates_angstrom'].__setitem__(0,99.),
                   lambda d:d['atoms'][0].update(atom_map=2),
                   lambda d:d['atoms'][0].update(lineage={'fake':True}),
                   lambda d:d['atoms'][0].update(mass_g_mol=1.)]
        for mutate in mutations:
            source=packed_ethanes();mutate(source)
            with self.subTest(mutate=mutate),self.assertRaises(ValueError):
                assign_compact_periodic_alkanes(source)

    def test_atom_local_and_aggregate_limits_fail_before_assignment(self):
        source=replicated_ethanes(10)
        for options in ({'max_atoms':79},{'max_atoms_per_molecule':7},
                        {'max_intramolecular_pairs':279}):
            with self.subTest(options=options),patch.object(periodic,'assign_alkane') as assign:
                with self.assertRaises(ValueError):assign_compact_periodic_alkanes(source,**options)
                assign.assert_not_called()
        self.assertEqual(assign_compact_periodic_alkanes(source,max_atoms=80,
            max_atoms_per_molecule=8,max_intramolecular_pairs=280).to_dict()['validation']['atom_count'],80)

    def test_budget_numbers_boolean_and_unbounded_expansion_are_rejected(self):
        source=packed_ethanes()
        for key,ceiling in (('max_atoms',8192),('max_atoms_per_molecule',512),
                            ('max_intramolecular_pairs',2_000_000)):
            for value in (True,0,-1,1.5,float('nan'),float('inf'),ceiling+1):
                with self.subTest(key=key,value=value),self.assertRaises(ValueError):
                    assign_compact_periodic_alkanes(source,**{key:value})

    def test_legacy_four_decane_canonical_assignment_is_byte_compatible(self):
        source=json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        record=assign_periodic_alkanes(source).to_dict()
        raw=json.dumps(record,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        self.assertEqual(hashlib.sha256(raw).hexdigest(),
                         '07cf4aff9366ad00b4da1972ded417c7cd283a5a8037e9cfd129218a0bfaa307')
        self.assertEqual(AssignedPeriodicAlkaneCell(raw.decode()).to_dict(),record)

    def test_unsupported_json_and_duplicate_fields_are_explicit(self):
        for raw in (True,'[]','{"source_cell":null,"source_cell":{}}'):
            with self.subTest(raw=raw),self.assertRaises(ValueError):CompactAssignedPeriodicAlkaneCell(raw)


if __name__=='__main__':unittest.main()
