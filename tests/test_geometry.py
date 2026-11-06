import copy
from dataclasses import FrozenInstanceError
import json
import math
import unittest

from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry


def hydrogen(atom_id, atom_map, isotope=0, mass=1.008):
    return GeometryAtom(atom_id, atom_map, 'H', isotope, 0, mass,
                        json.dumps({'atom_id': atom_id, 'atom_map': atom_map,
                                    'origin': 'analytic_fixture'}))


def dihydrogen():
    # Defined analytic coordinates, not a calculation or a measured structure.
    return MolecularGeometry((hydrogen('a', 1), hydrogen('b', 2, 2, 2.014)),
                             (GeometryBond(('a', 'b'), 1),),
                             ((0., 0., 0.), (0.74, 0., 0.)), '0'*64,
                             '{"kind":"analytic_fixture","method":"specified coordinates"}')


class GeometryContractTests(unittest.TestCase):
    def test_round_trip_preserves_order_isotope_lineage_and_positions(self):
        geometry = dihydrogen()
        restored = MolecularGeometry.from_mapping(geometry.to_dict())
        self.assertEqual(restored, geometry)
        self.assertEqual(restored.atoms[1].isotope, 2)
        self.assertEqual(restored.atoms[1].mass_g_mol, 2.014)
        self.assertEqual(restored.bonds[0].atom_ids, ('a', 'b'))
        self.assertAlmostEqual(math.dist(*restored.coordinates_angstrom), 0.74)

    def test_records_own_all_mutable_input_and_output_data(self):
        data = dihydrogen().to_dict()
        geometry = MolecularGeometry.from_mapping(data)
        data['coordinates_angstrom'][1][0] = 9
        data['atoms'][0]['lineage']['origin'] = 'modified'
        out = geometry.to_dict()
        out['provenance']['method'] = 'modified'
        out['atoms'][0]['lineage']['origin'] = 'modified'
        self.assertEqual(geometry.coordinates_angstrom[1][0], .74)
        self.assertEqual(geometry.atoms[0].lineage['origin'], 'analytic_fixture')
        self.assertEqual(geometry.provenance['method'], 'specified coordinates')
        with self.assertRaises(FrozenInstanceError):
            geometry.source_graph_sha256 = 'f'*64

    def test_constructor_normalizes_list_storage_without_retaining_aliases(self):
        atoms = [hydrogen('a', 1)]
        positions = [[0, 0, 0]]
        geometry = MolecularGeometry(atoms, [], positions, '0'*64, '{}')
        atoms.clear()
        positions[0][0] = 5
        self.assertEqual(len(geometry.atoms), 1)
        self.assertEqual(geometry.coordinates_angstrom, ((0., 0., 0.),))

    def test_duplicate_atom_maps_or_ids_are_rejected(self):
        for atom in (hydrogen('b', 1), hydrogen('a', 2)):
            with self.subTest(atom=atom), self.assertRaisesRegex(ValueError, 'unique'):
                MolecularGeometry((hydrogen('a', 1), atom), (), ((0,0,0), (1,0,0)), '0'*64, '{}')

    def test_missing_duplicate_and_disconnected_bonds_are_rejected(self):
        atoms = dihydrogen().atoms
        variants = [(), (GeometryBond(('a','missing'),1),),
                    (GeometryBond(('a','b'),1),GeometryBond(('b','a'),1))]
        for bonds in variants:
            with self.subTest(bonds=bonds), self.assertRaises(ValueError):
                MolecularGeometry(atoms,bonds,((0,0,0),(1,0,0)),'0'*64,'{}')
        with self.assertRaises(ValueError):
            GeometryBond(('a','a'), 1)

    def test_numerical_poisoning_is_rejected(self):
        for value in (True, float('nan'), float('inf'), -(10**1000), '1.0'):
            with self.subTest(value_type=type(value)), self.assertRaises(ValueError):
                MolecularGeometry((hydrogen('a',1),),(),((value,0,0),),'0'*64,'{}')
        for field in ('atom_map','isotope','formal_charge','mass_g_mol'):
            data = dihydrogen().to_dict()
            data['atoms'][0][field] = True
            with self.subTest(field=field), self.assertRaises(ValueError):
                MolecularGeometry.from_mapping(data)

    def test_invalid_dimensions_units_and_unexpected_fields(self):
        for key, value in [('coordinates_angstrom',[[0,0,0]]),
                           ('coordinates_angstrom',[[0,0],[1,0,0]]),
                           ('coordinate_units','nm'),('schema_version',True),
                           ('source_graph_sha256','not-a-digest')]:
            data = dihydrogen().to_dict()
            data[key] = value
            with self.subTest(key=key,value=value), self.assertRaises(ValueError):
                MolecularGeometry.from_mapping(data)
        data = dihydrogen().to_dict()
        data['unverified_success'] = True
        with self.assertRaises(ValueError):
            MolecularGeometry.from_mapping(data)

    def test_metadata_is_finite_unambiguous_object_json(self):
        for raw in ('[]','{"x":NaN}','{"x":1e400}','{"x":1,"x":2}'):
            with self.subTest(raw=raw), self.assertRaises(ValueError):
                MolecularGeometry((hydrogen('a',1),),(),((0,0,0),),'0'*64,raw)
            with self.subTest(lineage=raw), self.assertRaises(ValueError):
                GeometryAtom('a',1,'H',0,0,1.008,raw)

    def test_lineage_cannot_conflict_with_its_atom_identity(self):
        for raw in ('{"atom_id":"other"}','{"atom_map":true}',
                    '{"element":"C"}','{"formal_charge":1}','{"isotope":2}'):
            with self.subTest(raw=raw), self.assertRaisesRegex(ValueError, 'conflicts'):
                GeometryAtom('a',1,'H',0,0,1.008,raw)

    def test_serialization_does_not_promote_untrusted_provenance(self):
        geometry = MolecularGeometry((hydrogen('a',1),),(),((0,0,0),),'0'*64,
                                     '{"caller_supplied_claim":"physical validation passed"}')
        record = geometry.to_dict()
        self.assertEqual(record['kind'], 'molecular_geometry')
        self.assertNotIn('material_validation', record)
        self.assertNotIn('equilibrium_status', record)
        self.assertIn('caller_supplied_claim', record['provenance'])


if __name__ == '__main__':
    unittest.main()
