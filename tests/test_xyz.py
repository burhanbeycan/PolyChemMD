import copy
import math
import shlex
import unittest

from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry
from polychemmd.io import cell_xyz, molecular_xyz


def molecule():
    return MolecularGeometry(
        (GeometryAtom('oxygen', 1, 'O', 18, -1, 17.999, '{}'),
         GeometryAtom('hydrogen', 2, 'H', 0, 0, 1.008, '{}')),
        (GeometryBond(('oxygen', 'hydrogen'), 1),),
        ((-0.012345678901234567, 0., 2.5), (0.9446543210987654, 0., 2.5)),
        '0'*64, '{"source_kind":"analytic_fixture"}')


def cell():
    return {
        'kind': 'orthorhombic_rigid_packed_cell', 'coordinate_units': 'angstrom',
        'atom_count': 2,
        'box': {'lengths_angstrom': [10, 11, 12], 'periodic': [True]*3,
                'origin_angstrom': [0, 0, 0]},
        'atoms': [
            {'element': 'Cl', 'coordinates_angstrom': [10.25, -0.5, 3.125],
             'wrapped_coordinates_angstrom': [.25, 10.5, 3.125]},
            {'element': 'Na', 'coordinates_angstrom': [-0.5, 11.25, 12.5],
             'wrapped_coordinates_angstrom': [9.5, .25, .5]},
        ],
    }


def parse_xyz(text):
    lines = text.splitlines()
    count = int(lines[0])
    metadata = dict(token.split('=', 1) for token in shlex.split(lines[1]))
    rows = [(parts[0], tuple(map(float, parts[1:])))
            for parts in (line.split() for line in lines[2:])]
    return count, metadata, rows


class XyzExportTests(unittest.TestCase):
    def test_molecule_order_coordinates_precision_and_nonperiodic_header(self):
        geometry = molecule()
        count, header, rows = parse_xyz(molecular_xyz(geometry))
        self.assertEqual(count, 2)
        self.assertEqual([row[0] for row in rows], ['O', 'H'])
        self.assertEqual(tuple(row[1] for row in rows), geometry.coordinates_angstrom)
        self.assertEqual(header['pbc'], 'F F F')
        self.assertEqual(header['Properties'], 'species:S:1:pos:R:3')
        self.assertEqual(header['units'], 'angstrom')
        self.assertEqual(header['companion_JSON'], 'topology isotopes atom_identity')
        self.assertEqual(header['visualization_only'], 'T')
        self.assertNotIn('Lattice', header)

    def test_cell_lattice_and_default_unwrapped_positions(self):
        record = cell()
        before = copy.deepcopy(record)
        count, header, rows = parse_xyz(cell_xyz(record))
        self.assertEqual(count, len(rows))
        self.assertEqual(header['Lattice'], '10 0 0 0 11 0 0 0 12')
        self.assertEqual(header['pbc'], 'T T T')
        self.assertEqual(header['coordinates'], 'unwrapped')
        self.assertEqual(rows, [('Cl', (10.25, -.5, 3.125)), ('Na', (-.5, 11.25, 12.5))])
        self.assertEqual(record, before)

    def test_wrapped_view_selects_coordinates_without_reordering(self):
        count, header, rows = parse_xyz(cell_xyz(cell(), wrapped=True))
        self.assertEqual(header['coordinates'], 'wrapped')
        self.assertEqual(rows, [('Cl', (.25, 10.5, 3.125)), ('Na', (9.5, .25, .5))])

    def test_coordinate_validation_claims_are_not_exported_as_authority(self):
        record = cell()
        record['validation'] = {'scientific_accuracy': 'guaranteed by untrusted input'}
        text = cell_xyz(record)
        self.assertNotIn('guaranteed', text)
        self.assertNotIn('scientific_accuracy', text)

    def test_invalid_counts_units_and_kinds_are_rejected(self):
        for field, value in [('atom_count',True), ('atom_count',0), ('atom_count',3),
                             ('atoms',{}), ('coordinate_units','nm'), ('kind','unknown')]:
            record = cell()
            record[field] = value
            with self.subTest(field=field,value=value), self.assertRaises(ValueError):
                cell_xyz(record)
        for invalid in (None, [], 'cell'):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                cell_xyz(invalid)

    def test_invalid_species_and_nonfinite_coordinate_values_are_rejected(self):
        for species in ('Xx', 'C\nH', '', True, []):
            record = cell()
            record['atoms'][0]['element'] = species
            with self.subTest(species=species), self.assertRaises(ValueError):
                cell_xyz(record)
        for coordinate in (True, float('nan'), float('inf'), 10**1000, '1.0'):
            record = cell()
            record['atoms'][0]['coordinates_angstrom'][0] = coordinate
            with self.subTest(coordinate_type=type(coordinate)), self.assertRaises(ValueError):
                cell_xyz(record)
        record = cell()
        record['atoms'][0]['coordinates_angstrom'] = [1,2]
        with self.assertRaises(ValueError):
            cell_xyz(record)

    def test_invalid_box_axes_origin_lengths_and_wrapped_mode(self):
        for key, value in [('lengths_angstrom',[1,2]), ('lengths_angstrom',[0,2,3]),
                           ('lengths_angstrom',[True,2,3]), ('lengths_angstrom',[math.inf,2,3]),
                           ('periodic',[1,1,1]), ('periodic',[True,True,False]),
                           ('origin_angstrom',[1,0,0])]:
            record = cell()
            record['box'][key] = value
            with self.subTest(key=key,value=value), self.assertRaises(ValueError):
                cell_xyz(record)
        for wrapped in ('true',1,None):
            with self.subTest(wrapped=wrapped), self.assertRaises(ValueError):
                cell_xyz(cell(), wrapped=wrapped)
        record = cell()
        record['atoms'][0]['wrapped_coordinates_angstrom'][0] = 10.
        with self.assertRaises(ValueError):
            cell_xyz(record,wrapped=True)

    def test_non_geometry_cannot_use_molecular_export(self):
        with self.assertRaises(ValueError):
            molecular_xyz(molecule().to_dict())


if __name__ == '__main__':
    unittest.main()
