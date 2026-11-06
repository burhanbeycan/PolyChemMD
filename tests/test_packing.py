"""Independent geometry examples; no force-field or material benchmark data."""

import itertools
import json
import math
import random
import unittest
from unittest.mock import patch

from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry
from polychemmd.builders.packing import (
    PackingError, _fits_self_image_bound, _minimum_image_distance,
    _rotate, _rotation_matrix, _wrap, pack_cell,
)


def geometry(points=((0.0, 0.0, 0.0),), *, mass=4.0026):
    atoms = tuple(GeometryAtom(f"source:{i}", i+1, "He", 0, 0, mass,
                               json.dumps({"atom_id": f"source:{i}", "atom_map": i+1,
                                           "origin": "analytic_geometry_fixture"}))
                  for i in range(len(points)))
    bonds = tuple(GeometryBond((atoms[0].atom_id, atom.atom_id), 1.0) for atom in atoms[1:])
    return MolecularGeometry(atoms, bonds, tuple(points), "0"*64,
                             '{"source_kind":"analytic_geometry_fixture","physical_validation":false}')


def tetrahedron():
    return geometry(((0.0, 0.0, 0.0), (1.5, 0.0, 0.0),
                     (0.0, 1.5, 0.0), (0.0, 0.0, 1.5)))


def determinant(rows):
    a, b, c = rows
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


def signed_volume(points):
    origin = points[0]
    return determinant(tuple(tuple(a-b for a, b in zip(p, origin)) for p in points[1:4]))


class _IdentityAndFixedTranslations:
    def __init__(self, translations):
        self.translations = iter(translations)
        self.gaussian_index = 0

    def gauss(self, mean, sigma):
        result = 1.0 if self.gaussian_index % 4 == 0 else 0.0
        self.gaussian_index += 1
        return result

    def random(self):
        return next(self.translations)


class PeriodicGeometryTests(unittest.TestCase):
    def test_minimum_image_distances_across_one_and_three_faces(self):
        self.assertAlmostEqual(_minimum_image_distance((.1, 5, 5), (9.9, 5, 5), (10, 10, 10)), .2)
        self.assertAlmostEqual(_minimum_image_distance((.1, .2, .3), (9.9, 9.8, 9.7), (10, 10, 10)),
                               math.sqrt(.2**2+.4**2+.6**2))
        self.assertAlmostEqual(_minimum_image_distance((1, 2, 3), (4, 6, 3), (20, 20, 20)), 5.0)

    def test_wrapping_keeps_image_flags_and_half_open_box(self):
        point = (-.25, 10.0, 21.5)
        box = (10.0, 10.0, 10.0)
        wrapped, flags = _wrap(point, box)
        self.assertEqual(wrapped, (9.75, 0.0, 1.5))
        self.assertEqual(flags, (-1, 1, 2))
        self.assertEqual(tuple(x+n*l for x, n, l in zip(wrapped, flags, box)), point)

    def test_extent_bound_excludes_nonzero_periodic_self_images(self):
        box, cutoff = (10, 10, 10), 1.0
        points = ((-4.5, 0, 0), (4.5, 0, 0))
        self.assertTrue(_fits_self_image_bound(points, box, cutoff))
        self.assertFalse(_fits_self_image_bound(((-4.51, 0, 0), (4.5, 0, 0)), box, cutoff))
        for image in itertools.product((-1, 0, 1), repeat=3):
            if image == (0, 0, 0):
                continue
            for first in points:
                for second in points:
                    shifted = tuple(v+i*l for v, i, l in zip(second, image, box))
                    self.assertGreaterEqual(math.dist(first, shifted), cutoff)

    def test_rotation_is_proper_and_preserves_chiral_orientation(self):
        rng = random.Random(127)
        points = tetrahedron().coordinates_angstrom
        for _ in range(50):
            matrix = _rotation_matrix(rng)
            self.assertAlmostEqual(determinant(matrix), 1.0, places=12)
            for i, row in enumerate(matrix):
                for j, other in enumerate(matrix):
                    self.assertAlmostEqual(sum(a*b for a, b in zip(row, other)), float(i == j), places=12)
            transformed = tuple(_rotate(p, matrix) for p in points)
            self.assertAlmostEqual(signed_volume(transformed), signed_volume(points), places=12)


class CellPackingTests(unittest.TestCase):
    def pack(self, geometries=None, **kwargs):
        options = {"box_lengths_angstrom": (8, 9, 10), "seed": 991,
                   "min_distance_angstrom": 1.4}
        options.update(kwargs)
        return pack_cell(geometries or (tetrahedron(),)*4, **options).to_dict()

    def test_complete_independent_pair_and_image_audit(self):
        result = self.pack()
        box = result['box']['lengths_angstrom']
        minimum = math.inf
        for first, second in itertools.combinations(result['atoms'], 2):
            if first['molecule_index'] == second['molecule_index']:
                continue
            # Independent calculation directly from the 27 neighbour images.
            a, b = first['wrapped_coordinates_angstrom'], second['wrapped_coordinates_angstrom']
            distance = min(math.dist(a, tuple(v+i*l for v, i, l in zip(b, image, box)))
                           for image in itertools.product((-1, 0, 1), repeat=3))
            self.assertGreaterEqual(distance, 1.4)
            minimum = min(minimum, distance)
        self.assertAlmostEqual(result['packing']['minimum_inter_molecular_distance_angstrom'], minimum)
        self.assertTrue(result['validation']['intermolecular_minimum_image_audit'])

    def test_rigid_molecule_distances_and_handedness_preserved(self):
        source = tetrahedron()
        result = self.pack((source,)*3)
        for index in range(3):
            points = tuple(a['coordinates_angstrom'] for a in result['atoms'] if a['molecule_index'] == index)
            for i, j in itertools.combinations(range(4), 2):
                self.assertAlmostEqual(math.dist(points[i], points[j]),
                                       math.dist(source.coordinates_angstrom[i], source.coordinates_angstrom[j]), places=12)
            self.assertAlmostEqual(signed_volume(points), signed_volume(source.coordinates_angstrom), places=12)

    def test_global_identity_is_separate_from_source_lineage(self):
        source = tetrahedron()
        result = self.pack((source, source))
        self.assertEqual(len({a['cell_atom_id'] for a in result['atoms']}), 8)
        self.assertEqual([a['atom_map'] for a in result['atoms']], list(range(1, 9)))
        for record in result['atoms']:
            original = source.atoms[record['source_atom_map']-1]
            self.assertEqual(record['source_atom_id'], original.atom_id)
            self.assertEqual(record['lineage'], json.loads(original.lineage_json))
            self.assertEqual(record['lineage']['atom_map'], record['source_atom_map'])
        for molecule in result['molecules']:
            self.assertEqual(molecule['source_geometry'], source.to_dict())
        global_ids = {a['cell_atom_id'] for a in result['atoms']}
        self.assertEqual(len(result['bonds']), 6)
        self.assertTrue(all(set(bond['cell_atom_ids']).issubset(global_ids) for bond in result['bonds']))

    def test_wrapped_coordinates_reconstruct_unwrapped_result(self):
        result = self.pack()
        box = result['box']['lengths_angstrom']
        for atom in result['atoms']:
            for unwrapped, wrapped, image, length in zip(atom['coordinates_angstrom'], atom['wrapped_coordinates_angstrom'], atom['image_flags'], box):
                self.assertGreaterEqual(wrapped, 0)
                self.assertLess(wrapped, length)
                self.assertIs(type(image), int)
                self.assertAlmostEqual(wrapped + image*length, unwrapped, places=12)

    def test_seed_reproducibility_and_global_rng_isolation(self):
        state = random.getstate()
        first, second = self.pack(), self.pack()
        self.assertEqual(random.getstate(), state)
        self.assertEqual(first, second)
        self.assertNotEqual(first['atoms'], self.pack(seed=992)['atoms'])

    def test_density_is_constructed_mass_over_prescribed_volume(self):
        result = self.pack((geometry(mass=4.0026),)*2, box_lengths_angstrom=(10, 10, 10))
        expected = 2 * 4.0026 / (6.02214076e23 * 1000 * 1e-24)
        self.assertAlmostEqual(result['construction_density_g_cm3'], expected, places=15)
        self.assertIn('constructed, not predicted', result['density_interpretation'])
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertEqual(result['phase'], 'not_assigned')
        self.assertFalse(result['simulation_performed'])
        self.assertFalse(result['validation']['force_field_assigned'])

    def test_periodic_contact_is_rejected_even_when_unwrapped_points_are_far(self):
        fake = _IdentityAndFixedTranslations((.01, .5, .5, .99, .5, .5))
        with patch('polychemmd.builders.packing.random.Random', return_value=fake):
            with self.assertRaisesRegex(PackingError, 'contact rejections=1'):
                self.pack((geometry(),)*2, box_lengths_angstrom=(10, 10, 10),
                          min_distance_angstrom=.5, max_attempts_per_molecule=1)

    def test_single_molecule_self_images_and_impossible_extent_fail_explicitly(self):
        with self.assertRaisesRegex(PackingError, 'own periodic image'):
            self.pack((geometry(),), box_lengths_angstrom=(1, 2, 2))
        long_geometry = geometry(((0.0, 0.0, 0.0), (20.0, 0.0, 0.0)))
        with self.assertRaisesRegex(PackingError, 'self-image extent rejections=3'):
            self.pack((long_geometry,), box_lengths_angstrom=(10, 10, 10),
                      min_distance_angstrom=1, max_attempts_per_molecule=3)

    def test_budget_includes_final_independent_audit_and_returns_no_partial_result(self):
        # Two one-atom molecules need one candidate check and one final check.
        with self.assertRaisesRegex(PackingError, 'max_distance_checks exhausted'):
            self.pack((geometry(),)*2, max_distance_checks=1)
        one = self.pack((geometry(),), max_distance_checks=1)
        self.assertIsNone(one['packing']['minimum_inter_molecular_distance_angstrom'])
        self.assertEqual(one['packing']['distance_checks_used'], 0)

    def test_invalid_inputs_and_numeric_ranges(self):
        for kwargs in ({'seed': True}, {'seed': -1}, {'max_atoms': 3},
                       {'max_attempts_per_molecule': 0}, {'max_distance_checks': False},
                       {'min_distance_angstrom': True}, {'min_distance_angstrom': float('nan')},
                       {'box_lengths_angstrom': (1, 2)}, {'box_lengths_angstrom': '111'},
                       {'box_lengths_angstrom': (True, 2, 3)},
                       {'box_lengths_angstrom': (1e-300,)*3, 'min_distance_angstrom':1e-301},
                       {'box_lengths_angstrom': (1e200,)*3}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                self.pack(**kwargs)
        for geometries in ((), 'bad', [object()], None, iter((geometry(),))):
            with self.subTest(geometries=geometries), self.assertRaises(ValueError):
                pack_cell(geometries, box_lengths_angstrom=(10, 10, 10), seed=0)
        with self.assertRaisesRegex(ValueError, 'molar mass'):
            self.pack((geometry(mass=1e308),)*2)
        with self.assertRaisesRegex(ValueError, 'coincident'):
            self.pack((geometry(((0, 0, 0), (0, 0, 0))),))

    def test_isotope_charge_lineage_and_defensive_serialization(self):
        atom = GeometryAtom('D:0', 9, 'H', 2, 1, 2.01410177812,
                            '{"atom_id":"D:0","atom_map":9,"isotope":2,"origin":"analytic_fixture"}')
        source = MolecularGeometry((atom,), (), ((0, 0, 0),), '1'*64, '{"source_kind":"analytic_fixture"}')
        record = pack_cell((source, source), box_lengths_angstrom=(10, 10, 10), seed=37)
        first = record.to_dict()
        self.assertEqual(first['formal_charge'], 2)
        self.assertTrue(all(a['isotope'] == 2 and a['formal_charge'] == 1 for a in first['atoms']))
        self.assertEqual(first, json.loads(json.dumps(first, allow_nan=False)))
        self.assertIn('python_implementation', first['software'])
        self.assertIn('polychemmd', first['software'])
        first['atoms'][0]['coordinates_angstrom'][0] = 1000
        self.assertNotEqual(first, record.to_dict())

    def test_extreme_finite_coordinates_fail_with_explicit_value_errors(self):
        # Each input is finite, but a proper rotation aligns the vector with x
        # and its magnitude is larger than the largest representable float.
        source = geometry(((1.1e308, 1.1e308, 1.1e308),
                           (-1.1e308, -1.1e308, -1.1e308)))
        matrix = ((1/math.sqrt(3),)*3,
                  (1/math.sqrt(2), -1/math.sqrt(2), 0),
                  (1/math.sqrt(6), 1/math.sqrt(6), -2/math.sqrt(6)))
        self.assertAlmostEqual(determinant(matrix), 1.0)
        with patch('polychemmd.builders.packing._rotation_matrix', return_value=matrix):
            with self.assertRaisesRegex(ValueError, 'rotated coordinates'):
                self.pack((source,), max_attempts_per_molecule=1)
        skewed = geometry(((1.7e308, 0, 0), (-1.7e308, 0, 0), (-1.7e308, 1, 0)))
        with self.assertRaisesRegex(ValueError, 'molecular extent'):
            self.pack((skewed,))
        with self.assertRaisesRegex(ValueError, 'wrapped coordinates'):
            _wrap((1e308, 0, 0), (1e-308, 10, 10))


if __name__ == '__main__':
    unittest.main()
