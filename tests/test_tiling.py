"""Independent geometric fixtures; ordered tiling is not a phase benchmark."""

from itertools import combinations, product
import json
import math
from pathlib import Path
import random
import unittest

from polychemmd.builders.packing import PackingError
from polychemmd.builders.tiling import _Budget, _contact_audit, tile_molecular_geometry
from polychemmd.geometry import MolecularGeometry
from test_packing import geometry, signed_volume, tetrahedron


ROOT = Path(__file__).resolve().parents[1]


def source_decane():
    record = json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
    return MolecularGeometry.from_mapping(record['molecules'][0]['source_geometry'])


class TilingConstructionTests(unittest.TestCase):
    def tile(self, source=None, **options):
        defaults = {'repeats': (2, 2, 2), 'spacing_angstrom': (4, 5, 6),
                    'min_distance_angstrom': 1.0}
        defaults.update(options)
        return tile_molecular_geometry(tetrahedron() if source is None else source, **defaults).to_dict()

    def test_all_small_cell_pairs_independently_pass_27_image_audit(self):
        result = self.tile()
        box = result['box']['lengths_angstrom']
        for a, b in combinations(result['atoms'], 2):
            if a['molecule_index'] == b['molecule_index']:
                continue
            nearest = min(math.dist(a['wrapped_coordinates_angstrom'],
                                    tuple(v+i*l for v, i, l in zip(b['wrapped_coordinates_angstrom'], image, box)))
                          for image in product((-1, 0, 1), repeat=3))
            self.assertGreaterEqual(nearest, 1.0)
        self.assertTrue(result['validation']['intermolecular_minimum_image_audit'])
        self.assertIsNone(result['packing']['minimum_inter_molecular_distance_angstrom'])
        self.assertIn('not_measured', result['packing']['global_minimum_distance_status'])

    def test_proper_rotation_preserves_all_pair_distances_and_handedness(self):
        source = tetrahedron()
        matrix = ((0, -1, 0), (1, 0, 0), (0, 0, 1))
        result = self.tile(source, rotation_matrix=matrix)
        for index in range(result['molecule_count']):
            points = [row['coordinates_angstrom'] for row in result['atoms'] if row['molecule_index'] == index]
            for i, j in combinations(range(len(points)), 2):
                self.assertAlmostEqual(math.dist(points[i], points[j]),
                                       math.dist(source.coordinates_angstrom[i], source.coordinates_angstrom[j]), places=12)
            self.assertAlmostEqual(signed_volume(points), signed_volume(source.coordinates_angstrom), places=12)

    def test_global_and_source_identity_maps_bonds_and_images_survive(self):
        source = tetrahedron()
        result = self.tile(source, spacing_angstrom=(3, 4, 5))
        self.assertEqual([row['atom_map'] for row in result['atoms']], list(range(1, 33)))
        self.assertEqual(len({row['cell_atom_id'] for row in result['atoms']}), 32)
        box = result['box']['lengths_angstrom']
        for row in result['atoms']:
            atom = source.atoms[row['source_atom_map']-1]
            self.assertEqual(row['source_atom_id'], atom.atom_id)
            self.assertEqual(row['lineage'], atom.lineage)
            for v, w, image, length in zip(row['coordinates_angstrom'], row['wrapped_coordinates_angstrom'], row['image_flags'], box):
                self.assertGreaterEqual(w, 0); self.assertLess(w, length)
                self.assertAlmostEqual(v, w+image*length)
        self.assertEqual(len(result['bonds']), 8*len(source.bonds))
        self.assertTrue(all(m['source_geometry'] == source.to_dict() for m in result['molecules']))

    def test_density_and_scope_are_construction_inputs_only(self):
        result = self.tile(geometry(), repeats=(2, 1, 1), spacing_angstrom=(4, 5, 6))
        expected = 2*4.0026/(6.02214076e23*240e-24)
        self.assertAlmostEqual(result['construction_density_g_cm3'], expected, places=15)
        self.assertEqual(result['phase'], 'not_assigned')
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertFalse(result['simulation_performed'])
        self.assertEqual(result['independent_material_cell_count'], 1)
        self.assertFalse(result['repeated_molecules_are_independent_material_cells'])
        self.assertIsNone(result['packing']['seed'])

    def test_repeatability_global_rng_isolation_and_defensive_serialization(self):
        state = random.getstate()
        first, second = self.tile(), self.tile()
        self.assertEqual(first, second)
        self.assertEqual(random.getstate(), state)
        record = tile_molecular_geometry(geometry(), repeats=(2, 1, 1), spacing_angstrom=(3, 3, 3))
        modified = record.to_dict(); modified['atoms'][0]['coordinates_angstrom'][0] = 123
        self.assertNotEqual(modified, record.to_dict())

    def test_dense_decane_fixture_is_admissible_large_ordered_starting_cell(self):
        result = tile_molecular_geometry(source_decane()).to_dict()
        self.assertEqual(result['atom_count'], 2304)
        self.assertEqual(result['molecule_count'], 72)
        self.assertAlmostEqual(result['construction_density_g_cm3'], .6987592396767457, places=13)
        self.assertGreater(min(result['box']['lengths_angstrom']), 2*12.)
        self.assertGreater(result['packing']['minimum_checked_candidate_distance_angstrom'], 1.5)
        self.assertLess(result['packing']['contact_distance_checks_used'], result['atom_count']*20)
        self.assertIn('ordered repeated', result['construction_interpretation'])

    def test_small_tiling_passes_independent_frozen_parameter_source_validation(self):
        from polychemmd.forcefields.periodic import assign_periodic_alkanes
        cell = tile_molecular_geometry(source_decane(), repeats=(2, 2, 2)).to_dict()
        assigned = assign_periodic_alkanes(cell).to_dict()
        self.assertEqual(len(assigned['atoms']), 256)
        self.assertEqual(assigned['source_cell'], cell)
        self.assertTrue(assigned['validation']['source_graph_bonds_and_lineage_preserved'])
        self.assertAlmostEqual(math.fsum(row['charge_e'] for row in assigned['atoms']), 0., places=12)

    def test_earlier_denser_y_spacing_fails_fixed_geometric_floor(self):
        with self.assertRaisesRegex(PackingError, 'contact floor'):
            tile_molecular_geometry(source_decane(), spacing_angstrom=(14, 4.9, 4.9))

    def test_single_molecule_self_image_bound_is_not_skipped(self):
        with self.assertRaisesRegex(PackingError, 'self-image bound'):
            self.tile(geometry(((0, 0, 0), (5, 0, 0))), repeats=(1, 1, 1),
                      spacing_angstrom=(5, 6, 6), min_distance_angstrom=1)

    def test_translation_resolution_loss_is_rejected(self):
        with self.assertRaisesRegex(PackingError, 'pair distances'):
            self.tile(geometry(((0, 0, 0), (1, 0, 0))), repeats=(1, 1, 1),
                      spacing_angstrom=(1e100, 5, 5))

    def test_budget_includes_rigidity_and_contact_work(self):
        with self.assertRaisesRegex(PackingError, 'max_distance_checks'):
            self.tile(max_distance_checks=5)
        result = self.tile(geometry(), repeats=(2, 1, 1), spacing_angstrom=(2, 2, 2),
                           min_distance_angstrom=1.5, max_distance_checks=1)
        self.assertEqual(result['packing']['distance_checks_used'], 1)
        with self.assertRaisesRegex(PackingError, 'max_distance_checks'):
            self.tile(geometry(), repeats=(3, 1, 1), spacing_angstrom=(2, 2, 2),
                      min_distance_angstrom=1.5, max_distance_checks=1)

    def test_boolean_shape_nonfinite_and_derived_range_errors(self):
        for options in ({'repeats': (True, 2, 2)}, {'repeats': (0, 1, 1)}, {'repeats': (2, 2)},
                        {'repeats': (10**100, 1, 1)}, {'spacing_angstrom': (1e308,)*3},
                        {'spacing_angstrom': (1e-300,)*3}, {'spacing_angstrom': (True, 4, 4)},
                        {'spacing_angstrom': (float('nan'), 4, 4)}, {'max_atoms': False},
                        {'max_atoms': 8193}, {'max_distance_checks': False}, {'min_distance_angstrom': True}):
            with self.subTest(options=options), self.assertRaises(ValueError):
                self.tile(geometry(), **options)
        with self.assertRaises(ValueError):
            tile_molecular_geometry(object())

    def test_nonproper_or_malformed_rotations_fail(self):
        for matrix in (((-1, 0, 0), (0, 1, 0), (0, 0, 1)), ((2, 0, 0), (0, 1, 0), (0, 0, 1)),
                       ((True, 0, 0), (0, 1, 0), (0, 0, 1)), ((float('inf'), 0, 0), (0, 1, 0), (0, 0, 1)),
                       ((1, 0, 0), (0, 1, 0)), 'identity'):
            with self.subTest(matrix=matrix), self.assertRaises(ValueError):
                self.tile(rotation_matrix=matrix)


class LinkedBinTests(unittest.TestCase):
    def test_one_and_two_bin_periodicity_deduplicates_actual_checks(self):
        for box, points in (((.8, .8, .8), [(.1, .1, .1), (.5, .5, .5)]),
                            ((2, .8, .8), [(.25, .25, .25), (1.5, .25, .25)])):
            budget = _Budget(1)
            bins, checks, minimum = _contact_audit([[points[0]], [points[1]]], box, .5, budget)
            self.assertEqual(checks, 1)
            self.assertEqual(budget.used, 1)
            self.assertGreaterEqual(minimum, .5)
            self.assertIn(1, bins)

    def test_seeded_synthetic_contact_classification_matches_independent_image_scan(self):
        rng = random.Random(952710)
        for index in range(64):
            box = tuple(rng.uniform(1, 8) for _ in range(3))
            cutoff = rng.uniform(.1, min(box)*.8)
            groups = [[tuple(rng.random()*length for length in box) for _ in range(2)] for _ in range(4)]
            expected_collision = any(
                min(math.dist(a, tuple(v+i*l for v, i, l in zip(b, image, box)))
                    for image in product((-1, 0, 1), repeat=3)) < cutoff
                for first, second in combinations(groups, 2) for a in first for b in second)
            budget = _Budget(100)
            try:
                _contact_audit(groups, box, cutoff, budget)
                collision = False
            except PackingError:
                collision = True
            self.assertEqual(collision, expected_collision, index)
            self.assertLessEqual(budget.used, 24)

    def test_contact_across_periodic_face_rejected(self):
        with self.assertRaisesRegex(PackingError, 'contact floor'):
            _contact_audit([[(.1, 1, 1)], [(9.9, 1, 1)]], (10, 10, 10), .5, _Budget(10))

    def test_intermolecular_floor_does_not_exclude_source_bonded_points(self):
        _, checks, _ = _contact_audit([[(1, 1, 1), (1.1, 1, 1)]], (10, 10, 10), 1., _Budget(1))
        self.assertEqual(checks, 0)

    def test_tiny_finite_contact_floor_never_overflows_bin_ratio(self):
        result = tile_molecular_geometry(geometry(), repeats=(2, 1, 1), spacing_angstrom=(2, 2, 2),
                                         min_distance_angstrom=5e-324).to_dict()
        self.assertEqual(result['packing']['linked_bin_counts'], [2, 2, 2])
        json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
