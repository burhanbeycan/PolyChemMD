from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from polychemmd.analysis.structure import (
    _eigensystem_symmetric3, _gyration, _nematic, analyze_chain_structure,
)
from polychemmd.builders import tile_molecular_geometry
from polychemmd.forcefields import assign_periodic_alkanes, assign_compact_periodic_alkanes
from polychemmd.geometry import MolecularGeometry
from test_periodic_assignment import analytic_ethane, packed_ethanes, translate_molecule


ROOT = Path(__file__).resolve().parents[1]


def engine_state(assigned, *, rotation=None, translation=(0., 0., 0.)):
    record = assigned.to_dict(); lengths = record['box']['lengths_angstrom']
    rotation = rotation or [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    rows = []
    for index, atom in enumerate(record['atoms']):
        old = atom['coordinates_angstrom']
        point = [math.fsum(a*b for a, b in zip(row, old))+offset for row, offset in zip(rotation, translation)]
        image = [math.floor(v/l) for v, l in zip(point, lengths)]
        rows.append({'engine_id': index+1, 'cell_atom_id': atom['cell_atom_id'],
                     'molecule_index': atom['molecule_index'], 'charge_e': atom['charge_e'],
                     'coordinates_angstrom': point, 'wrapped_coordinates_angstrom': [v-i*l for v, i, l in zip(point, image, lengths)],
                     'image_flags': image, 'forces_kj_mol_angstrom': [0., 0., 0.]})
    return {'step': 200, 'atom_count': len(rows),
            'box': {key: record['box'][key] for key in ('origin_angstrom', 'lengths_angstrom', 'periodic')},
            'atoms': rows}


class StructureKernelTests(unittest.TestCase):
    def test_mass_weighted_two_point_radius_and_tensor(self):
        # Weights 1/4 and 3/4, center x=1.5, variance=3/4 Å².
        result = _gyration([[0, 0, 0], [2, 0, 0]], [1, 3])
        self.assertAlmostEqual(result['radius_of_gyration_angstrom'], math.sqrt(.75))
        self.assertEqual(result['gyration_tensor_angstrom2'], [[.75, 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        self.assertEqual(result['relative_shape_anisotropy'], 1.)

    def test_discrete_rod_matches_known_variance(self):
        # Five equally weighted sites 0,a,...,4a have variance 2*a².
        result = _gyration([[2*i, 0, 0] for i in range(5)], [1]*5)
        self.assertAlmostEqual(result['radius_of_gyration_angstrom'], math.sqrt(8))
        self.assertAlmostEqual(sum(result['gyration_eigenvalues_angstrom2']), 8.)

    def test_rotated_tensor_has_signed_cross_terms_without_sqrt(self):
        result = _gyration([[-1., 1., 0.], [1., -1., 0.]], [1, 1])
        self.assertEqual(result['gyration_tensor_angstrom2'][0][1], -1.)
        self.assertAlmostEqual(result['radius_of_gyration_angstrom'], math.sqrt(2))
        self.assertEqual(result['gyration_eigenvalues_angstrom2'], [2., 0., 0.])

    def test_isotropic_point_shape_has_zero_anisotropy(self):
        points = [[1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1]]
        result = _gyration(points, [1]*6)
        self.assertAlmostEqual(result['radius_of_gyration_angstrom'], 1.)
        self.assertAlmostEqual(result['relative_shape_anisotropy'], 0.)

    def test_parallel_and_antiparallel_axes_have_unit_order(self):
        result = _nematic([[1., 2., 3.], [-2., -4., -6.]])
        self.assertAlmostEqual(result['nematic_order_parameter'], 1.)
        self.assertEqual(result['orientation_status'], 'computed')
        self.assertAlmostEqual(abs(sum(a*b for a, b in zip(result['nematic_director'], [1/math.sqrt(14), 2/math.sqrt(14), 3/math.sqrt(14)]))), 1.)

    def test_three_orthogonal_axes_have_zero_order_and_no_director(self):
        result = _nematic([[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
        self.assertEqual(result['orientation_tensor'], [[0., 0., 0.], [0., 0., 0.], [0., 0., 0.]])
        self.assertEqual(result['nematic_order_parameter'], 0.)
        self.assertIsNone(result['nematic_director'])
        self.assertEqual(result['orientation_status'], 'inconclusive_degenerate_director')

    def test_planar_axis_degeneracy_does_not_invent_a_director(self):
        result = _nematic([[1., 0., 0.], [0., 1., 0.]])
        self.assertEqual(result['nematic_order_parameter'], .25)
        self.assertIsNone(result['nematic_director'])

    def test_subnormal_orientation_normalizes_without_false_superunit_order(self):
        result = _nematic([[5e-324, 5e-324, 5e-324]])
        self.assertAlmostEqual(result['nematic_order_parameter'], 1.)
        self.assertAlmostEqual(math.hypot(*result['nematic_director']), 1.)

    def test_zero_endpoint_is_inconclusive_without_dropping_chain(self):
        result = _nematic([[1., 0., 0.], [0., 0., 0.]])
        self.assertEqual(result['orientation_status'], 'inconclusive_degenerate_endpoints')
        self.assertIsNone(result['nematic_order_parameter'])
        self.assertIsNone(result['orientation_tensor'])

    def test_eigenvalues_and_vectors_match_nondiagonal_analytic_matrix(self):
        values, vectors, count = _eigensystem_symmetric3([[2., 1., 0.], [1., 2., 0.], [0., 0., 4.]])
        for actual, expected in zip(values, [4., 3., 1.]): self.assertAlmostEqual(actual, expected)
        self.assertLessEqual(count, 64)
        self.assertAlmostEqual(abs(vectors[1][0]), 1/math.sqrt(2))

    def test_numerical_kernels_reject_nonfinite_or_unrepresentable_statistics(self):
        for masses in ([True, 1], [1, -1], [1, math.inf]):
            with self.subTest(masses=masses), self.assertRaises(ValueError): _gyration([[0, 0, 0], [1, 0, 0]], masses)
        with self.assertRaises(ValueError): _gyration([[1e308, 0, 0], [-1e308, 0, 0]], [1, 1])
        with self.assertRaises(ValueError): _gyration([[0, 0, 0], [1e200, 0, 0]], [1, 1])
        with self.assertRaises(ValueError): _gyration([[0, 0, 0], [1e-200, 0, 0]], [1, 1])
        with self.assertRaises(ValueError): _nematic([[math.nan, 0, 0]])
        with self.assertRaises(ValueError): _eigensystem_symmetric3([[1, 2, 0], [0, 1, 0], [0, 0, 1]])


class ChainStructureTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.legacy = assign_periodic_alkanes(packed_ethanes())
        cls.compact = assign_compact_periodic_alkanes(packed_ethanes())

    def test_legacy_compact_same_snapshot_dimensions_and_ethane_backbone_oracle(self):
        first, second = analyze_chain_structure(self.legacy), analyze_chain_structure(self.compact)
        self.assertEqual(first['per_chain'], second['per_chain'])
        self.assertEqual(first['orientation_tensor'], second['orientation_tensor'])
        self.assertAlmostEqual(first['mean_end_to_end_distance_angstrom'], 1.53)
        self.assertAlmostEqual(first['mean_backbone_radius_of_gyration_angstrom'], .765)
        self.assertEqual(first['coordinate_source'], 'assigned_initial_geometry')
        self.assertEqual(first['equilibrium_status'], 'not_assessed')
        self.assertFalse(first['production_ready'])
        self.assertFalse(first['configuration_statistics_are_equilibrium_averages'])
        self.assertEqual(first['independent_material_cell_count'], 1)

    def test_rigid_rotation_translation_and_lattice_images_preserve_dimensions(self):
        base = analyze_chain_structure(self.compact, engine_state(self.compact))
        moved = analyze_chain_structure(self.compact, engine_state(self.compact,
            rotation=[[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]], translation=(123.2, -52.4, .7)))
        for name in ('mean_mass_weighted_radius_of_gyration_angstrom', 'mean_backbone_radius_of_gyration_angstrom',
                     'mean_end_to_end_distance_angstrom', 'nematic_order_parameter'):
            self.assertAlmostEqual(base[name], moved[name], places=12)

    def test_boundary_crossing_whole_molecule_translation_is_not_a_long_chain(self):
        source = packed_ethanes(); center = source['molecules'][0]['placement']['translated_centroid_angstrom']
        moved = translate_molecule(source, 0, [-center[0]+.01, 0, 0])
        self.assertTrue(any(atom['image_flags'] != [0, 0, 0] for atom in moved['atoms']))
        review = analyze_chain_structure(assign_periodic_alkanes(moved))
        self.assertAlmostEqual(review['mean_end_to_end_distance_angstrom'], 1.53)
        self.assertAlmostEqual(review['mean_backbone_radius_of_gyration_angstrom'], .765)

    def test_full_engine_state_hash_and_source_identity_are_explicit(self):
        state = engine_state(self.compact); state['atoms'].reverse()
        before = deepcopy(state); review = analyze_chain_structure(self.compact, state)
        self.assertEqual(state, before)
        expected = hashlib.sha256(json.dumps(state, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()
        self.assertEqual(review['state_canonical_sha256'], expected)
        self.assertEqual(review['coordinate_source'], 'supplied_engine_state')
        self.assertEqual(review['step'], 200)

    def test_missing_atom_duplicate_id_molecule_charge_mass_and_type_changes_fail(self):
        mutations = [lambda state: state['atoms'].pop(),
            lambda state: state['atoms'].__setitem__(1, deepcopy(state['atoms'][0])),
            lambda state: state['atoms'][0].__setitem__('engine_id', True),
            lambda state: state['atoms'][0].__setitem__('molecule_index', 99),
            lambda state: state['atoms'][0].__setitem__('charge_e', 0.),
            lambda state: state['atoms'][0].__setitem__('mass_g_mol', 999.),
            lambda state: state['atoms'][0].__setitem__('element', 'N')]
        for mutate in mutations:
            state = engine_state(self.compact); mutate(state)
            with self.subTest(mutation=mutations.index(mutate)), self.assertRaises(ValueError): analyze_chain_structure(self.compact, state)

    def test_inconsistent_or_individually_shifted_image_flag_fails(self):
        state = engine_state(self.compact); state['atoms'][0]['image_flags'][0] += 1
        with self.assertRaisesRegex(ValueError, 'image'): analyze_chain_structure(self.compact, state)
        state = engine_state(self.compact)
        state['atoms'][0]['image_flags'][0] += 1
        state['atoms'][0]['coordinates_angstrom'][0] += state['box']['lengths_angstrom'][0]
        with self.assertRaisesRegex(ValueError, 'bonded image'): analyze_chain_structure(self.compact, state)

    def test_nonfinite_optional_forces_and_invalid_box_or_summary_fail(self):
        state = engine_state(self.compact); state['atoms'][0]['forces_kj_mol_angstrom'][0] = math.inf
        with self.assertRaises(ValueError): analyze_chain_structure(self.compact, state)
        state = engine_state(self.compact); state['box']['periodic'] = [True, True, False]
        with self.assertRaises(ValueError): analyze_chain_structure(self.compact, state)
        with self.assertRaisesRegex(ValueError, 'full'): analyze_chain_structure(self.compact, {'atom_count': 16})
        with self.assertRaises(ValueError): analyze_chain_structure(self.compact.to_dict())

    def test_ordered_actual_decane_motif_reports_unit_order_without_phase_claim(self):
        source = json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        geometry = MolecularGeometry.from_mapping(source['molecules'][0]['source_geometry'])
        cell = tile_molecular_geometry(geometry, repeats=(2, 2, 2)).to_dict()
        result = analyze_chain_structure(assign_compact_periodic_alkanes(cell))
        self.assertAlmostEqual(result['nematic_order_parameter'], 1., places=13)
        self.assertEqual(result['molecule_count'], 8)
        self.assertEqual(result['phase_classification'], 'not_performed')
        self.assertFalse(result['molecules_are_independent_material_cells'])
        self.assertGreater(result['mean_mass_weighted_radius_of_gyration_angstrom'], 0.)

    def test_optional_packages_are_not_imported(self):
        code = ('import sys;sys.path.insert(0,sys.argv[1]);'
                'import polychemmd.analysis.structure;'
                'assert not any(n in sys.modules for n in ("rdkit","numpy","lammps"))')
        result = subprocess.run([sys.executable, '-S', '-c', code, str(ROOT/'src')],
                                cwd=tempfile.gettempdir(), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__': unittest.main()
