import math
import unittest

from polychemmd.analysis import analyze_density, density_from_mass_volume
from polychemmd.analysis.density import AVOGADRO_CONSTANT_MOL_INV


class DensityAnalysisTests(unittest.TestCase):
    def test_mass_volume_conversion_independent_water_example(self):
        # One water molecule in a volume chosen analytically for 1 g/cm^3.
        masses = [15.999, 1.008, 1.008]
        volume = sum(masses) / AVOGADRO_CONSTANT_MOL_INV / 1e-24
        self.assertAlmostEqual(density_from_mass_volume(masses, volume), 1.0, places=14)
        self.assertAlmostEqual(density_from_mass_volume([12.0], 1000.0), 0.01992646880608616, places=15)

    def test_mass_volume_requires_positive_finite_inputs(self):
        for masses, volume in [([], 1), ([0], 1), ([True], 1), ([1], 0), ([1], math.inf), ([math.nan], 1), ([10**400], 1)]:
            with self.subTest(masses=masses, volume=volume), self.assertRaises(ValueError):
                density_from_mass_volume(masses, volume)

    def test_exact_block_means_and_conditional_standard_error(self):
        series = [value for value in (.9, 1.1, .9, 1.1, .9, 1.1, .9, 1.1) for _ in range(2)]
        result = analyze_density(series, block_size=2)
        self.assertEqual(result['block_means_g_cm3'], [.9, 1.1]*4)
        self.assertAlmostEqual(result['mean_g_cm3'], 1.)
        expected_se = math.sqrt(.08 / 7) / math.sqrt(8)
        self.assertAlmostEqual(result['conditional_standard_error_g_cm3'], expected_se)
        self.assertFalse(result['conditional_normal_interval']['calibrated'])
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertEqual(result['polymer_density_validation'], 'not_performed')

    def test_discard_and_incomplete_tail_are_explicit(self):
        result = analyze_density([9., 9.]+[.99,1.01]*8+[8.], block_size=2, discard=2)
        self.assertEqual(result['n_used'], 16)
        self.assertEqual(result['n_omitted_trailing'], 1)
        self.assertEqual(result['n_discarded_prefix'], 2)
        self.assertEqual(result['mean_g_cm3'], 1.)
        self.assertEqual(result['sampling_status'], 'inconclusive')
        self.assertIsNone(result['conditional_normal_interval'])

    def test_constant_density_is_inconclusive_not_equilibrated(self):
        result = analyze_density([1.]*80, block_size=10)
        self.assertEqual(result['sampling_status'], 'inconclusive')
        self.assertIsNone(result['conditional_standard_error_g_cm3'])
        self.assertIsNone(result['lag_one_block_correlation'])
        self.assertTrue(result['flags'])

    def test_split_half_drift_requests_review(self):
        result = analyze_density([.8,.81,.8,.81,1.2,1.21,1.2,1.21], block_size=1)
        self.assertEqual(result['sampling_status'], 'requires_review')
        self.assertGreater(result['split_half']['relative_difference'], .02)

    def test_coarsened_blocks_do_not_claim_effective_independence(self):
        result = analyze_density([.9,.91,.92,.93,1.0,1.01,1.02,1.03], block_size=1)
        sensitivity = result['block_length_sensitivity']
        self.assertEqual(sensitivity['n_coarsened_blocks'], 4)
        self.assertIsNotNone(sensitivity['coarsened_standard_error_g_cm3'])
        self.assertIn('not established', sensitivity['interpretation'])

    def test_short_or_malformed_input_is_rejected(self):
        cases = [([], 1, {}), ([1.]*7, 1, {}), ([1.]*8, True, {}),
                 ([1.]*8, 1, {'discard': True}), ([1.]*8, 1, {'discard': 8}),
                 ([1.]*8, 1, {'min_blocks': 4}), ([1.]*8, 1, {'drift_relative_limit': 0}),
                 ([1.]*7+[float('nan')], 1, {}), ([1.]*7+[True], 1, {}),
                 ('11111111', 1, {})]
        for samples, block_size, kwargs in cases:
            with self.subTest(samples=samples, block_size=block_size, kwargs=kwargs), self.assertRaises(ValueError):
                analyze_density(samples, block_size, **kwargs)

    def test_output_is_json_finite_or_numeric_range_is_explicitly_rejected(self):
        import json
        for series in ([1e308, 1e307]*4, [1e-320, 2e-320]*4):
            with self.subTest(series=series):
                try:
                    result = analyze_density(series, block_size=1)
                except ValueError as exc:
                    self.assertIn('numeric range', str(exc))
                else:
                    json.dumps(result, allow_nan=False)


if __name__ == '__main__':
    unittest.main()
