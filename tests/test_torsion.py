"""Analytic scientific checks for the independent Fourier descriptor core."""

import json
import math
import unittest

from polychemmd.descriptors import FourierTerm, TorsionPotential
from polychemmd.descriptors.torsion import CONVENTION


class TorsionTests(unittest.TestCase):
    def test_periodic_energy_and_phase(self):
        potential = TorsionPotential((FourierTerm(2.1, 3, 0.7), FourierTerm(0.4, 2, -0.2)))
        for angle in (-8.2, -0.3, 0.0, 1.1, 8.4):
            self.assertAlmostEqual(potential.energy(angle), potential.energy(angle + math.tau), places=12)

    def test_known_threefold_potential(self):
        result = TorsionPotential((FourierTerm(2.0, 3),)).describe()
        self.assertAlmostEqual(result["energy_min_grid_kj_mol"], 0.0, places=12)
        self.assertAlmostEqual(result["energy_max_grid_kj_mol"], 4.0, places=12)
        self.assertAlmostEqual(result["barrier_grid_kj_mol"], 4.0, places=12)
        self.assertEqual(result["n_minima_grid"], 3)

    def test_all_fourier_terms_are_preserved(self):
        potential = TorsionPotential((FourierTerm(2.0, 1), FourierTerm(1.0, 2)))
        result = potential.describe()
        self.assertEqual(len(result["terms"]), 2)
        self.assertAlmostEqual(potential.energy(math.pi), 2.0)
        self.assertAlmostEqual(result["barrier_grid_kj_mol"], 4.5)
        self.assertEqual(result["n_minima_grid"], 2)
        self.assertNotEqual(potential.energy(math.pi), TorsionPotential((FourierTerm(2.0, 1),)).energy(math.pi))

    def test_term_order_does_not_change_energy(self):
        terms = (FourierTerm(1.2, 2, 0.2), FourierTerm(0.6, 1, 1.2), FourierTerm(-0.1, 4))
        self.assertAlmostEqual(TorsionPotential(terms).energy(0.41), TorsionPotential(tuple(reversed(terms))).energy(0.41))

    def test_signed_coefficients_are_explicitly_supported(self):
        result = TorsionPotential((FourierTerm(-2.0, 1),)).describe()
        self.assertAlmostEqual(result["energy_min_grid_kj_mol"], -4.0)
        self.assertAlmostEqual(result["energy_max_grid_kj_mol"], 0.0)
        self.assertAlmostEqual(result["barrier_grid_kj_mol"], 4.0)
        self.assertEqual(result["n_minima_grid"], 1)
        self.assertEqual(result["minima_angles_grid_rad"], [0.0])

    def test_flat_potential_has_no_isolated_minima(self):
        result = TorsionPotential((FourierTerm(0.0, 3),)).describe(n_grid=36)
        self.assertEqual(result["n_minima_grid"], 0)
        self.assertEqual(result["barrier_grid_kj_mol"], 0.0)
        self.assertAlmostEqual(result["shannon_grid_entropy_nats"], math.log(36))
        self.assertAlmostEqual(result["shannon_grid_entropy_fraction_of_uniform"], 1.0)
        self.assertTrue(all(abs(p - 1 / 36) < 1e-14 for p in result["boltzmann_bin_probabilities"]))

    def test_two_bin_minimum_plateau_counts_once(self):
        result = TorsionPotential((FourierTerm(1.0, 1, math.pi / 36),)).describe(n_grid=36)
        self.assertEqual(result["n_minima_grid"], 1)
        self.assertGreater(result["energy_min_grid_kj_mol"], 0.0)

    def test_cyclic_boundary_plateau_counts_once(self):
        result = TorsionPotential((FourierTerm(-1.0, 1, -math.pi / 36),)).describe(n_grid=36)
        self.assertEqual(result["n_minima_grid"], 1)

    def test_boltzmann_normalization_and_temperature(self):
        potential = TorsionPotential((FourierTerm(3.0, 1),))
        cold = potential.describe(temperature_k=100)
        hot = potential.describe(temperature_k=1000)
        self.assertAlmostEqual(math.fsum(cold["boltzmann_bin_probabilities"]), 1.0)
        self.assertAlmostEqual(math.fsum(hot["boltzmann_bin_probabilities"]), 1.0)
        self.assertGreater(hot["shannon_grid_entropy_nats"], cold["shannon_grid_entropy_nats"])
        self.assertGreater(hot["mean_relative_energy_grid_kj_mol"], cold["mean_relative_energy_grid_kj_mol"])

    def test_stable_weights_for_large_energy_separation(self):
        result = TorsionPotential((FourierTerm(1e6, 2),)).describe(temperature_k=1.0)
        self.assertAlmostEqual(math.fsum(result["boltzmann_bin_probabilities"]), 1.0)
        self.assertTrue(all(math.isfinite(p) and p >= 0 for p in result["boltzmann_bin_probabilities"]))
        json.dumps(result, allow_nan=False)

    def test_mapping_and_explicit_convention(self):
        potential = TorsionPotential.from_mapping({"convention": CONVENTION, "terms": [{"k_kj_mol": 2, "periodicity": 3}]})
        self.assertEqual(potential.terms, (FourierTerm(2.0, 3, 0.0),))
        self.assertEqual(potential.describe()["convention"], CONVENTION)

    def test_invalid_term_values(self):
        for coefficient in (True, "2", math.inf, math.nan):
            with self.subTest(coefficient=coefficient), self.assertRaises(ValueError):
                FourierTerm(coefficient, 1)
        for periodicity in (True, False, 1.0, 0, -1, "3"):
            with self.subTest(periodicity=periodicity), self.assertRaises(ValueError):
                FourierTerm(1, periodicity)
        for phase in (True, math.inf, math.nan, "pi"):
            with self.subTest(phase=phase), self.assertRaises(ValueError):
                FourierTerm(1, 1, phase)

    def test_invalid_potential_mappings(self):
        cases = [
            {}, {"terms": []}, {"terms": "bad"}, {"terms": ["bad"]},
            {"terms": [{"k_kj_mol": 1, "periodicity": 1, "phase_deg": 180}]},
            {"terms": [{"k_kj_mol": 1, "periodicity": 1}], "convention": "half_prefactor"},
            {"terms": [{"k_kj_mol": 1, "periodicity": 1}], "temperature": 300},
        ]
        for data in cases:
            with self.subTest(data=data), self.assertRaises(ValueError):
                TorsionPotential.from_mapping(data)

    def test_invalid_grid_temperature_or_angle(self):
        potential = TorsionPotential((FourierTerm(1, 3),))
        for grid in (True, 12, 36.0, -1):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                potential.describe(n_grid=grid)
        with self.assertRaisesRegex(ValueError, "four bins"):
            TorsionPotential((FourierTerm(1, 20),)).describe(n_grid=36)
        for temperature in (True, 0, -1, math.nan, math.inf, 5e-324):
            with self.subTest(temperature=temperature), self.assertRaises(ValueError):
                potential.describe(temperature_k=temperature)
        for angle in (True, math.inf, math.nan):
            with self.subTest(angle=angle), self.assertRaises(ValueError):
                potential.energy(angle)

    def test_nonfinite_evaluation_is_rejected(self):
        potential = TorsionPotential((FourierTerm(1e308, 1),))
        with self.assertRaisesRegex(ValueError, "floating-point"):
            potential.energy(0.0)

    def test_terms_are_immutable(self):
        original = [FourierTerm(1, 1)]
        potential = TorsionPotential(original)
        original.append(FourierTerm(2, 2))
        self.assertEqual(len(potential.terms), 1)


if __name__ == "__main__":
    unittest.main()
