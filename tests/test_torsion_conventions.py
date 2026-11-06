import copy
import json
import math
from pathlib import Path
import unittest

from polychemmd.descriptors.conventions import (
    KCAL_TO_KJ, export_lammps_fourier, import_lammps_fourier,
    import_lammps_opls,
)
from polychemmd.descriptors.torsion import FourierTerm, TorsionPotential


def fourier_record():
    return {
        "style": "fourier", "units": "real", "m": 3,
        "terms": [
            {"k_kcal_mol": -0.4, "periodicity": 3, "phase_deg": -120.0},
            {"k_kcal_mol": 1.7, "periodicity": 1, "phase_deg": 47.0},
            {"k_kcal_mol": 0.0, "periodicity": 5, "phase_deg": 725.0},
        ],
    }


class LammpsFourierConventionTests(unittest.TestCase):
    def test_against_independently_written_lammps_equation(self):
        data = fourier_record()
        potential = import_lammps_fourier(data)
        for phi in (-7.8, -math.pi, -.13, .37, 1.19, 2.4, 7.3):
            expected_kcal = sum(term['k_kcal_mol'] *
                (1 + math.cos(term['periodicity']*phi - term['phase_deg']*math.pi/180))
                for term in data['terms'])
            self.assertAlmostEqual(potential.energy(phi), 4.184*expected_kcal, places=12)

    def test_real_unit_and_no_implicit_half_factor(self):
        data = {'style': 'fourier', 'units': 'real', 'm': 1,
                'terms': [{'k_kcal_mol': 1., 'periodicity': 1, 'phase_deg': 0.}]}
        potential = import_lammps_fourier(data)
        self.assertEqual(KCAL_TO_KJ, 4.184)
        self.assertAlmostEqual(potential.energy(0), 8.368, places=13)
        self.assertAlmostEqual(potential.energy(math.pi), 0., places=13)

    def test_degrees_are_not_interpreted_as_radians(self):
        data = {'style': 'fourier', 'units': 'real', 'm': 1,
                'terms': [{'k_kcal_mol': 2., 'periodicity': 1, 'phase_deg': 90.}]}
        potential = import_lammps_fourier(data)
        self.assertAlmostEqual(potential.energy(math.pi/2), 16.736, places=12)

    def test_all_terms_order_and_unwrapped_phases_preserved(self):
        data = fourier_record()
        potential = import_lammps_fourier(data)
        self.assertEqual(len(potential.terms), 3)
        self.assertEqual([t.periodicity for t in potential.terms], [3, 1, 5])
        self.assertEqual(potential.terms[2].k_kj_mol, 0.)
        self.assertAlmostEqual(potential.terms[2].phase_rad, 725*math.pi/180, places=14)
        self.assertEqual(data, fourier_record())

    def test_roundtrip_preserves_coefficients_and_profile(self):
        source = fourier_record()
        exported = export_lammps_fourier(import_lammps_fourier(source))
        self.assertEqual(exported['m'], source['m'])
        for before, after in zip(source['terms'], exported['terms']):
            self.assertEqual(before['periodicity'], after['periodicity'])
            self.assertAlmostEqual(before['k_kcal_mol'], after['k_kcal_mol'], places=14)
            self.assertAlmostEqual(before['phase_deg'], after['phase_deg'], places=12)
        restored = import_lammps_fourier(exported)
        for phi in (.1, 1.27, -3.8):
            self.assertAlmostEqual(import_lammps_fourier(source).energy(phi), restored.energy(phi), places=13)

    def test_canonical_roundtrip_preserves_negative_terms(self):
        source = TorsionPotential((FourierTerm(-1.3, 2, 12.4), FourierTerm(2.8, 5, -.9)))
        restored = import_lammps_fourier(export_lammps_fourier(source))
        for phi in (0., .75, 1.2, -2.7):
            self.assertAlmostEqual(source.energy(phi), restored.energy(phi), places=12)

    def test_zero_periodicity_is_explicitly_unsupported_even_if_zero_energy(self):
        for coefficient in (0., 2.):
            data = fourier_record()
            data['terms'][0].update(periodicity=0, k_kcal_mol=coefficient)
            with self.assertRaisesRegex(ValueError, 'constant-offset'):
                import_lammps_fourier(data)

    def test_count_types_and_count_mismatch(self):
        for count in (True, 3., 0, -2, 2):
            with self.subTest(count=count):
                data = fourier_record()
                data['m'] = count
                with self.assertRaises(ValueError):
                    import_lammps_fourier(data)

    def test_schema_style_and_unit_errors(self):
        for key, value in [('style', 'harmonic'), ('units', 'metal'), ('terms', '1 2 3')]:
            data = fourier_record()
            data[key] = value
            with self.assertRaises(ValueError):
                import_lammps_fourier(data)
        for key in ('style', 'm'):
            data = fourier_record()
            del data[key]
            with self.assertRaises(ValueError):
                import_lammps_fourier(data)
        data = fourier_record()
        data['source_units'] = 'real'
        with self.assertRaises(ValueError):
            import_lammps_fourier(data)

    def test_term_field_numeric_and_boolean_errors(self):
        for key, value in [('k_kcal_mol', True), ('k_kcal_mol', float('inf')),
                           ('k_kcal_mol', 10**400), ('phase_deg', float('nan')),
                           ('periodicity', True), ('periodicity', 1.), ('periodicity', -1)]:
            data = fourier_record()
            data['terms'][0][key] = value
            with self.subTest(key=key, value=str(value)), self.assertRaises(ValueError):
                import_lammps_fourier(data)
        data = fourier_record()
        data['terms'][0]['phase_rad'] = 0
        with self.assertRaisesRegex(ValueError, 'unknown'):
            import_lammps_fourier(data)

    def test_unrepresentable_conversions_fail_explicitly(self):
        data = fourier_record()
        data['terms'][0]['k_kcal_mol'] = 1e308
        with self.assertRaisesRegex(ValueError, 'floating-point'):
            import_lammps_fourier(data)
        source = TorsionPotential((FourierTerm(1., 1, 1e308),))
        with self.assertRaisesRegex(ValueError, 'floating-point'):
            export_lammps_fourier(source)
        source = TorsionPotential((FourierTerm(5e-324, 1, 0.),))
        with self.assertRaisesRegex(ValueError, 'floating-point'):
            export_lammps_fourier(source)

    def test_invalid_mapping_and_export_type(self):
        for value in (None, [], 'fourier', {1: 2}):
            with self.assertRaises(ValueError):
                import_lammps_fourier(value)
        with self.assertRaises(ValueError):
            export_lammps_fourier({'terms': []})

    def test_checked_in_example_matches_schema(self):
        path = Path(__file__).resolve().parents[1]/'examples'/'torsion_lammps_fourier.json'
        source = json.loads(path.read_text())
        potential = import_lammps_fourier(source)
        self.assertEqual(len(potential.terms), source['m'])


class LammpsOplsConventionTests(unittest.TestCase):
    def test_all_four_terms_and_half_factor_against_equation(self):
        coefficients = [1.2, -.7, 2.3, 0.]
        data = {'style': 'opls', 'units': 'real', 'k_kcal_mol': coefficients}
        potential = import_lammps_opls(data)
        self.assertEqual(len(potential.terms), 4)
        self.assertEqual([t.periodicity for t in potential.terms], [1, 2, 3, 4])
        for phi in (-3.2, -.73, 0., .4, 1.7, 2.8):
            a, b, c, d = coefficients
            expected = .5*(a*(1+math.cos(phi)) + b*(1-math.cos(2*phi)) +
                           c*(1+math.cos(3*phi)) + d*(1-math.cos(4*phi)))
            self.assertAlmostEqual(potential.energy(phi), 4.184*expected, places=12)

    def test_hand_calculated_quadrature_points(self):
        potential = import_lammps_opls({'style': 'opls', 'units': 'real', 'k_kcal_mol': [1., 2., 3., 4.]})
        self.assertAlmostEqual(potential.energy(0), 16.736, places=12)
        self.assertAlmostEqual(potential.energy(math.pi), 0., places=12)
        self.assertAlmostEqual(potential.energy(math.pi/2), 16.736, places=12)

    def test_zero_coefficients_are_not_removed(self):
        potential = import_lammps_opls({'style': 'opls', 'units': 'real', 'k_kcal_mol': [0., 0., 0., 0.]})
        self.assertEqual(len(potential.terms), 4)
        self.assertEqual(potential.energy(.8), 0.)

    def test_opls_schema_and_numeric_errors(self):
        good = {'style': 'opls', 'units': 'real', 'k_kcal_mol': [1., 2., 3., 4.]}
        for key, value in [('style', 'fourier'), ('units', 'metal'), ('k_kcal_mol', [1, 2, 3]),
                           ('k_kcal_mol', [True, 0., 0., 0.]),
                           ('k_kcal_mol', [float('nan'), 0., 0., 0.]),
                           ('k_kcal_mol', [1e308, 0., 0., 0.])]:
            data = copy.deepcopy(good)
            data[key] = value
            with self.assertRaises(ValueError):
                import_lammps_opls(data)
        with self.assertRaises(ValueError):
            import_lammps_opls(dict(good, phase_deg=0))
