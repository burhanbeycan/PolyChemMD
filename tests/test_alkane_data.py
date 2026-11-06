import math
import unittest

from polychemmd.descriptors.torsion import TorsionPotential
from polychemmd.forcefields.alkane_data import MODEL_ID, load_alkane_parameters
from polychemmd.forcefields.contracts import ParameterSet


class AlkaneSourceConversionTests(unittest.TestCase):
    def definitions(self, kind):
        data = load_alkane_parameters()
        return {entry['type_id']: entry['values'] for entry in data['parameter_set']['parameters'][kind]}

    def test_nm_and_half_prefactor_bond_conversion_preserves_reference_energy(self):
        values = self.definitions('bonds')
        # Independent calculation in the upstream nm/half-prefactor convention.
        for name, r0_nm, upstream_k in [('bond_CC', .15290, 224262.4),
                                       ('bond_CH', .10900, 284512.0)]:
            with self.subTest(name=name):
                record = values[name]
                r_angstrom = r0_nm*10 + .1
                expected = .5*upstream_k*(r_angstrom/10-r0_nm)**2
                actual = record['k_kj_mol_angstrom2']*(r_angstrom-record['r0_angstrom'])**2
                self.assertAlmostEqual(actual, expected, places=11)

    def test_exact_rounded_source_angle_constants_are_not_replaced(self):
        values = self.definitions('angles')
        for name, theta_deg, upstream_k in [('angle_CCC',112.7,488.273),
                                           ('angle_CCH',110.7,313.8),
                                           ('angle_HCH',107.8,276.144)]:
            with self.subTest(name=name):
                record = values[name]
                self.assertEqual(record['k_kj_mol_rad2'], upstream_k/2)
                self.assertAlmostEqual(record['theta0_rad'],theta_deg*math.pi/180,places=15)
                expected = .5*upstream_k*.1**2
                self.assertAlmostEqual(record['k_kj_mol_rad2']*.1**2,expected,places=13)

    def test_full_fourier_matches_independent_upstream_rb_polynomial(self):
        values = self.definitions('dihedrals')
        upstream = {'torsion_CCCC':(2.9288,-1.4644,.2092,-1.6736,0,0),
                    'torsion_CCCH':(.6276,1.8828,0,-2.5104,0,0),
                    'torsion_HCCH':(.6276,1.8828,0,-2.5104,0,0)}
        for name,coefficients in upstream.items():
            potential = TorsionPotential.from_mapping(values[name])
            for phi in (-2*math.pi,-2.4,-math.pi/3,0,.413,math.pi/2,math.pi,5.1):
                with self.subTest(name=name,phi=phi):
                    cosine=math.cos(phi-math.pi)
                    expected=math.fsum(c*cosine**n for n,c in enumerate(coefficients))
                    self.assertAlmostEqual(potential.energy(phi),expected,places=12)

    def test_fixed_charge_lj_and_opls_14_conventions_remain_explicit(self):
        data=load_alkane_parameters()
        parameters=ParameterSet.from_mapping(data['parameter_set'])
        self.assertEqual(data['charge_assignments'],{'opls_135':-.18,'opls_136':-.12,'opls_140':.06})
        self.assertEqual(parameters.nonbonded.mixing_rule,'geometric')
        self.assertEqual((parameters.nonbonded.lj_14_scale,parameters.nonbonded.coulomb_14_scale),(.5,.5))
        values=self.definitions('atoms')
        self.assertAlmostEqual(values['opls_135']['epsilon_kj_mol'],.066*4.184,places=14)
        self.assertAlmostEqual(values['opls_140']['epsilon_kj_mol'],.030*4.184,places=14)
        self.assertEqual((values['opls_135']['sigma_angstrom'],values['opls_140']['sigma_angstrom']),(3.5,2.5))

    def test_frozen_source_is_attributed_and_each_load_is_independent(self):
        data=load_alkane_parameters()
        self.assertEqual(data['model_id'],MODEL_ID)
        self.assertEqual(data['source_snapshot']['commit'],'820db9e499e4edc7f0dc734e6baeb0e9c3f87859')
        for group in data['parameter_set']['parameters'].values():
            for entry in group:
                self.assertEqual(entry['source']['kind'],'literature')
                self.assertIn('/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/',entry['source']['reference'])
        data['charge_assignments']['opls_135']=99
        self.assertEqual(load_alkane_parameters()['charge_assignments']['opls_135'],-.18)


if __name__ == '__main__':
    unittest.main()
