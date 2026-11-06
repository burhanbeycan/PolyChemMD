import copy
import json
import math
from pathlib import Path
import unittest

from polychemmd.forcefields import (
    MissingParametersError, ParameterSet, TopologyRequirements, check_coverage,
)
from polychemmd.forcefields.contracts import ParameterDefinition, ParameterSource, NonbondedConvention

EXAMPLE = Path(__file__).resolve().parents[1] / 'examples' / 'parameter_contract_demo.json'


class ParameterContractTests(unittest.TestCase):
    def setUp(self):
        data = json.loads(EXAMPLE.read_text())
        self.parameters = data['parameter_set']
        self.topology = data['topology']

    def test_roundtrip_preserves_conventions_and_provenance(self):
        parsed = ParameterSet.from_mapping(self.parameters)
        self.assertEqual(parsed.to_dict(), self.parameters)
        self.assertEqual(TopologyRequirements.from_mapping(self.topology).to_dict(), self.topology)

    def test_complete_is_coverage_only_not_material_validation(self):
        report = check_coverage(self.parameters, self.topology)
        self.assertTrue(report.complete)
        self.assertIs(report.require_complete(), report)
        data = report.to_dict()
        self.assertEqual(data['topology_counts']['atoms'], 4)
        self.assertFalse(data['production_ready'])
        self.assertEqual(data['chemical_typing'], 'not_performed')
        self.assertEqual(data['physical_accuracy'], 'not_assessed')

    def test_every_interaction_kind_is_checked(self):
        for kind in ('atoms', 'bonds', 'angles', 'dihedrals', 'impropers'):
            with self.subTest(kind=kind):
                data = copy.deepcopy(self.parameters)
                data['parameters'][kind] = []
                report = check_coverage(data, self.topology)
                self.assertFalse(report.complete)
                self.assertTrue(report.to_dict()['missing_types'][kind])
                with self.assertRaises(MissingParametersError) as caught:
                    report.require_complete()
                self.assertIs(caught.exception.report, report)

    def test_missing_is_different_from_unused(self):
        self.topology['atoms'][0]['type_id'] = 'MISSING'
        extra = copy.deepcopy(self.parameters['parameters']['atoms'][0])
        extra['type_id'] = 'UNUSED'
        self.parameters['parameters']['atoms'].append(extra)
        report = check_coverage(self.parameters, self.topology).to_dict()
        self.assertEqual(report['missing_types']['atoms'], ['MISSING'])
        self.assertEqual(report['unused_types']['atoms'], ['UNUSED'])
        self.assertEqual(report['missing_interactions'], [{'kind':'atoms', 'id':'1', 'type_id':'MISSING'}])

    def test_all_registries_are_explicit_even_when_empty(self):
        del self.parameters['parameters']['impropers']
        with self.assertRaisesRegex(ValueError, 'missing'):
            ParameterSet.from_mapping(self.parameters)
        del self.topology['angles']
        with self.assertRaisesRegex(ValueError, 'missing'):
            TopologyRequirements.from_mapping(self.topology)

    def test_no_default_one_four_scaling_or_double_counting(self):
        cases = [
            lambda x: x['nonbonded'].pop('one_four'),
            lambda x: x['nonbonded']['one_four'].pop('coulomb_scale'),
            lambda x: x['nonbonded']['one_four'].update(handling='explicit_plus_scaled'),
            lambda x: x['nonbonded'].update(excluded_bonded_distances=[1,2,3]),
            lambda x: x['nonbonded'].update(explicit_pairs_14=[]),
        ]
        for mutate in cases:
            with self.subTest(mutate=mutate):
                data = copy.deepcopy(self.parameters)
                mutate(data)
                with self.assertRaises(ValueError):
                    ParameterSet.from_mapping(data)

    def test_scaling_is_finite_in_unit_interval_and_not_boolean(self):
        for value in (True, -0.1, 1.1, math.nan, math.inf, 10**400):
            with self.subTest(value=value):
                self.parameters['nonbonded']['one_four']['lj_scale'] = value
                with self.assertRaises(ValueError):
                    ParameterSet.from_mapping(self.parameters)

    def test_mixed_units_and_wrong_functional_forms_are_rejected(self):
        self.parameters['units']['energy'] = 'kcal/mol'
        with self.assertRaisesRegex(ValueError, 'units'):
            ParameterSet.from_mapping(self.parameters)
        self.parameters['units']['energy'] = 'kJ/mol'
        self.parameters['parameters']['bonds'][0]['form'] = 'half_k_delta_r_squared'
        with self.assertRaisesRegex(ValueError, 'form'):
            ParameterSet.from_mapping(self.parameters)

    def test_wrong_or_additional_numeric_field_is_rejected(self):
        values = self.parameters['parameters']['angles'][0]['values']
        values['theta0_deg'] = 90
        with self.assertRaisesRegex(ValueError, 'unknown'):
            ParameterSet.from_mapping(self.parameters)

    def test_invalid_parameter_magnitudes(self):
        for kind, field, value in [('atoms','mass_g_mol',0), ('atoms','sigma_angstrom',-1),
                                   ('atoms','epsilon_kj_mol',-1), ('bonds','r0_angstrom',True),
                                   ('angles','theta0_rad',4), ('impropers','chi0_rad',4)]:
            data = copy.deepcopy(self.parameters)
            data['parameters'][kind][0]['values'][field] = value
            with self.subTest(kind=kind, field=field):
                with self.assertRaises(ValueError):
                    ParameterSet.from_mapping(data)

    def test_duplicate_parameter_types_and_topology_ids_are_rejected(self):
        self.parameters['parameters']['atoms'].append(copy.deepcopy(self.parameters['parameters']['atoms'][0]))
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            ParameterSet.from_mapping(self.parameters)
        self.topology['atoms'].append({'id':1, 'type_id':'X'})
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            TopologyRequirements.from_mapping(self.topology)

    def test_empty_topology_is_invalid(self):
        for key in self.topology:
            self.topology[key] = []
        with self.assertRaisesRegex(ValueError, 'contain atoms'):
            TopologyRequirements.from_mapping(self.topology)

    def test_provenance_is_required_at_set_and_entry_level(self):
        del self.parameters['parameters']['bonds'][0]['source']
        with self.assertRaisesRegex(ValueError, 'missing'):
            ParameterSet.from_mapping(self.parameters)

    def test_values_are_defensively_copied(self):
        parameters = ParameterSet.from_mapping(self.parameters)
        atom = parameters.definitions[0]
        values = atom.values
        values['sigma_angstrom'] = 100
        self.assertEqual(atom.values['sigma_angstrom'], 1.0)

    def test_full_fourier_profile_is_retained(self):
        terms = self.parameters['parameters']['dihedrals'][0]['values']['terms']
        terms.append({'k_kj_mol':-2.0,'periodicity':1,'phase_rad':.4})
        parameters = ParameterSet.from_mapping(self.parameters)
        entry = next(x for x in parameters.definitions if x.kind == 'dihedrals')
        self.assertEqual(entry.values['terms'], terms)

    def test_direct_constructors_cannot_bypass_validation(self):
        with self.assertRaises(ValueError):
            TopologyRequirements(())
        with self.assertRaises(ValueError):
            NonbondedConvention('geometric', True, 1.)
        with self.assertRaises(ValueError):
            ParameterSource('qm_fit', '')
        source = ParameterSource('analytic_demo', 'direct constructor test')
        with self.assertRaises(ValueError):
            ParameterDefinition('atoms', 'X', 'lj12_6', '{"epsilon_kj_mol":-1}', source)
