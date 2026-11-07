"""Actual-record rejection checks and explicitly synthetic review fixtures."""
from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import random
import unittest

from polychemmd.analysis.readiness import (
    PreparationReadinessSpec, analyze_preparation_readiness,
)
from polychemmd.analysis.sampling import DensitySamplingSpec
from polychemmd.forcefields.periodic import assign_periodic_alkanes
from polychemmd.protocols.alkane import BAR_PER_ATM, _THERMO_COLUMNS
from polychemmd.protocols.numerics import NumericalQualificationSpec, compare_numerical_records
from test_numerics import synthetic_records


ROOT = Path(__file__).resolve().parents[1]


def digest(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(',', ':'), allow_nan=False).encode()).hexdigest()


def numerical_fixture(assigned, prep):
    """Synthetic numerical records, never an executed physical calculation."""
    assignment = assigned.to_dict()
    spec = NumericalQualificationSpec(cutoffs_angstrom=(6, 8, 9, 9.5, 10),
                                     reciprocal_nmax=12, refined_reciprocal_nmax=16)
    ids = [atom['cell_atom_id'] for atom in assignment['atoms']]
    rows = synthetic_records(spec=spec, count=len(ids), ids=ids)
    geometry = {'box': assignment['box'], 'atoms': [{key: atom[key] for key in (
        'cell_atom_id', 'molecule_index', 'charge_e', 'coordinates_angstrom',
        'wrapped_coordinates_angstrom', 'image_flags')} for atom in assignment['atoms']]}
    for row in rows:
        row.update(assignment_canonical_sha256=digest(assignment), geometry_canonical_sha256=digest(geometry),
                   model_id=assignment['model_id'])
        row['engine_run']['engine_identity'] = deepcopy(prep['engine_run']['engine_identity'])
    review = compare_numerical_records(rows, spec)
    return {'schema_version': '1', 'kind': 'periodic_numerical_qualification',
        'status': 'completed_numerical_review', 'spec': spec.to_dict(),
        'spec_canonical_sha256': digest(spec.to_dict()), 'review': review,
        'stages': [{'status': row['status'], 'record': row} for row in rows],
        **{key: deepcopy(review[key]) for key in ('assignment_canonical_sha256',
            'geometry_canonical_sha256', 'model_id', 'atom_count', 'engine_identity',
            'numerical_status', 'selected_settings')}}


def synthetic_long_record(assigned, original):
    """Schema-consistent synthetic NPT series; not a molecular trajectory."""
    report = deepcopy(original)
    report['protocol_spec'].update(nvt_steps=1, npt_steps=2047)
    report['md_integration_steps_executed'] = 2048
    report['final_state']['step'] = 2048
    final = report['final_state']
    final.update(temperature_k=300., pressure_bar=1.)
    masses = [atom['mass_g_mol'] for atom in assigned.to_dict()['atoms']]
    mass = math.fsum(masses)
    base = final['density_g_cm3']
    rng = random.Random(901)
    rows = []
    for i in range(2048):
        density = base+rng.uniform(-1e-5, 1e-5)
        rows.append({'step': i+1, 'temperature_k': 300+rng.uniform(-.5, .5),
            'pressure_bar': 1+rng.uniform(-.5, .5), 'density_g_cm3': density,
            'volume_angstrom3': mass/(.602214076*density)})
    rows[-1].update(temperature_k=300., pressure_bar=1., density_g_cm3=base,
                    volume_angstrom3=final['volume_angstrom3'])
    nvt = [dict(rows[0], step=0), deepcopy(rows[0])]
    report['stages']['nvt'] = {'status': 'executed', 'steps': 1,
                              'elapsed_time_fs': .25, 'thermo_rows': nvt}
    report['stages']['npt'] = {'status': 'executed', 'steps': 2047,
                              'elapsed_time_fs': 2047*.25, 'thermo_rows': rows}
    def raw(row):
        volume = row['volume_angstrom3']; length = volume**(1/3)
        return [row['step'], row['temperature_k'], row['pressure_bar']/BAR_PER_ATM,
                volume, mass/(.602214129*volume), 1., 2., .1, .2, .3, .4, .5, -.5,
                length, length, length]
    report['engine_run']['thermo_blocks'][1] = {'columns': list(_THERMO_COLUMNS), 'rows': [raw(row) for row in nvt]}
    report['engine_run']['thermo_blocks'][2] = {'columns': list(_THERMO_COLUMNS), 'rows': [raw(row) for row in rows]}
    return report


class ReadinessSpecTests(unittest.TestCase):
    def test_roundtrip_explicit_force_policy_and_no_thermal_zero_force_default(self):
        spec = PreparationReadinessSpec()
        self.assertEqual(PreparationReadinessSpec.from_mapping(spec.to_dict()), spec)
        self.assertIsNone(spec.max_force_norm_kj_mol_angstrom)
        self.assertIsNone(spec.max_rms_force_kj_mol_angstrom)
        self.assertTrue(spec.require_minimizer_force_convergence)

    def test_invalid_numbers_targets_unknown_fields_and_nested_mismatch(self):
        for values in ({'phase_target': 'verified_liquid'}, {'unknown': 1},
            {'temperature_k': True}, {'pressure_bar': float('nan')}, {'timestep_fs': 1.1},
            {'require_minimizer_force_convergence': 1}, {'max_nematic_order': 1.1},
            {'max_force_norm_kj_mol_angstrom': -1}, {'max_rms_force_kj_mol_angstrom': True},
            {'sampling_spec': {'target_temperature_k': 301.}}):
            with self.subTest(values=values), self.assertRaises(ValueError):
                PreparationReadinessSpec.from_mapping(values)


class PreparationReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cell = json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        cls.assigned = assign_periodic_alkanes(cell)
        cls.prep = json.loads((ROOT/'evidence/restart_v06/preparation/protocol.json').read_text())

    def analyze(self, report=None, spec=None, numerical_review=None):
        return analyze_preparation_readiness(self.assigned, report if report is not None else self.prep,
            spec or PreparationReadinessSpec(), numerical_review=numerical_review)

    def test_actual_short_record_retains_budget_and_sampling_failures(self):
        result = self.analyze()
        self.assertEqual(result['status'], 'requires_preparation')
        self.assertIn('minimizer_policy', result['flags'])
        self.assertIn('npt_density_sampling', result['flags'])
        self.assertIn('starting_configuration_numerical_review', result['flags'])
        self.assertGreater(result['force_summary']['max_atom_force_norm_kj_mol_angstrom'], .01)
        self.assertTrue(result['gates']['finite_final_forces_and_requested_development_bounds'])
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertEqual(result['minimizer_force_norms']['status'], 'not_recorded_by_existing_protocol')

    def test_complete_numerical_review_applies_only_to_start_not_final_state(self):
        result = self.analyze(numerical_review=numerical_fixture(self.assigned, self.prep))
        self.assertTrue(result['gates']['starting_configuration_numerical_review'])
        self.assertTrue(result['numerical_qualification_scope']['full_stage_comparison_recomputed'])
        self.assertIn('not_established', result['numerical_qualification_scope']['final_configuration_numerics'])
        self.assertFalse(result['production_ready'])

    def test_numerical_review_forgery_missing_stages_and_wrong_controls_fail(self):
        base = numerical_fixture(self.assigned, self.prep)
        cases = []
        changed = deepcopy(base); changed['stages'].pop(); cases.append(changed)
        changed = deepcopy(base); changed['review']['numerical_status'] = 'requires_review'; cases.append(changed)
        changed = deepcopy(base); changed['selected_settings']['cutoff_angstrom'] = 17.; cases.append(changed)
        changed = deepcopy(base); changed['assignment_canonical_sha256'] = 'a'*64; cases.append(changed)
        changed = deepcopy(base); changed['geometry_canonical_sha256'] = 'b'*64; cases.append(changed)
        changed = deepcopy(base); changed['stages'][0]['record']['forces_kj_mol_angstrom'] = []; cases.append(changed)
        for review in cases:
            with self.subTest(keys=review.keys()), self.assertRaises(ValueError):
                self.analyze(numerical_review=review)

    def test_changed_producer_engine_is_not_adopted(self):
        report = deepcopy(self.prep)
        review = numerical_fixture(self.assigned, report)
        report['engine_run']['engine_identity']['entrypoint_sha256'] = 'a'*64
        with self.assertRaisesRegex(ValueError, 'engine identities'):
            self.analyze(report, numerical_review=review)

    def test_source_identity_or_unit_proof_mismatch_is_rejected(self):
        for case in ('assignment', 'units', 'step', 'duration', 'native_density', 'raw_step', 'boundary'):
            report = deepcopy(self.prep)
            if case == 'assignment': report['assignment_canonical_sha256'] = 'a'*64
            if case == 'units': report['independent_output_unit_check']['header'] = 'units = metal'
            if case == 'step': report['final_state']['step'] += 1
            if case == 'duration': report['stages']['npt']['elapsed_time_fs'] += 1
            if case == 'native_density': report['engine_run']['thermo_blocks'][-1]['rows'][-1][4] *= 2
            if case == 'raw_step': report['engine_run']['thermo_blocks'][-1]['rows'][-1][0] -= 1
            if case == 'boundary': report['stages']['npt']['thermo_rows'][0]['temperature_k'] += 1
            with self.subTest(case=case), self.assertRaises(ValueError): self.analyze(report)

    def test_summary_missing_final_atoms_and_tampered_images_forces_rejected(self):
        for case in ('summary', 'image', 'nan', 'shape'):
            report = deepcopy(self.prep)
            if case == 'summary': del report['final_state']['atoms']
            if case == 'image': report['final_state']['atoms'][0]['image_flags'][0] += 1
            if case == 'nan': report['final_state']['atoms'][0]['forces_kj_mol_angstrom'][0] = float('nan')
            if case == 'shape': report['final_state']['atoms'][0]['forces_kj_mol_angstrom'] = [1, 2]
            with self.subTest(case=case), self.assertRaises(ValueError): self.analyze(report)

    def test_minimizer_flag_is_not_inferred_from_small_thermal_forces(self):
        report = deepcopy(self.prep)
        for atom in report['final_state']['atoms']: atom['forces_kj_mol_angstrom'] = [0., 0., 0.]
        result = self.analyze(report)
        self.assertFalse(result['gates']['minimizer_policy'])
        self.assertEqual(result['force_summary']['rms_atom_force_norm_kj_mol_angstrom'], 0)
        report['stages']['minimization']['status'] = 'force_tolerance_met'
        with self.assertRaisesRegex(ValueError, 'criterion'): self.analyze(report)

    def test_optional_force_screen_uses_vector_norm_and_explicit_policy(self):
        report = deepcopy(self.prep)
        for atom in report['final_state']['atoms']: atom['forces_kj_mol_angstrom'] = [3., 4., 0.]
        result = self.analyze(report, PreparationReadinessSpec(max_force_norm_kj_mol_angstrom=4.9,
            max_rms_force_kj_mol_angstrom=4.9))
        self.assertEqual(result['force_summary']['max_atom_force_norm_kj_mol_angstrom'], 5.)
        self.assertEqual(result['force_summary']['rms_atom_force_norm_kj_mol_angstrom'], 5.)
        self.assertFalse(result['gates']['finite_final_forces_and_requested_development_bounds'])

    def test_failed_engine_outcome_retained_without_fabricated_statistics(self):
        report = deepcopy(self.prep)
        report.update(status='failed', failure_phase='engine', failure_reason='explicit synthetic timeout')
        report.pop('final_state'); report['engine_run'] = None
        result = self.analyze(report)
        self.assertEqual(result['status'], 'requires_preparation')
        self.assertIn('timeout', result['failure_reason'])
        self.assertIsNone(result['force_summary'])
        self.assertIsNone(result['density_sampling_review'])

    def test_passing_local_diagnostics_still_cannot_certify_liquid_or_production(self):
        report = synthetic_long_record(self.assigned, self.prep)
        sampling = DensitySamplingSpec(timestep_fs=.25, min_duration_ps=.5, block_size=32)
        spec = PreparationReadinessSpec(require_minimizer_force_convergence=False,
            max_nematic_order=1., sampling_spec=sampling)
        result = self.analyze(report, spec, numerical_fixture(self.assigned, report))
        self.assertEqual(result['status'], 'diagnostics_passed_but_unqualified', result['flags'])
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['phase_identification'], 'not_performed')
        self.assertEqual(result['independent_cell_count'], 1)


if __name__ == '__main__': unittest.main()
