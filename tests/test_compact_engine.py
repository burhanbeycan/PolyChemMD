"""Compact representation equivalence and bounded single-point contracts.

Synthetic engine fixtures below test auditing, never stand in for real MD.
"""

from pathlib import Path
import json
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.builders.packing import pack_cell
from polychemmd.engines.lammps import LammpsError, LammpsRunResult
from polychemmd.engines.periodic_alkane import (
    _COMPONENTS, _THERMO, PeriodicSinglepointError, build_periodic_alkane_input,
    prepare_periodic_alkane_singlepoint, run_periodic_alkane_singlepoint,
)
from polychemmd.forcefields.periodic import (
    assign_compact_periodic_alkanes, assign_periodic_alkanes,
)
from polychemmd.forcefields.periodic_reference import evaluate_periodic_alkanes
from test_periodic_assignment import analytic_ethane, packed_ethanes
from test_periodic_lammps import dump_text


ROOT = Path(__file__).resolve().parents[1]


class CompactEngineEquivalenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cell = packed_ethanes()
        cls.legacy = assign_periodic_alkanes(cls.cell)
        cls.compact = assign_compact_periodic_alkanes(cls.cell)

    def test_small_legacy_and_compact_export_identical_scientific_inputs(self):
        for settings in ({}, {"cutoff_angstrom": 11., "reciprocal_nmax": 20}):
            self.assertEqual(build_periodic_alkane_input(self.legacy, **settings),
                             build_periodic_alkane_input(self.compact, **settings))

    def test_actual_decane_fixture_exports_identically_without_global_pairs(self):
        cell = json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        compact = assign_compact_periodic_alkanes(cell)
        self.assertNotIn('nonbonded_pairs', compact.to_dict())
        self.assertEqual(build_periodic_alkane_input(assign_periodic_alkanes(cell)),
                         build_periodic_alkane_input(compact))

    def test_independent_small_reference_has_identical_energies_and_forces(self):
        legacy = evaluate_periodic_alkanes(self.legacy, reciprocal_nmax=3)
        compact = evaluate_periodic_alkanes(self.compact, reciprocal_nmax=3)
        self.assertEqual(legacy['energy_components_kj_mol'], compact['energy_components_kj_mol'])
        self.assertEqual(legacy['forces_kj_mol_angstrom'], compact['forces_kj_mol_angstrom'])
        self.assertEqual(legacy['pair_counts'], compact['pair_counts'])
        self.assertFalse(compact['pair_work']['global_pair_array_materialized_by_kernel'])
        self.assertNotIn('nonbonded_pairs', self.compact.to_dict())

    def test_reference_pair_budget_fails_before_bonded_or_reciprocal_work(self):
        with patch('polychemmd.forcefields.periodic_reference._differentiate') as bonded, \
             patch('polychemmd.forcefields.periodic_reference._reciprocal_ewald') as reciprocal:
            for budget in (119, True, 0, 130817):
                with self.subTest(budget=budget), self.assertRaises(ValueError):
                    evaluate_periodic_alkanes(self.compact, reciprocal_nmax=2,
                                             max_pair_evaluations=budget)
            bonded.assert_not_called()
            reciprocal.assert_not_called()

    def test_large_compact_export_succeeds_but_reference_remains_explicitly_small(self):
        cell = pack_cell([analytic_ethane()]*65, box_lengths_angstrom=(80.,80.,80.),
                         seed=9091, min_distance_angstrom=1.5).to_dict()
        assigned = assign_compact_periodic_alkanes(cell)
        script, data = build_periodic_alkane_input(assigned)
        self.assertIn('520 atoms', data)
        self.assertIn('kspace_style ewald', script)
        with patch('polychemmd.forcefields.periodic_reference._numeric_pairs') as pairs, \
             patch('polychemmd.forcefields.periodic_reference._reciprocal_ewald') as reciprocal:
            with self.assertRaisesRegex(ValueError, '1 to 512 atoms'):
                evaluate_periodic_alkanes(assigned, reciprocal_nmax=2)
            pairs.assert_not_called()
            reciprocal.assert_not_called()


class SinglepointExecutionContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assigned = assign_compact_periodic_alkanes(packed_ethanes())
        cls.record = cls.assigned.to_dict()

    def fixture_run(self, spec, *, mutation=None):
        # Deliberately synthetic process output, not a simulation.
        spec.run_directory.mkdir()
        text = dump_text(self.record['atoms'], self.record['box']['lengths_angstrom'])
        if mutation == 'image':
            lines=text.splitlines(); row=lines[-1].split(); row[5]=str(int(row[5])+1)
            lines[-1]=' '.join(row); text='\n'.join(lines)+'\n'
        (spec.run_directory/'forces.dump').write_text(text)
        (spec.run_directory/'singlepoint.data').write_text(
            'Synthetic native fixture, units = '+('lj' if mutation == 'units' else 'real')+'\n')
        thermo = [0., 0., 1., 2., 3., 4., 5., 6., 21.]
        if mutation == 'sum': thermo[-1] = 20.
        if mutation == 'overflow': thermo[2] = 1e308
        return LammpsRunResult({'simulation_performed': True,
            'thermo_blocks': [{'columns': list(_THERMO), 'rows': [thermo]}],
            'engine_version_banner': 'synthetic process fixture, not LAMMPS'})

    def test_prepared_readback_is_reviewable_without_execution_or_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)/'not-created'
            spec = prepare_periodic_alkane_singlepoint(self.assigned, run_directory=target,
                                                       native_readback=True)
            self.assertFalse(target.exists())
            self.assertIn('write_data singlepoint.data', spec.input_script)
            self.assertIn('singlepoint.data', spec.expected_artifacts)
            self.assertEqual(spec.units_style, 'real')

    def test_singlepoint_does_not_silently_invoke_a_large_reference_or_claim_parity(self):
        with tempfile.TemporaryDirectory() as temp:
            target = Path(temp)/'run'
            with patch('polychemmd.engines.periodic_alkane.run_lammps', side_effect=self.fixture_run), \
                 patch('polychemmd.forcefields.periodic_reference.evaluate_periodic_alkanes') as reference:
                report = run_periodic_alkane_singlepoint(self.assigned, run_directory=target)
                reference.assert_not_called()
            self.assertEqual(json.loads((target/'periodic_singlepoint_record.json').read_text()), report)
        self.assertEqual(report['status'], 'completed_singlepoint')
        self.assertEqual(report['energy_components_kj_mol']['total'], 21.*4.184)
        self.assertEqual(len(report['forces_kj_mol_angstrom']), 16)
        self.assertEqual(report['independent_reference_status'], 'not_performed')
        self.assertEqual(report['numerical_convergence_status'], 'not_established')
        self.assertEqual(report['equilibrium_status'], 'not_assessed')
        self.assertFalse(report['bulk_production_model_validated'])
        self.assertEqual(report['md_integration_steps'], 0)

    def test_failed_identity_units_sum_and_canonical_overflow_are_retained(self):
        for mutation in ('image', 'units', 'sum', 'overflow'):
            with self.subTest(mutation=mutation), tempfile.TemporaryDirectory() as temp:
                target=Path(temp)/'run'
                def fake(spec): return self.fixture_run(spec, mutation=mutation)
                with patch('polychemmd.engines.periodic_alkane.run_lammps', side_effect=fake):
                    with self.assertRaises(PeriodicSinglepointError) as raised:
                        run_periodic_alkane_singlepoint(self.assigned, run_directory=target)
                report=json.loads((target/'periodic_singlepoint_record.json').read_text())
                self.assertEqual(report['status'], 'failed')
                self.assertEqual(report['failure_phase'], 'singlepoint_output_audit')
                self.assertEqual(raised.exception.record, report)

    def test_engine_failure_and_unavailable_cases_do_not_fabricate_observables(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'absent'
            with patch('polychemmd.engines.periodic_alkane.run_lammps',
                       side_effect=LammpsError('unavailable synthetic fixture')):
                with self.assertRaises(PeriodicSinglepointError) as raised:
                    run_periodic_alkane_singlepoint(self.assigned, run_directory=target)
            self.assertFalse(target.exists())
            self.assertFalse(raised.exception.record['simulation_performed'])
            self.assertNotIn('energy_components_kj_mol', raised.exception.record)

    def test_timeout_budget_is_checked_before_engine_or_directory(self):
        with patch('polychemmd.engines.periodic_alkane.run_lammps') as engine:
            for timeout in (True, 0., -1., float('nan'), 120.0001):
                with self.subTest(timeout=timeout), self.assertRaises(ValueError):
                    run_periodic_alkane_singlepoint(self.assigned, run_directory='not-created',
                                                    timeout_seconds=timeout)
            engine.assert_not_called()


if __name__ == '__main__':
    unittest.main()
