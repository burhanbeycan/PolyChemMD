"""Hand-controlled numerical sweeps and synthetic execution failure contracts.

These fixtures are analytic/synthetic records, not computed material evidence.
"""
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.builders.packing import pack_cell
from polychemmd.forcefields.periodic import assign_compact_periodic_alkanes
from polychemmd.engines.lammps import LammpsError
from polychemmd.engines.periodic_thermo import plan_periodic_alkane_thermo_singlepoint
from polychemmd.protocols.numerics import (
    NumericalQualificationError, NumericalQualificationSpec, compare_numerical_records,
    plan_numerical_qualification, run_numerical_qualification,
)
from test_periodic_assignment import analytic_ethane


def synthetic_records(spec=None, count=2, ids=None):
    spec=spec or NumericalQualificationSpec()
    settings=[{'cutoff_angstrom':cut,'alpha_per_angstrom':spec.ewald_alpha_per_angstrom,
               'reciprocal_nmax':spec.reciprocal_nmax,'lj_tail_correction':False,'surface_boundary':'tin_foil','relative_engine_accuracy_request':1e-8,'realspace_coulomb_table':False}
              for cut in spec.cutoffs_angstrom]
    settings.append(dict(settings[-1],reciprocal_nmax=spec.refined_reciprocal_nmax))
    values=[100.,90.,88.,87.98,87.97,87.97]
    pressures=[200.,150.,125.,124.,123.8,123.8]
    forces=[10.,9.,8.,7.99,7.98,7.98]
    result=[]
    for index,setting in enumerate(settings):
        energy={'bond':1.,'angle':2.,'dihedral':3.,'lj':values[index],
                'coulomb_realspace':5.,'coulomb_kspace':6.}
        energy['total']=sum(energy.values())
        pressure=pressures[index]
        result.append({'kind':'periodic_alkane_static_thermo_singlepoint','status':'completed_static_thermo_singlepoint','simulation_performed':True,
            'md_integration_steps':0,'trajectory_production_performed':False,
            'pressure_interpretation':'static_virial_not_equilibrium_average',
            'initial_velocities':'zeroed_and_output_verified','temperature_k':0.,'kinetic_energy_kj_mol':0.,
            'atom_coordinate_molecule_charge_image_identity':'passed','pressure_compute':'compute polychemmd_virial all pressure NULL virial',
            'independent_output_unit_check':{'confirmed_units_style':'real','header':'synthetic fixture units = real'},
            'assignment_canonical_sha256':'a'*64,'geometry_canonical_sha256':'b'*64,
            'model_id':'explicit_synthetic_equation_fixture','atom_count':count,
            'atom_ids':ids or [f'atom:{i}' for i in range(count)],
            'energy_units':'kJ/mol','force_units':'kJ/(mol angstrom)','pressure_units':'bar',
            'volume_units':'angstrom^3','density_units':'g/cm^3',
            'energy_components_kj_mol':energy,'ewald_settings':setting,
            'forces_kj_mol_angstrom':[[forces[index],0.,0.] for _ in range(count)],
            'pressure_bar':pressure,'pressure_tensor_bar':[[pressure if i==j else 0. for j in range(3)] for i in range(3)],
            'volume_angstrom3':64000.,'density_g_cm3':.1,
            'engine_run':{'status':'completed','simulation_performed':True,'engine_identity':{'kind':'synthetic_test_engine_not_LAMMPS'}}})
    return result


class NumericalSpecTests(unittest.TestCase):
    def test_exact_roundtrip_defaults_and_immutable_cutoffs(self):
        values=[10,12,14,16,17]
        spec=NumericalQualificationSpec.from_mapping({'cutoffs_angstrom':values})
        values[-1]=99
        self.assertEqual(spec.cutoffs_angstrom[-1],17)
        self.assertEqual(NumericalQualificationSpec.from_mapping(spec.to_dict()),spec)
        self.assertEqual(spec.max_engine_runs,6)

    def test_bad_limits_booleans_nonfinite_and_unknown_fields_rejected(self):
        for data in ({'unknown':1},{'cutoffs_angstrom':[10,10,12]}, {'cutoffs_angstrom':[10,12]},
            {'cutoffs_angstrom':[10,True,17]},{'cutoff_pressure_limit_bar':float('nan')},
            {'total_wall_seconds':241},{'per_run_timeout_seconds':121},
            {'per_run_timeout_seconds':60,'total_wall_seconds':30},
            {'max_engine_runs':5},{'min_consecutive_cutoff_pairs':1},
            {'reciprocal_nmax':20},{'refined_reciprocal_nmax':16}):
            with self.subTest(data=data),self.assertRaises(ValueError):
                NumericalQualificationSpec.from_mapping(data)


class NumericalReviewTests(unittest.TestCase):
    def test_two_finest_pairs_pass_without_erasing_coarser_failures(self):
        review=compare_numerical_records(synthetic_records(),NumericalQualificationSpec())
        self.assertEqual(review['numerical_status'],'fixed_geometry_diagnostics_passed')
        self.assertEqual([r['status'] for r in review['cutoff_comparisons']],
                         ['requires_review','requires_review','diagnostics_passed','diagnostics_passed'])
        self.assertEqual(review['selected_settings'],{'cutoff_angstrom':17.,'ewald_alpha_per_angstrom':.3,'reciprocal_nmax':20})
        self.assertAlmostEqual(review['cutoff_comparisons'][-2]['absolute_total_energy_difference_per_atom_kj_mol'],.01)
        self.assertFalse(review['production_ready'])
        self.assertEqual(review['equilibrium_status'],'not_assessed')
        self.assertEqual(review['bulk_property_convergence'],'not_established')

    def test_one_finest_passing_pair_is_insufficient(self):
        records=synthetic_records()
        records[3]['pressure_bar']=150.
        records[3]['pressure_tensor_bar']=[[150. if i==j else 0. for j in range(3)] for i in range(3)]
        review=compare_numerical_records(records,NumericalQualificationSpec())
        self.assertEqual(review['numerical_status'],'requires_review')
        self.assertIsNone(review['selected_settings'])

    def test_tensor_anisotropy_cannot_hide_behind_scalar_trace(self):
        records=synthetic_records(); p=records[4]['pressure_bar']
        records[4]['pressure_tensor_bar']=[[p+100,0.,0.],[0.,p-100,0.],[0.,0.,p]]
        review=compare_numerical_records(records,NumericalQualificationSpec())
        last=review['cutoff_comparisons'][-1]
        self.assertTrue(last['gates']['scalar_static_pressure'])
        self.assertFalse(last['gates']['static_pressure_tensor'])
        self.assertEqual(review['numerical_status'],'requires_review')

    def test_reciprocal_change_retains_failure_and_bonded_invariant_violation(self):
        for field in ('coulomb_kspace','bond'):
            records=synthetic_records();records[-1]['energy_components_kj_mol'][field]+=1
            records[-1]['energy_components_kj_mol']['total']+=1
            review=compare_numerical_records(records,NumericalQualificationSpec())
            self.assertEqual(review['numerical_status'],'requires_review')
            self.assertIsNone(review['selected_settings'])
            if field=='bond':self.assertFalse(review['reciprocal_comparison']['gates']['unchanged_components'])

    def test_changed_source_order_model_engine_controls_or_units_rejected(self):
        changes=[('geometry_canonical_sha256','c'*64),('assignment_canonical_sha256','c'*64),
                 ('atom_ids',['atom:1','atom:0']),('model_id','other'),('pressure_units','atm')]
        for field,value in changes:
            records=synthetic_records();records[2][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):
                compare_numerical_records(records,NumericalQualificationSpec())
        for changed in ('engine','tail','cutoff'):
            records=synthetic_records()
            if changed=='engine':records[1]['engine_run']['engine_identity']['kind']='other_engine'
            if changed=='tail':records[1]['ewald_settings']['lj_tail_correction']=True
            if changed=='cutoff':records[1]['ewald_settings']['cutoff_angstrom']=11.
            with self.subTest(changed=changed),self.assertRaises(ValueError):
                compare_numerical_records(records,NumericalQualificationSpec())

    def test_missing_and_failed_records_cannot_be_selected_away(self):
        for records in (synthetic_records()[-2:],synthetic_records()[:-1],[]):
            # A truncated ordered prefix is insufficient; a cherry-picked suffix
            # fails its requested-settings identity check rather than passing.
            try:review=compare_numerical_records(records,NumericalQualificationSpec())
            except ValueError:continue
            self.assertEqual(review['numerical_status'],'incomplete_records')
            self.assertIsNone(review['selected_settings'])
        records=synthetic_records();records[1]={'status':'failed','failure_reason':'retained synthetic timeout'}
        review=compare_numerical_records(records,NumericalQualificationSpec())
        self.assertEqual(review['numerical_status'],'incomplete_records')
        self.assertIn('timeout',review['case_outcomes'][1]['failure_reason'])

    def test_failed_partial_record_still_requires_finite_json_and_valid_status_reason(self):
        for record in ({'status':'failed','failure_reason':float('nan')},
                       {'status':7,'failure_reason':'bad type'},
                       {'status':'failed','failure_reason':123},
                       {'status':'failed','engine_run':{'elapsed_seconds':float('inf')}}):
            with self.subTest(record=record),self.assertRaises(ValueError):
                compare_numerical_records([record],NumericalQualificationSpec())

    def test_nonmonotone_changes_are_reported_without_fabricating_error_bounds(self):
        records=synthetic_records();records[4]['energy_components_kj_mol']['lj']=87.95
        records[4]['energy_components_kj_mol']['total']=104.95
        records[5]['energy_components_kj_mol']=deepcopy(records[4]['energy_components_kj_mol'])
        review=compare_numerical_records(records,NumericalQualificationSpec())
        self.assertFalse(review['absolute_cutoff_energy_changes_monotone_nonincreasing'])
        self.assertTrue(review['nonmonotonic_changes_are_retained'])
        self.assertIn('no exact-model error',review['cutoff_comparisons'][-1]['interpretation'])

    def test_nonfinite_vectors_scalar_trace_or_component_sum_rejected(self):
        for kind in ('force','sum','trace','boolean'):
            records=synthetic_records()
            if kind=='force':records[-1]['forces_kj_mol_angstrom'][0][0]=float('inf')
            if kind=='sum':records[-1]['energy_components_kj_mol']['total']+=10
            if kind=='trace':records[-1]['pressure_tensor_bar'][0][0]+=10
            if kind=='boolean':records[-1]['pressure_bar']=True
            with self.subTest(kind=kind),self.assertRaises(ValueError):
                compare_numerical_records(records,NumericalQualificationSpec())

    def test_false_static_labels_do_not_hide_temperature_kinetic_or_native_unit_changes(self):
        for field,value in (('temperature_k',1000.),('kinetic_energy_kj_mol',999.),
                            ('md_integration_steps',False),('kind','other')):
            records=synthetic_records();records[2][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):
                compare_numerical_records(records,NumericalQualificationSpec())
        records=synthetic_records();records[2]['independent_output_unit_check']['header']='units = lj'
        with self.assertRaises(ValueError):compare_numerical_records(records,NumericalQualificationSpec())
        records=synthetic_records();records[2]['engine_run']['status']='failed'
        with self.assertRaises(ValueError):compare_numerical_records(records,NumericalQualificationSpec())

    def test_unrepresentable_finite_difference_fails_explicitly(self):
        records=synthetic_records()
        for index,value in ((3,-1e308),(4,1e308)):
            records[index]['energy_components_kj_mol']['lj']=value
            records[index]['energy_components_kj_mol']['total']=value
        with self.assertRaisesRegex(ValueError,'difference'):
            compare_numerical_records(records,NumericalQualificationSpec())


class NumericalExecutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assigned=assign_compact_periodic_alkanes(pack_cell([analytic_ethane()]*2,
            box_lengths_angstrom=(40.,40.,40.),seed=985,min_distance_angstrom=1.5).to_dict())

    def fake_runner(self, assigned, **kwargs):
        directory=Path(kwargs['run_directory']); root=directory.parent
        manifest=json.loads((root/'numerical_qualification.json').read_text())
        self.assertEqual(manifest['status'],'running')
        self.assertTrue((root/'common/system.data').is_file())
        plan=plan_periodic_alkane_thermo_singlepoint(assigned,**kwargs)
        requests=manifest['stages'];index=next(i for i,item in enumerate(requests) if item['run_directory']==str(directory))
        result=synthetic_records(count=len(plan['atom_ids']),ids=plan['atom_ids'])[index]
        for name in ('assignment_canonical_sha256','geometry_canonical_sha256','model_id','atom_count','input_files'):
            result[name]=deepcopy(plan[name])
        # Match native fixed control metadata to the actual reviewable plan.
        result['ewald_settings']=deepcopy(plan['ewald_settings'])
        directory.mkdir();return result

    def test_plan_validates_entire_sweep_without_engine_or_files(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'absent'
            with patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint') as run:
                plan=plan_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
                run.assert_not_called()
            self.assertFalse(target.exists())
            self.assertEqual(len(plan['stages']),6)
            self.assertFalse(plan['simulation_performed'])
            self.assertEqual(hashlib.sha256(plan['staged_inputs']['common/system.data'].encode()).hexdigest(),plan['shared_input']['sha256'])

    def test_oversized_cutoff_rejected_before_engine_or_directory(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'absent'
            spec=NumericalQualificationSpec(cutoffs_angstrom=(10,12,14,16,20))
            with self.assertRaisesRegex(ValueError,'half'):
                plan_numerical_qualification(self.assigned,spec,run_directory=target)
            self.assertFalse(target.exists())

    def test_complete_declared_execution_preserves_script_and_common_input_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            with patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint',side_effect=self.fake_runner) as run:
                report=run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
            self.assertEqual(run.call_count,6)
            self.assertEqual(report['status'],'completed_numerical_review')
            self.assertEqual(report['numerical_status'],'fixed_geometry_diagnostics_passed')
            self.assertEqual(report['engine_identity'],{'kind':'synthetic_test_engine_not_LAMMPS'})
            self.assertEqual(json.loads((target/'numerical_qualification.json').read_text()),report)
            self.assertNotIn('staged_inputs',report)

    def test_failed_process_is_retained_without_missing_case_success(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            error=LammpsError('synthetic unavailable engine',record={'status':'failed','simulation_performed':False,'failure_reason':'unavailable'})
            with patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint',side_effect=error):
                with self.assertRaises(NumericalQualificationError) as caught:
                    run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
            report=caught.exception.record
            self.assertEqual(report['status'],'failed')
            self.assertFalse(report['simulation_performed'])
            self.assertEqual(report['engine_requests_started'],1)
            self.assertEqual(report['engine_execution_attempts_reported'],0)
            self.assertEqual(report['partial_review']['numerical_status'],'incomplete_records')
            self.assertEqual(report['stages'][0]['status'],'failed')
            self.assertEqual(json.loads((target/'numerical_qualification.json').read_text()),report)

    def test_changed_input_hash_fails_and_retains_partial_case(self):
        def changed(assigned,**kwargs):
            result=self.fake_runner(assigned,**kwargs);result['input_files']['system.data_sha256']='f'*64;return result
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            with patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint',side_effect=changed):
                with self.assertRaisesRegex(NumericalQualificationError,'hashes'):
                    run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
            report=json.loads((target/'numerical_qualification.json').read_text())
            self.assertEqual(report['status'],'failed')
            self.assertEqual(report['engine_requests_started'],1)

    def test_aggregate_budget_stops_before_a_second_process(self):
        clock=[0.]
        def slow(assigned,**kwargs):
            result=self.fake_runner(assigned,**kwargs);clock[0]=181.;return result
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            with patch('polychemmd.protocols.numerics.time.monotonic',side_effect=lambda:clock[0]), \
                 patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint',side_effect=slow) as run:
                with self.assertRaisesRegex(NumericalQualificationError,'budget'):
                    run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
            self.assertEqual(run.call_count,1)
            report=json.loads((target/'numerical_qualification.json').read_text())
            self.assertEqual(report['measured_total_wall_seconds'],181.)

    def test_postprocessing_and_final_serialization_overruns_cannot_return_success(self):
        from polychemmd.protocols import numerics
        original_compare=numerics.compare_numerical_records
        original_persist=numerics._persist
        for phase in ('postprocess','serialize'):
            clock=[0.]
            def compare(records,spec):
                result=original_compare(records,spec)
                if phase=='postprocess' and len(records)==6:clock[0]=181.
                return result
            def persist(directory,report):
                original_persist(directory,report)
                if phase=='serialize' and report['status']=='completed_numerical_review':clock[0]=181.
            with self.subTest(phase=phase),tempfile.TemporaryDirectory() as temp:
                target=Path(temp)/'run'
                with patch('polychemmd.protocols.numerics.time.monotonic',side_effect=lambda:clock[0]), \
                     patch('polychemmd.protocols.numerics.compare_numerical_records',side_effect=compare), \
                     patch('polychemmd.protocols.numerics._persist',side_effect=persist), \
                     patch('polychemmd.engines.periodic_thermo.run_periodic_alkane_thermo_singlepoint',side_effect=self.fake_runner):
                    with self.assertRaises(NumericalQualificationError):
                        run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
                report=json.loads((target/'numerical_qualification.json').read_text())
                self.assertEqual(report['status'],'failed')
                self.assertIsNone(report['selected_settings'])

    def test_stale_directory_is_never_mutated(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run';target.mkdir();sentinel=target/'numerical_qualification.json';sentinel.write_text('old')
            with self.assertRaises(NumericalQualificationError):
                run_numerical_qualification(self.assigned,NumericalQualificationSpec(),run_directory=target)
            self.assertEqual(sentinel.read_text(),'old')


if __name__=='__main__':unittest.main()
