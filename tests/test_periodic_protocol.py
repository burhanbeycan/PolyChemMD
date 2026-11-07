"""Protocol contracts use explicit synthetic frames, not scientific MD evidence."""

import json
import math
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

from polychemmd.analysis.density import density_from_mass_volume
from polychemmd.engines.lammps import LammpsError, LammpsRunResult
from polychemmd.protocols.alkane import (
    BAR_PER_ATM, PeriodicAlkaneProtocolSpec, PeriodicProtocolError,
    _DUMP_COLUMNS, _MARKERS, _THERMO_COLUMNS, _audit_outputs, _minimizer_status,
    parse_periodic_final_dump, plan_periodic_alkane_protocol,
    run_periodic_alkane_protocol,
)


def fixture_assignment():
    return {'model_id':'SYNTHETIC_PROTOCOL_PROCESS_FIXTURE_NOT_A_FORCE_FIELD','interactions':{'bonds':[]},
            'atoms':[{'cell_atom_id':'molecule:0/a','molecule_index':0,'charge_e':-.06,'mass_g_mol':1.008},
                     {'cell_atom_id':'molecule:0/b','molecule_index':0,'charge_e':.06,'mass_g_mol':1.008}]}


def fixture_dump(step=5):
    return '\n'.join(['ITEM: TIMESTEP',str(step),'ITEM: NUMBER OF ATOMS','2',
                       'ITEM: BOX BOUNDS pp pp pp','0 24','0 24','0 24',
                       'ITEM: ATOMS '+' '.join(_DUMP_COLUMNS),
                       '2 1 .06 2 3 4 0 0 0 -1 -2 -3',
                       '1 1 -.06 1 2 3 -1 0 1 1 2 3'])+'\n'


def setup_fixture(*args,**kwargs):
    return 'clear\nunits real\n# Synthetic setup contract fixture; no executable is run.\n','SYNTHETIC_DATA'


class ProtocolSpecTests(unittest.TestCase):
    def test_roundtrip_defaults_and_input_copy(self):
        spec=PeriodicAlkaneProtocolSpec()
        self.assertEqual(spec.nvt_steps+spec.npt_steps,200)
        self.assertEqual(PeriodicAlkaneProtocolSpec.from_mapping(spec.to_dict()),spec)
        out=spec.to_dict();out['temperature_k']=500
        self.assertEqual(spec.temperature_k,300)

    def test_invalid_scientific_numbers_and_integer_budgets(self):
        for field in ('temperature_k','pressure_bar','timestep_fs','timeout_seconds',
                      'cutoff_angstrom','ewald_alpha_per_angstrom'):
            for value in (False,0,-1,float('nan'),float('inf'),10**1000):
                with self.subTest(field=field,value_type=type(value)),self.assertRaises(ValueError):
                    PeriodicAlkaneProtocolSpec(**{field:value})
        for field in ('minimize_max_iterations','minimize_max_evaluations','nvt_steps',
                      'npt_steps','velocity_seed','reciprocal_nmax'):
            for value in (True,0,-1,1.5):
                with self.subTest(field=field,value=value),self.assertRaises(ValueError):
                    PeriodicAlkaneProtocolSpec(**{field:value})

    def test_upper_resource_bounds_and_unknown_fields(self):
        for overrides in ({'timestep_fs':2},{'timeout_seconds':121},
                          {'nvt_steps':6000,'npt_steps':6000},{'velocity_seed':900000001},
                          {'minimize_max_iterations':5001},{'minimize_max_evaluations':50001},
                          {'reciprocal_nmax':21}):
            with self.subTest(overrides=overrides),self.assertRaises(ValueError):
                PeriodicAlkaneProtocolSpec(**overrides)
        with self.assertRaises(ValueError):
            PeriodicAlkaneProtocolSpec.from_mapping({'assume_equilibrium':True})


class ProtocolPlanningTests(unittest.TestCase):
    def test_plan_creates_nothing_and_records_native_unit_conversions(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)/'not-created'
            assigned=SimpleNamespace(to_dict=fixture_assignment)
            spec=PeriodicAlkaneProtocolSpec(pressure_bar=BAR_PER_ATM)
            with patch('polychemmd.engines.periodic_alkane.build_periodic_alkane_setup',side_effect=setup_fixture),\
                 patch('polychemmd.protocols.alkane.run_lammps') as run:
                plan=plan_periodic_alkane_protocol(assigned,spec,run_directory=directory,executable='missing-engine')
                run.assert_not_called()
            self.assertFalse(directory.exists())
        self.assertFalse(plan['simulation_performed'])
        self.assertFalse(plan['engine_plan']['simulation_performed'])
        self.assertEqual(plan['scientific_controls']['engine_pressure_target_atm'],1)
        self.assertAlmostEqual(plan['scientific_controls']['minimize_force_tolerance_kj_mol_angstrom'],4.184e-6)
        self.assertIn('minimize 1e-6 1e-6 200 1000',plan['input_script'])
        self.assertIn('iso 1 1 1000.0',plan['input_script'])
        self.assertIn('(lx<=20)||(ly<=20)||(lz<=20)',plan['input_script'])
        self.assertIn('halt 1 v_polychemmd_box_invalid != 0 error hard',plan['input_script'])
        self.assertIn('halt 1 bondmax >= ',plan['input_script'])
        self.assertEqual(plan['engine_plan']['units_style_confidence'],'unresolved')
        self.assertEqual(plan['engine_plan']['declared_units_style'],'real')

    def test_minimizer_stop_status_is_separate_from_process_completion(self):
        for criterion,status,force in [('max iterations','budget_exhausted',False),
                                        ('max force evaluations','budget_exhausted',False),
                                        ('force tolerance','force_tolerance_met',True),
                                        ('energy tolerance','energy_tolerance_met_force_unverified',False),
                                        ('linesearch alpha is zero','inconclusive',False)]:
            record=_minimizer_status(f'Stopping criterion = {criterion}\nIterations, force evaluations = 5 17\n')
            self.assertEqual(record['status'],status)
            self.assertEqual(record['force_tolerance_met'],force)
            self.assertEqual((record['iterations'],record['force_evaluations']),(5,17))
        self.assertEqual(_minimizer_status('no stopping information')['status'],'inconclusive')


class FinalSnapshotTests(unittest.TestCase):
    def test_source_identity_charge_image_flags_and_force_units(self):
        frame=parse_periodic_final_dump(fixture_dump(),fixture_assignment(),5)
        self.assertEqual([a['cell_atom_id'] for a in frame['atoms']],['molecule:0/a','molecule:0/b'])
        self.assertEqual(frame['atoms'][0]['coordinates_angstrom'],[-23.,2.,27.])
        self.assertEqual(frame['atoms'][0]['forces_kj_mol_angstrom'],[4.184,8.368,12.552])
        self.assertEqual(frame['box']['lengths_angstrom'],[24.,24.,24.])

    def test_malformed_final_frames_fail(self):
        good=fixture_dump()
        cases=[good+good,good.replace('2 1 .06','1 1 .06'),good.replace('2 1 .06','2 2 .06'),
               good.replace('2 1 .06','2 1 .07'),good.replace('2 3 4','nan 3 4'),
               good.replace('2 3 4','24 3 4'),good.replace('0 0 0 -1','0.5 0 0 -1'),
               good.replace('BOX BOUNDS pp pp pp','BOX BOUNDS xy xz yz pp pp pp')]
        for text in cases:
            with self.subTest(prefix=text[:50]),self.assertRaises(ValueError):
                parse_periodic_final_dump(text,fixture_assignment(),5)
        with self.assertRaises(ValueError):
            parse_periodic_final_dump(good,fixture_assignment(),6)


class ExecutionContractTests(unittest.TestCase):
    def fake_run(self,run_spec,*,bad_step=False,bad_box=False):
        # This is a synthetic output-schema fixture, never actual MD evidence.
        directory=run_spec.run_directory;directory.mkdir()
        log='\n'.join([_MARKERS[0],'Stopping criterion = max iterations',
                       'Iterations, force evaluations = 200 300',*_MARKERS[1:]])+'\n'
        (directory/'log.lammps').write_text(log)
        (directory/'prepared.data').write_text('LAMMPS data file, timestep = 5, units = real\n')
        (directory/'final.dump').write_text(fixture_dump(step=6 if bad_step else 5))
        (directory/'checkpoint.restart').write_bytes(b'SYNTHETIC_NATIVE_CHECKPOINT_FIXTURE_NOT_ENGINE_EVIDENCE')
        checkpoint_lines=fixture_dump(step=6 if bad_step else 5).splitlines()
        checkpoint_lines[8]+=' vx vy vz'
        checkpoint_lines[9:]=[line+' 0.001 -0.002 0.003' for line in checkpoint_lines[9:]]
        (directory/'checkpoint.dump').write_text('\n'.join(checkpoint_lines)+'\n')
        density=(1.008+1.008)/(24**3)/.602214129
        def row(step):
            return [step,300.,1.,24**3,density,1.,2.,.1,.2,.3,.4,.5,-.5,19. if bad_box else 24.,24.,24.]
        blocks=[{'columns':list(_THERMO_COLUMNS),'rows':[row(0),row(200)]},
                {'columns':list(_THERMO_COLUMNS),'rows':[row(0),row(2)]},
                {'columns':list(_THERMO_COLUMNS),'rows':[row(2),row(5)]}]
        return LammpsRunResult({'status':'completed','evidence_kind':'synthetic_process_fixture',
                                'simulation_performed':True,'thermo_blocks':blocks,
                                'engine_identity':{'fixture':'not_a_real_engine'},
                                'engine_version_banner':'SYNTHETIC_PROCESS_FIXTURE'})

    def test_executed_pilot_does_not_promote_sampling_or_minimizer_status(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)/'run'
            with patch('polychemmd.engines.periodic_alkane.build_periodic_alkane_setup',side_effect=setup_fixture),\
                 patch('polychemmd.protocols.alkane.run_lammps',side_effect=self.fake_run):
                report=run_periodic_alkane_protocol(SimpleNamespace(to_dict=fixture_assignment),
                    PeriodicAlkaneProtocolSpec(nvt_steps=2,npt_steps=3),run_directory=directory)
            self.assertEqual(report['status'],'completed_preparation_pilot')
            self.assertEqual(report['stages']['minimization']['status'],'budget_exhausted')
            self.assertEqual(report['md_integration_steps_executed'],5)
            self.assertEqual(report['equilibrium_status'],'insufficient_sampling')
            self.assertFalse(report['production_ready'])
            self.assertEqual(report['density_prediction_validation'],'not_performed')
            self.assertAlmostEqual(report['final_state']['pressure_bar'],1.01325)
            self.assertAlmostEqual(report['final_state']['density_g_cm3'],
                                   density_from_mass_volume([1.008,1.008],24**3),places=15)
            self.assertEqual(report['final_state']['density_conversion_audit']['engine_abs_difference_g_cm3'],0)
            self.assertNotEqual(report['final_state']['density_g_cm3'],report['final_state']['engine_density_g_cm3'])
            self.assertTrue((directory/'protocol_record.json').is_file())
            self.assertEqual(report['checkpoint']['step'],5)
            self.assertEqual(report['checkpoint']['integrator']['fix_id'],'polychemmd_npt')

    def test_malformed_or_out_of_regime_actual_output_is_a_failed_protocol(self):
        for kwargs in ({'bad_step':True},{'bad_box':True}):
            with self.subTest(kwargs=kwargs),tempfile.TemporaryDirectory() as temporary:
                directory=Path(temporary)/'run'
                with patch('polychemmd.engines.periodic_alkane.build_periodic_alkane_setup',side_effect=setup_fixture),\
                     patch('polychemmd.protocols.alkane.run_lammps',side_effect=lambda spec:self.fake_run(spec,**kwargs)):
                    with self.assertRaises(PeriodicProtocolError) as caught:
                        run_periodic_alkane_protocol(SimpleNamespace(to_dict=fixture_assignment),
                            PeriodicAlkaneProtocolSpec(nvt_steps=2,npt_steps=3),run_directory=directory)
                self.assertEqual(caught.exception.record['status'],'failed')
                saved=json.loads((directory/'protocol_record.json').read_text())
                self.assertEqual(saved['failure_phase'],'protocol_output_audit')

    def test_unavailable_or_stale_engine_path_creates_no_misleading_result(self):
        with tempfile.TemporaryDirectory() as temporary:
            for existing in (False,True):
                directory=Path(temporary)/('existing' if existing else 'absent')
                if existing:directory.mkdir()
                with patch('polychemmd.engines.periodic_alkane.build_periodic_alkane_setup',side_effect=setup_fixture),\
                     patch('polychemmd.protocols.alkane.run_lammps',side_effect=LammpsError('fixture unavailable or stale')):
                    with self.assertRaises(PeriodicProtocolError) as caught:
                        run_periodic_alkane_protocol(SimpleNamespace(to_dict=fixture_assignment),
                            PeriodicAlkaneProtocolSpec(),run_directory=directory)
                self.assertEqual(caught.exception.record['failure_phase'],'engine')
                self.assertFalse((directory/'protocol_record.json').exists())
                self.assertEqual(directory.exists(),existing)

    def test_final_unwrapped_bond_image_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as temporary:
            directory=Path(temporary)/'run'
            run=self.fake_run(SimpleNamespace(run_directory=directory)).to_dict()
            assignment=fixture_assignment()
            # Source identities remain correct, but the fixture image flags
            # would stretch this asserted bond across multiple box images.
            assignment['interactions']['bonds']=[{'atom_ids':['molecule:0/a','molecule:0/b']}]
            with self.assertRaisesRegex(ValueError,'unwrapped bond'):
                _audit_outputs(run,PeriodicAlkaneProtocolSpec(nvt_steps=2,npt_steps=3),assignment,directory)


if __name__=='__main__':
    unittest.main()
