"""Static-thermo contracts; synthetic fixtures are not material simulations."""

from copy import deepcopy
import hashlib
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.engines.lammps import LammpsError, LammpsRunResult
from polychemmd.engines.periodic_alkane import build_periodic_alkane_input
from polychemmd.engines.periodic_thermo import (
    BAR_PER_ATM, PeriodicThermoError, _THERMO, _geometry_digest,
    parse_static_thermo_dump, plan_periodic_alkane_thermo_singlepoint,
    run_periodic_alkane_thermo_singlepoint,
)
from polychemmd.forcefields import assign_compact_periodic_alkanes, assign_periodic_alkanes
from test_periodic_assignment import packed_ethanes


def static_dump(record, *, mutation=None):
    box=record['box']['lengths_angstrom'];atoms=record['atoms']
    lines=['ITEM: TIMESTEP','0','ITEM: NUMBER OF ATOMS',str(len(atoms)),
           'ITEM: BOX BOUNDS pp pp pp',*(f'0 {value:.17g}' for value in box),
           'ITEM: ATOMS id mol q x y z ix iy iz fx fy fz vx vy vz']
    for index,atom in enumerate(atoms,1):
        fields=[str(index),str(atom['molecule_index']+1),format(atom['charge_e'],'.17g'),
                *(format(value,'.17g') for value in atom['wrapped_coordinates_angstrom']),
                *(str(value) for value in atom['image_flags']), '1','-2','3','0','0','0']
        if mutation and index==1:
            field,value=mutation;fields[field]=value
        lines.append(' '.join(fields))
    return '\n'.join(lines)+'\n'


class StaticDumpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.record=assign_compact_periodic_alkanes(packed_ethanes()).to_dict()

    def test_charge_molecule_image_positions_zero_velocities_and_force_units(self):
        lines=static_dump(self.record).splitlines()
        lines[9:]=reversed(lines[9:])
        forces=parse_static_thermo_dump('\n'.join(lines),self.record['atoms'],self.record['box']['lengths_angstrom'])
        self.assertEqual(len(forces),16)
        self.assertEqual(forces[0],[4.184,-8.368,12.552])

    def test_identity_charge_image_position_and_nonzero_velocity_fail_independently(self):
        for mutation in ((0,'2'),(1,'99'),(2,'.999'),(3,'999'),(6,'999'),
                         (9,'nan'),(10,'1e308'),(12,'1e-100'),(14,'inf')):
            with self.subTest(mutation=mutation),self.assertRaises(ValueError):
                parse_static_thermo_dump(static_dump(self.record,mutation=mutation),
                                        self.record['atoms'],self.record['box']['lengths_angstrom'])

    def test_frame_box_and_column_contract_is_strict(self):
        good=static_dump(self.record)
        for text in (good+good,good.replace('TIMESTEP\n0','TIMESTEP\n1'),
                     good.replace('pp pp pp','ff pp pp'),good.replace('0 24','0 25',1),
                     good.replace('q x y','q xs ys'),good.rsplit('\n',2)[0]):
            with self.subTest(text=text[:50]),self.assertRaises(ValueError):
                parse_static_thermo_dump(text,self.record['atoms'],self.record['box']['lengths_angstrom'])


class StaticPlanTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cell=packed_ethanes();cls.legacy=assign_periodic_alkanes(cell)
        cls.compact=assign_compact_periodic_alkanes(cell)

    def test_exact_recipe_is_no_execution_and_has_static_compute(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'absent'
            with patch('polychemmd.engines.periodic_thermo.run_lammps') as engine:
                plan=plan_periodic_alkane_thermo_singlepoint(self.compact,run_directory=target)
                engine.assert_not_called()
            self.assertFalse(target.exists())
        self.assertEqual(plan['status'],'planned')
        self.assertFalse(plan['simulation_performed'])
        self.assertIn('compute polychemmd_virial all pressure NULL virial',plan['input_script'])
        self.assertIn('thermo_modify press polychemmd_virial',plan['input_script'])
        self.assertIn('velocity all set 0.0 0.0 0.0',plan['input_script'])
        self.assertEqual(plan['engine_plan']['required_thermo_columns'],list(_THERMO))
        self.assertEqual(plan['input_files']['input.lammps_sha256'],hashlib.sha256(plan['input_script'].encode()).hexdigest())

    def test_representations_have_same_source_geometry_and_scientific_recipe(self):
        options={'run_directory':'not-created','reciprocal_nmax':2}
        legacy=plan_periodic_alkane_thermo_singlepoint(self.legacy,**options)
        compact=plan_periodic_alkane_thermo_singlepoint(self.compact,**options)
        self.assertEqual(legacy['geometry_canonical_sha256'],compact['geometry_canonical_sha256'])
        self.assertNotEqual(legacy['assignment_canonical_sha256'],compact['assignment_canonical_sha256'])
        self.assertEqual(legacy['input_script'],compact['input_script'])
        self.assertEqual(legacy['staged_inputs'],compact['staged_inputs'])
        self.assertEqual(build_periodic_alkane_input(self.legacy),build_periodic_alkane_input(self.compact))

    def test_geometry_fingerprint_changes_for_positions_images_charge_and_box(self):
        original=self.compact.to_dict();digest=_geometry_digest(original)
        variants=[]
        for field,change in (('coordinates_angstrom',[.1,.2,.3]),('image_flags',[99,0,0]),('charge_e',.123)):
            data=deepcopy(original);data['atoms'][0][field]=change;variants.append(data)
        data=deepcopy(original);data['box']['lengths_angstrom'][0]+=1;variants.append(data)
        for data in variants:self.assertNotEqual(_geometry_digest(data),digest)

    def test_invalid_resource_or_numerical_setting_cannot_launch(self):
        bad=[{'timeout_seconds':value} for value in (True,0,-1,float('nan'),120.001)]
        bad.extend(({'cutoff_angstrom':12.},{'reciprocal_nmax':True},{'ewald_alpha_per_angstrom':float('inf')}))
        with tempfile.TemporaryDirectory() as temp,patch('polychemmd.engines.periodic_thermo.run_lammps') as engine:
            for options in bad:
                target=Path(temp)/'absent'
                with self.subTest(options=options),self.assertRaises(ValueError):
                    run_periodic_alkane_thermo_singlepoint(self.compact,run_directory=target,**options)
                self.assertFalse(target.exists())
            engine.assert_not_called()


class StaticExecutionTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.assigned=assign_compact_periodic_alkanes(packed_ethanes())
        cls.record=cls.assigned.to_dict()

    def fake_run(self,spec,mutation=None):
        # These deterministic numbers test an output contract; no engine ran.
        spec.run_directory.mkdir()
        dump_mutations={'charge':(2,'.123'),'velocity':(12,'1e-100'),'image':(6,'999')}
        (spec.run_directory/'forces.dump').write_text(static_dump(self.record,mutation=dump_mutations.get(mutation)))
        (spec.run_directory/'singlepoint.data').write_text('Synthetic fixture; units = '+('lj' if mutation=='units' else 'real')+'\n')
        volume=math.prod(self.record['box']['lengths_angstrom'])
        density=math.fsum(atom['mass_g_mol'] for atom in self.record['atoms'])/volume/.602214129
        values=[0.,0.,0.,2.,1.,2.,3.,4.,5.,6.,volume,density,1.,2.,3.,4.,5.,6.,21.]
        edits={'pressure':(3,2.5),'temperature':(1,1e-30),'kinetic':(2,1e-30),
               'volume':(10,volume+1.),'density':(11,density*1.01),'sum':(18,20.),'overflow':(12,1e308)}
        if mutation in edits:index,value=edits[mutation];values[index]=value
        return LammpsRunResult({'simulation_performed':True,
            'thermo_blocks':[{'columns':list(_THERMO),'rows':[values]}],
            'engine_version_banner':'synthetic contract fixture; not an actual engine'})

    def test_native_units_tensor_trace_and_exact_density_are_explicit(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            with patch('polychemmd.engines.periodic_thermo.run_lammps',side_effect=self.fake_run):
                result=run_periodic_alkane_thermo_singlepoint(self.assigned,run_directory=target)
            self.assertEqual(json.loads((target/'periodic_thermo_record.json').read_text()),result)
        self.assertEqual(result['status'],'completed_static_thermo_singlepoint')
        self.assertEqual(result['pressure_bar'],2*BAR_PER_ATM)
        self.assertEqual(result['pressure_tensor_bar'],[[BAR_PER_ATM,4*BAR_PER_ATM,5*BAR_PER_ATM],
            [4*BAR_PER_ATM,2*BAR_PER_ATM,6*BAR_PER_ATM],[5*BAR_PER_ATM,6*BAR_PER_ATM,3*BAR_PER_ATM]])
        # Four carbons and twelve hydrogens: 4*12.011 + 12*1.008 = 60.140 g/mol.
        self.assertAlmostEqual(result['density_g_cm3'],60.140/(.602214076*15600.),places=14)
        self.assertNotEqual(result['density_g_cm3'],result['engine_density_g_cm3'])
        self.assertEqual(result['energy_components_kj_mol']['total'],21*4.184)
        self.assertEqual(result['initial_velocities'],'zeroed_and_output_verified')
        self.assertEqual(result['pressure_interpretation'],'static_virial_not_equilibrium_average')
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['independent_reference_status'],'not_performed')
        self.assertEqual(result['md_integration_steps'],0)

    def test_failed_pressure_state_mass_or_unit_audits_retain_failed_records(self):
        for mutation in ('pressure','temperature','kinetic','volume','density','sum','overflow','units','charge','velocity','image'):
            with self.subTest(mutation=mutation),tempfile.TemporaryDirectory() as temp:
                target=Path(temp)/'run'
                with patch('polychemmd.engines.periodic_thermo.run_lammps',side_effect=lambda spec:self.fake_run(spec,mutation)):
                    with self.assertRaises(PeriodicThermoError) as raised:
                        run_periodic_alkane_thermo_singlepoint(self.assigned,run_directory=target)
                saved=json.loads((target/'periodic_thermo_record.json').read_text())
                self.assertEqual(saved,raised.exception.record)
                self.assertEqual(saved['failure_phase'],'static_thermo_output_audit')
                self.assertEqual(saved['status'],'failed')
                self.assertNotIn('pressure_bar',saved)

    def test_unavailable_or_stale_directory_creates_no_observables(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'absent'
            with patch('polychemmd.engines.periodic_thermo.run_lammps',side_effect=LammpsError('unavailable synthetic fixture')):
                with self.assertRaises(PeriodicThermoError) as raised:
                    run_periodic_alkane_thermo_singlepoint(self.assigned,run_directory=target)
            self.assertFalse(target.exists())
            self.assertFalse(raised.exception.record['simulation_performed'])
            self.assertNotIn('pressure_bar',raised.exception.record)

    def test_plan_script_and_staged_input_are_the_actual_run_request(self):
        with tempfile.TemporaryDirectory() as temp:
            target=Path(temp)/'run'
            plan=plan_periodic_alkane_thermo_singlepoint(self.assigned,run_directory=target,reciprocal_nmax=2)
            with patch('polychemmd.engines.periodic_thermo.run_lammps',side_effect=self.fake_run) as engine:
                result=run_periodic_alkane_thermo_singlepoint(self.assigned,run_directory=target,reciprocal_nmax=2)
            actual=engine.call_args.args[0]
        self.assertEqual(plan['input_script'],actual.input_script)
        self.assertEqual(plan['staged_inputs'],dict(actual.input_files))
        self.assertEqual(plan['input_files'],result['input_files'])


if __name__=='__main__':
    unittest.main()
