import contextlib
import hashlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.cli import main

ROOT = Path(__file__).resolve().parents[1]
GEOMETRY = ROOT/'evidence'/'geometry_v03'/'conformer.json'


class ForceFieldCliTests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status=main(args)
        return status,out.getvalue(),err.getvalue()

    def test_assignment_has_exact_geometry_input_hash_and_supported_chemistry(self):
        status,text,error=self.invoke(['parameterize-alkane',str(GEOMETRY)])
        self.assertEqual(status,0,error)
        result=json.loads(text)
        self.assertEqual(result['record_kind'],'alkane_parameter_assignment')
        self.assertEqual(result['input_sha256'],hashlib.sha256(GEOMETRY.read_bytes()).hexdigest())
        payload=result['payload']
        self.assertEqual(payload['validation']['molecular_formula'],'C10H22')
        self.assertTrue(payload['validation']['parameters_assigned'])
        self.assertFalse(payload['validation']['production_ready'])
        self.assertEqual(payload['source_snapshot']['tag'],'v2025.3')

    def test_reference_evaluation_reports_force_units_and_no_equilibrium_claim(self):
        status,text,error=self.invoke(['alkane-energy',str(GEOMETRY)])
        self.assertEqual(status,0,error)
        result=json.loads(text)
        self.assertEqual(result['record_kind'],'alkane_energy_force_evaluation')
        payload=result['payload']
        self.assertEqual(len(payload['forces_kj_mol_angstrom']),32)
        self.assertEqual(payload['force_units'],'kJ/(mol angstrom)')
        self.assertEqual(payload['boundary_conditions'],'isolated_nonperiodic')
        self.assertFalse(payload['finite_difference_forces'])
        self.assertFalse(payload['trajectory_integration_performed'])
        self.assertEqual(payload['material_validation'],'not_performed')

    def test_singlepoint_plan_contains_reviewable_inputs_but_runs_nothing(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir=Path(directory)/'not-created'
            with patch('polychemmd.engines.alkane.run_lammps') as run:
                status,text,error=self.invoke(['alkane-lammps',str(GEOMETRY),'--run-dir',str(run_dir),'--plan'])
                run.assert_not_called()
            self.assertEqual(status,0,error)
            self.assertFalse(run_dir.exists())
        payload=json.loads(text)['payload']
        self.assertEqual(payload['engine_plan']['status'],'planned')
        self.assertFalse(payload['engine_plan']['simulation_performed'])
        self.assertIn('run 0',payload['input_script'])
        special=next(line.split() for line in payload['input_script'].splitlines()
                     if line.startswith('special_bonds '))
        self.assertEqual(special[:2],['special_bonds','lj/coul'])
        self.assertEqual([float(v) for v in special[2:]],[0.,0.,.5])
        self.assertIn('system.data',payload['staged_inputs'])

    def test_summary_cannot_overwrite_engine_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir=Path(directory)/'run'
            status,text,error=self.invoke(['alkane-lammps',str(GEOMETRY),'--run-dir',str(run_dir),
                                           '--plan','--output',str(run_dir/'forces.dump')])
            self.assertEqual(status,2)
            self.assertEqual(text,'')
            self.assertIn('preserve its artifacts',error)
            self.assertFalse(run_dir.exists())

    def test_all_new_commands_protect_source_before_calculation(self):
        for command in ('parameterize-alkane','alkane-energy','alkane-lammps'):
            with self.subTest(command=command),tempfile.TemporaryDirectory() as directory:
                source=Path(directory)/'geometry.json'
                source.write_bytes(GEOMETRY.read_bytes())
                args=[command,str(source),'--output',str(source)]
                if command=='alkane-lammps':
                    args+=['--run-dir',str(Path(directory)/'run'),'--plan']
                status,text,error=self.invoke(args)
                self.assertEqual(status,2)
                self.assertEqual(text,'')
                self.assertIn('must not overwrite',error)
                self.assertEqual(source.read_bytes(),GEOMETRY.read_bytes())

    def test_unsupported_chemistry_does_not_write_an_assignment(self):
        source=ROOT/'evidence'/'geometry_v03'/'pilot-pvdf-dp8.json'
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)/'assignment.json'
            status,text,error=self.invoke(['parameterize-alkane',str(source),'--output',str(target)])
            self.assertEqual(status,2)
            self.assertEqual(text,'')
            self.assertIn('only explicit C and H',error)
            self.assertFalse(target.exists())

    def test_invalid_cutoff_and_missing_engine_never_become_success(self):
        status,text,error=self.invoke(['alkane-energy',str(GEOMETRY),'--cutoff-angstrom','nan'])
        self.assertEqual(status,2)
        self.assertEqual(text,'')
        with tempfile.TemporaryDirectory() as directory:
            target=Path(directory)/'result.json'
            status,text,error=self.invoke(['alkane-lammps',str(GEOMETRY),'--run-dir',str(Path(directory)/'run'),
                                           '--executable',str(Path(directory)/'missing-lmp'),'--output',str(target)])
            self.assertEqual(status,2)
            self.assertEqual(text,'')
            self.assertIn('unavailable',error)
            self.assertFalse(target.exists())


if __name__=='__main__':
    unittest.main()
