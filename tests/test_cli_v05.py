import contextlib
import hashlib
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.cli import main
from test_periodic_assignment import packed_ethanes


ROOT=Path(__file__).resolve().parents[1]
CELL=ROOT/'evidence'/'geometry_v03'/'cell.json'
SPEC=ROOT/'examples'/'periodic_protocol.json'


class PeriodicCliTests(unittest.TestCase):
    def invoke(self,args):
        out,err=io.StringIO(),io.StringIO()
        with contextlib.redirect_stdout(out),contextlib.redirect_stderr(err):
            status=main(args)
        return status,out.getvalue(),err.getvalue()

    def test_periodic_assignment_has_exact_input_hash_and_inter_molecular_pairs(self):
        status,text,error=self.invoke(['parameterize-cell',str(CELL)])
        self.assertEqual(status,0,error)
        result=json.loads(text)
        self.assertEqual(result['record_kind'],'periodic_alkane_parameter_assignment')
        self.assertEqual(result['input_sha256'],hashlib.sha256(CELL.read_bytes()).hexdigest())
        data=result['payload']
        self.assertEqual(len(data['atoms']),128)
        self.assertEqual(len(data['nonbonded_pairs']),128*127//2)
        self.assertEqual(sum(not p['same_molecule'] for p in data['nonbonded_pairs']),6144)
        self.assertEqual(len({a['cell_atom_id'] for a in data['atoms']}),128)

    def test_periodic_reference_cli_is_finite_sum_not_equilibrium_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'cell.json'
            source.write_text(json.dumps(packed_ethanes()))
            status,text,error=self.invoke(['periodic-energy',str(source),'--reciprocal-nmax','2'])
        self.assertEqual(status,0,error)
        data=json.loads(text)['payload']
        self.assertEqual(len(data['forces_kj_mol_angstrom']),16)
        self.assertEqual(data['reciprocal_nmax'],2)
        self.assertIn('periodic',data['boundary_conditions'])
        self.assertEqual(data['md_integration_steps'],0)
        self.assertEqual(data['equilibrium_status'],'not_assessed')
        self.assertEqual(data['material_property_validation'],'not_performed')

    def test_periodic_singlepoint_plan_runs_nothing_and_preserves_images(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir=Path(directory)/'not-created'
            with patch('polychemmd.engines.periodic_alkane.run_lammps') as run:
                status,text,error=self.invoke(['periodic-lammps',str(CELL),'--run-dir',str(run_dir),'--plan'])
                run.assert_not_called()
            self.assertEqual(status,0,error)
            self.assertFalse(run_dir.exists())
        payload=json.loads(text)['payload']
        self.assertEqual(payload['engine_plan']['status'],'planned')
        self.assertFalse(payload['engine_plan']['simulation_performed'])
        self.assertIn('boundary p p p',payload['input_script'])
        self.assertIn('kspace_style ewald',payload['input_script'])
        self.assertIn('system.data',payload['staged_inputs'])

    def test_protocol_plan_keeps_spec_hash_and_creates_no_engine_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            run_dir=Path(directory)/'not-created'
            with patch('polychemmd.protocols.alkane.run_lammps') as run:
                status,text,error=self.invoke(['periodic-prepare',str(CELL),'--protocol',str(SPEC),
                                               '--run-dir',str(run_dir),'--plan'])
                run.assert_not_called()
            self.assertEqual(status,0,error)
            self.assertFalse(run_dir.exists())
        result=json.loads(text)
        self.assertEqual(result['record_kind'],'periodic_preparation_plan')
        payload=result['payload']
        self.assertEqual(payload['status'],'planned')
        self.assertFalse(payload['simulation_performed'])
        self.assertEqual(payload['md_integration_steps_requested'],200)
        self.assertFalse(payload['production_ready'])
        self.assertEqual(payload['protocol_input_sha256'],hashlib.sha256(SPEC.read_bytes()).hexdigest())

    def test_summary_and_protocol_input_aliases_are_rejected_before_execution(self):
        for command in ('periodic-lammps','periodic-prepare'):
            with self.subTest(command=command),tempfile.TemporaryDirectory() as directory:
                run_dir=Path(directory)/'run'
                args=[command,str(CELL),'--run-dir',str(run_dir),'--plan','--output',str(run_dir/'log.lammps')]
                if command=='periodic-prepare':
                    args+=['--protocol',str(SPEC)]
                status,text,error=self.invoke(args)
                self.assertEqual(status,2)
                self.assertEqual(text,'')
                self.assertIn('preserve its artifacts',error)
                self.assertFalse(run_dir.exists())
        with tempfile.TemporaryDirectory() as directory:
            spec=Path(directory)/'protocol.json'
            spec.write_bytes(SPEC.read_bytes());before=spec.read_bytes()
            status,_,error=self.invoke(['periodic-prepare',str(CELL),'--protocol',str(spec),
                                        '--run-dir',str(Path(directory)/'run'),'--plan','--output',str(spec)])
            self.assertEqual(status,2)
            self.assertIn('must not overwrite',error)
            self.assertEqual(spec.read_bytes(),before)

    def test_nonperiodic_input_cannot_generate_an_assigned_result(self):
        with tempfile.TemporaryDirectory() as directory:
            data=packed_ethanes();data['box']['periodic'][0]=False
            source,target=Path(directory)/'input.json',Path(directory)/'result.json'
            source.write_text(json.dumps(data))
            status,text,error=self.invoke(['parameterize-cell',str(source),'--output',str(target)])
            self.assertEqual(status,2)
            self.assertEqual(text,'')
            self.assertFalse(target.exists())

    def test_all_periodic_modules_import_without_optional_packages(self):
        code='import sys;sys.path.insert(0,sys.argv[1]);import polychemmd.forcefields,polychemmd.engines,polychemmd.protocols;assert "rdkit" not in sys.modules and "numpy" not in sys.modules'
        result=subprocess.run([sys.executable,'-S','-c',code,str(ROOT/'src')],
                              cwd=tempfile.gettempdir(),capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__=='__main__':
    unittest.main()
