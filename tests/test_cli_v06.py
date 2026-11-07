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


ROOT = Path(__file__).resolve().parents[1]


class SamplingCampaignCliTests(unittest.TestCase):
    def invoke(self, arguments):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(arguments)
        return status, out.getvalue(), err.getvalue()

    def test_campaign_example_exposes_capacity_block_and_executes_nothing(self):
        source = ROOT/'examples/density_campaign.json'
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)/'plan.json'
            with patch('polychemmd.engines.lammps.run_lammps') as run:
                status, _, err = self.invoke(['density-plan', str(source), '--output', str(output)])
                run.assert_not_called()
            self.assertEqual(status, 0, err)
            result = json.loads(output.read_text())
        self.assertEqual(result['input_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(result['record_kind'], 'density_campaign_plan')
        self.assertFalse(result['payload']['simulation_performed'])
        self.assertFalse(result['payload']['execution_supported'])
        self.assertIn('incompatible', json.dumps(result['payload']))

    def test_real_pilot_input_is_insufficient_and_spec_hash_is_preserved(self):
        pilot = json.loads((ROOT/'evidence/periodic_v05/preparation_pilot/protocol.json').read_text())
        rows = pilot['stages']['npt']['thermo_rows']
        rows = [{key: row[key] for key in ('step', 'temperature_k', 'pressure_bar', 'density_g_cm3')}
                for row in rows]
        data = {'units': {'density': 'g/cm^3', 'temperature': 'K', 'pressure': 'bar'},
                'rows': rows, 'source_reference': 'periodic_v05/preparation_pilot/protocol.json',
                'source_kind': 'simulation'}
        spec = ROOT/'examples/density_sampling.json'
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'series.json'
            source.write_text(json.dumps(data))
            status, text, err = self.invoke(['density-sampling', str(source), '--protocol', str(spec)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)['payload']
        self.assertEqual(result['status'], 'insufficient_sampling')
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['protocol_input_sha256'], hashlib.sha256(spec.read_bytes()).hexdigest())

    def test_sampling_rejects_reduced_units_and_output_aliases(self):
        with tempfile.TemporaryDirectory() as directory:
            source, spec = Path(directory)/'series.json', ROOT/'examples/density_sampling.json'
            source.write_text(json.dumps({'units': 'LJ', 'rows': [], 'source_kind': 'simulation',
                                          'source_reference': 'reduced-unit example'}))
            status, text, err = self.invoke(['density-sampling', str(source), '--protocol', str(spec)])
            self.assertEqual(status, 2)
            self.assertEqual(text, '')
            self.assertIn('explicit', err)
            before = source.read_bytes()
            status, _, err = self.invoke(['density-sampling', str(source), '--protocol', str(spec),
                                          '--output', str(source)])
            self.assertEqual(status, 2)
            self.assertIn('overwrite', err)
            self.assertEqual(source.read_bytes(), before)

    def test_resume_plan_boundary_preserves_both_supplementary_input_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            checkpoint = path/'manifest.json'
            checkpoint.write_text('{"checkpoint":{"example":"mocked boundary test"}}')
            protocol = ROOT/'examples/periodic_continuation.json'
            with patch('polychemmd.protocols.restart.plan_periodic_alkane_continuation',
                       return_value={'status': 'planned', 'simulation_performed': False}) as plan:
                status, text, err = self.invoke(['periodic-resume', str(ROOT/'evidence/geometry_v03/cell.json'),
                    '--checkpoint', str(checkpoint), '--checkpoint-dir', str(path/'absent-checkpoint-dir'),
                    '--protocol', str(protocol), '--run-dir', str(path/'absent-run-dir'), '--plan'])
            self.assertEqual(status, 0, err)
            self.assertEqual(plan.call_args.args[1], {'example': 'mocked boundary test'})
            result = json.loads(text)['payload']
            self.assertEqual(result['checkpoint_input_sha256'], hashlib.sha256(checkpoint.read_bytes()).hexdigest())
            self.assertEqual(result['protocol_input_sha256'], hashlib.sha256(protocol.read_bytes()).hexdigest())
            self.assertFalse((path/'absent-run-dir').exists())

    def test_resume_cannot_overwrite_checkpoint_artifacts(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)
            checkpoint_dir = path/'producer'
            status, _, err = self.invoke(['periodic-resume', str(ROOT/'evidence/geometry_v03/cell.json'),
                '--checkpoint', str(path/'unused.json'), '--checkpoint-dir', str(checkpoint_dir),
                '--protocol', str(ROOT/'examples/periodic_continuation.json'),
                '--run-dir', str(path/'run'), '--plan', '--output', str(checkpoint_dir/'checkpoint.restart')])
            self.assertEqual(status, 2)
            self.assertIn('preserve', err)
            self.assertFalse(checkpoint_dir.exists())

    def test_new_modules_import_without_optional_dependencies(self):
        code = ('import sys;sys.path.insert(0,sys.argv[1]);'
                'import polychemmd.analysis,polychemmd.protocols;'
                'assert "rdkit" not in sys.modules and "numpy" not in sys.modules and "lammps" not in sys.modules')
        result = subprocess.run([sys.executable, '-S', '-c', code, str(ROOT/'src')],
                                cwd=tempfile.gettempdir(), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
