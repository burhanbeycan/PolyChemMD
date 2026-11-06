import contextlib
import hashlib
import io
import json
import os
from pathlib import Path
import tempfile
import unittest

from polychemmd.cli import main

ROOT = Path(__file__).resolve().parents[1]


class CliTests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(args)
        return status, out.getvalue(), err.getvalue()

    def test_validation_exposes_limits_and_exact_input_hash(self):
        path = ROOT / 'examples' / 'pe_record.json'
        status, text, err = self.invoke(['validate', str(path)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)
        self.assertEqual(result['input_sha256'], hashlib.sha256(path.read_bytes()).hexdigest())
        self.assertEqual(result['payload']['molecular_graph_validation'], 'not_performed')
        self.assertFalse(result['payload']['simulation_performed'])

    def test_nonfinite_json_has_no_result(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'bad.json'
            path.write_text('{"x":NaN}')
            status, out, err = self.invoke(['validate', str(path)])
        self.assertEqual(status, 2)
        self.assertEqual(out, '')
        self.assertIn('non-finite', err)

    def test_duplicate_json_fields_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'bad.json'
            path.write_text('{"terms":[],"terms":[]}')
            status, out, err = self.invoke(['torsion-describe', str(path)])
        self.assertEqual(status, 2)
        self.assertEqual(out, '')
        self.assertIn('duplicate JSON field', err)

    def test_input_and_output_aliases_are_rejected_without_mutation(self):
        for command, example in [('validate', 'pe_record.json'), ('torsion-describe', 'torsion_demo.json')]:
            with self.subTest(command=command), tempfile.TemporaryDirectory() as directory:
                path = Path(directory) / 'input.json'
                raw = (ROOT/'examples'/example).read_bytes()
                path.write_bytes(raw)
                hardlink = Path(directory)/'hardlink.json'
                os.link(path, hardlink)
                for output in (path, path.parent/'.'/path.name, hardlink):
                    status, out, err = self.invoke([command, str(path), '--output', str(output)])
                    self.assertEqual(status, 2)
                    self.assertIn('must not overwrite', err)
                    self.assertEqual(path.read_bytes(), raw)
                    self.assertEqual(out, '')

    def test_large_integer_is_a_validation_error_not_a_crash(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'input.json'
            data = json.loads((ROOT/'examples'/'pe_record.json').read_text())
            data['state']['temperature_k'] = 10**400
            path.write_text(json.dumps(data))
            status, out, err = self.invoke(['validate', str(path)])
        self.assertEqual(status, 2)
        self.assertEqual(out, '')
        self.assertIn('finite number', err)

    def test_catalogue_cannot_claim_planned_properties_implemented(self):
        status, out, err = self.invoke(['catalogue'])
        self.assertEqual(status, 0, err)
        families = json.loads(out)['families']
        self.assertGreater(len(families), 10)
        self.assertTrue(all(not row['implemented'] and row['status'] == 'planned' for row in families))

    def test_torsion_output_is_labeled_analytic(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'result.json'
            status, out, err = self.invoke(['torsion-describe', str(ROOT/'examples'/'torsion_demo.json'), '--n-grid', '360', '--output', str(path)])
            self.assertEqual(status, 0, err)
            self.assertEqual(out, '')
            result = json.loads(path.read_text())
        self.assertEqual(result['record_kind'], 'analytic_torsion_descriptor')
        self.assertFalse(result['payload']['material_validated'])
        self.assertFalse(result['payload']['simulation_performed'])

    def test_invalid_torsion_does_not_write_new_result(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'result.json'
            status, out, err = self.invoke(['torsion-describe', str(ROOT/'examples'/'torsion_demo.json'), '--temperature-k', '-1', '--output', str(path)])
            self.assertEqual(status, 2)
            self.assertFalse(path.exists())
            self.assertEqual(out, '')

    def test_doctor_is_only_a_presence_check(self):
        status, out, err = self.invoke(['doctor'])
        self.assertEqual(status, 0, err)
        data = json.loads(out)
        self.assertEqual(data['md_adapter_status'], 'implemented_process_adapter')
        self.assertEqual(data['polymer_md_protocol_status'], 'bounded_periodic_preparation_only')
        self.assertEqual(data['polymer_equilibrium_validation'], 'not_implemented')
        self.assertIn('Presence check only', data['scope'])
