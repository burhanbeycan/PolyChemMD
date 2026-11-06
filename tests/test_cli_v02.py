import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest

from polychemmd.cli import main

ROOT = Path(__file__).resolve().parents[1]


class AdvancedCliTests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(args)
        return status, out.getvalue(), err.getvalue()

    @unittest.skipUnless(importlib.util.find_spec('rdkit'), 'optional chemistry dependency absent')
    def test_build_linear_is_graph_evidence(self):
        status, text, err = self.invoke(['build-linear', str(ROOT/'examples'/'linear_pe.json')])
        self.assertEqual(status, 0, err)
        result = json.loads(text)
        self.assertEqual(result['record_kind'], 'linear_molecular_graph')
        self.assertEqual(result['payload']['evidence_level'], 'graph_validated')
        self.assertFalse(result['payload']['validation']['geometry_generated'])

    def test_torsion_conversion_is_not_parameter_validation(self):
        status, text, err = self.invoke(['torsion-import', str(ROOT/'examples'/'torsion_lammps_fourier.json'), '--source', 'lammps-fourier'])
        self.assertEqual(status, 0, err)
        result = json.loads(text)
        self.assertEqual(result['payload']['parameter_physical_validation'], 'not_performed')
        self.assertFalse(result['payload']['simulation_performed'])

    def test_coverage_demo_is_not_production_ready(self):
        status, text, err = self.invoke(['parameter-check', str(ROOT/'examples'/'parameter_contract_demo.json'), '--require-complete'])
        self.assertEqual(status, 0, err)
        result = json.loads(text)
        self.assertFalse(result['payload']['production_ready'])

    def test_engine_plan_runs_nothing_and_creates_no_run_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'not-created'
            status, text, err = self.invoke(['lammps-run', str(ROOT/'examples'/'lammps_smoke.in'), '--run-dir', str(path), '--plan'])
            self.assertEqual(status, 0, err)
            self.assertFalse(path.exists())
        self.assertFalse(json.loads(text)['payload']['simulation_performed'])

    def test_engine_summary_cannot_replace_thermo_log(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'run'
            status, text, err = self.invoke(['lammps-run', str(ROOT/'examples'/'lammps_smoke.in'), '--run-dir', str(path), '--plan', '--output', str(path/'log.lammps')])
            self.assertEqual(status, 2)
            self.assertFalse(path.exists())
        self.assertIn('preserve its artifacts', err)

    def test_density_requires_explicit_physical_units(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'density.json'
            data = {'units':'lj', 'densities':[.8442]*64, 'source_reference':'analytic check', 'source_kind':'analytic_demo'}
            path.write_text(json.dumps(data))
            status, text, err = self.invoke(['density-analyze', str(path), '--block-size', '4'])
        self.assertEqual(status, 2)
        self.assertIn('reduced LJ units', err)

    def test_density_diagnostics_do_not_certify_equilibrium(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'density.json'
            data = {'units':'g/cm^3', 'densities':[.9,.91,.89,.92]*16, 'source_reference':'analytic numerical check, no measured polymer', 'source_kind':'analytic_demo'}
            path.write_text(json.dumps(data))
            status, text, err = self.invoke(['density-analyze', str(path), '--block-size', '2'])
        self.assertEqual(status, 0, err)
        result = json.loads(text)['payload']
        self.assertEqual(result['equilibrium_status'], 'not_assessed')
        self.assertEqual(result['polymer_density_validation'], 'not_performed')
