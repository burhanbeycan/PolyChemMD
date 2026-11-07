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
from polychemmd.builders import tile_molecular_geometry
from polychemmd.geometry import MolecularGeometry


ROOT = Path(__file__).resolve().parents[1]
CELL = ROOT/'evidence/geometry_v03/cell.json'
PREPARATION = ROOT/'evidence/restart_v06/preparation/protocol.json'


class QualificationCliTests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(args)
        return status, out.getvalue(), err.getvalue()

    def test_static_thermo_plan_executes_nothing_and_retains_input_identity(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory)/'absent'
            with patch('polychemmd.engines.periodic_thermo.run_lammps') as run:
                status, text, err = self.invoke(['periodic-thermo', str(CELL), '--run-dir', str(target), '--plan'])
                run.assert_not_called()
            self.assertEqual(status, 0, err)
            self.assertFalse(target.exists())
        result = json.loads(text)
        self.assertEqual(result['record_kind'], 'periodic_static_thermo_plan')
        self.assertEqual(result['input_sha256'], hashlib.sha256(CELL.read_bytes()).hexdigest())
        self.assertEqual(result['payload']['status'], 'planned')
        self.assertIn('pressure NULL virial', result['payload']['input_script'])

    def test_numerical_plan_retains_protocol_hash_and_declared_run_count(self):
        spec = {'cutoffs_angstrom': [6., 8., 10.], 'reciprocal_nmax': 2,
                'refined_reciprocal_nmax': 4, 'max_engine_runs': 4}
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory); protocol = path/'spec.json'; protocol.write_text(json.dumps(spec))
            with patch('polychemmd.engines.periodic_thermo.run_lammps') as run:
                status, text, err = self.invoke(['numerical-qualify', str(CELL), '--protocol', str(protocol),
                                                '--run-dir', str(path/'absent'), '--plan'])
                run.assert_not_called()
            self.assertEqual(status, 0, err)
            self.assertFalse((path/'absent').exists())
        result = json.loads(text)['payload']
        self.assertEqual(result['status'], 'planned')
        self.assertEqual(result['expected_engine_runs'], 4)
        self.assertFalse(result['production_ready'])
        self.assertEqual(result['protocol_input_sha256'], hashlib.sha256(json.dumps(spec).encode()).hexdigest())

    def test_ordered_structure_is_one_static_cell_with_reported_chain_dimensions(self):
        source = json.loads(CELL.read_text())['payload']['molecules'][0]['source_geometry']
        cell = tile_molecular_geometry(MolecularGeometry.from_mapping(source), repeats=[2, 2, 2]).to_dict()
        from polychemmd.forcefields import assign_compact_periodic_alkanes
        compact = assign_compact_periodic_alkanes(cell).to_dict()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'assigned.json'; path.write_text(json.dumps(compact))
            status, text, err = self.invoke(['chain-structure', str(path)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)['payload']
        self.assertAlmostEqual(result['nematic_order_parameter'], 1.)
        self.assertGreater(result['mean_mass_weighted_radius_of_gyration_angstrom'], 0.)
        self.assertGreater(result['mean_end_to_end_distance_angstrom'], 0.)
        self.assertEqual(result['independent_material_cell_count'], 1)
        self.assertFalse(result['configuration_statistics_are_equilibrium_averages'])
        self.assertEqual(result['material_property_validation'], 'not_performed')

    def test_evolved_state_cli_preserves_state_hash(self):
        status, text, err = self.invoke(['chain-structure', str(CELL), '--state', str(PREPARATION)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)['payload']
        self.assertEqual(result['state_input_sha256'], hashlib.sha256(PREPARATION.read_bytes()).hexdigest())
        self.assertEqual(result['step'], 200)
        self.assertFalse(result['production_ready'])

    def test_actual_short_preparation_is_unqualified_and_thermal_forces_are_not_minimizer_residuals(self):
        with tempfile.TemporaryDirectory() as directory:
            spec = Path(directory)/'spec.json'; spec.write_text('{}')
            status, text, err = self.invoke(['preparation-review', str(CELL), '--record', str(PREPARATION),
                                            '--protocol', str(spec)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)['payload']
        self.assertEqual(result['status'], 'requires_preparation')
        self.assertFalse(result['gates']['minimizer_policy'])
        self.assertFalse(result['gates']['npt_density_sampling'])
        self.assertTrue(result['gates']['finite_final_forces_and_requested_development_bounds'])
        self.assertEqual(result['record_input_sha256'], hashlib.sha256(PREPARATION.read_bytes()).hexdigest())
        self.assertFalse(result['production_ready'])

    def test_secondary_input_aliases_and_engine_directory_outputs_fail_before_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory); spec = path/'spec.json'; spec.write_text('{}'); before = spec.read_bytes()
            status, _, err = self.invoke(['numerical-qualify', str(CELL), '--protocol', str(spec),
                '--run-dir', str(path/'run'), '--plan', '--output', str(spec)])
            self.assertEqual(status, 2)
            self.assertIn('overwrite', err)
            self.assertEqual(spec.read_bytes(), before)
            status, _, err = self.invoke(['periodic-thermo', str(CELL), '--run-dir', str(path/'run'),
                                          '--plan', '--output', str(path/'run/log.lammps')])
            self.assertEqual(status, 2)
            self.assertIn('preserve', err)
            self.assertFalse((path/'run').exists())

    def test_new_exports_import_without_optional_packages_in_any_layer_order(self):
        for order in ('forcefields,analysis,protocols,engines', 'analysis,engines,protocols,forcefields'):
            code = ('import sys,importlib;sys.path.insert(0,sys.argv[1]);'
                    '[importlib.import_module("polychemmd."+n) for n in sys.argv[2].split(",")];'
                    'assert not any(n in sys.modules for n in ("numpy","rdkit","lammps"))')
            result = subprocess.run([sys.executable, '-S', '-c', code, str(ROOT/'src'), order],
                                    cwd=tempfile.gettempdir(), capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
