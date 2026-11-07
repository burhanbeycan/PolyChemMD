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
from polychemmd.geometry import MolecularGeometry
from polychemmd.builders import tile_molecular_geometry
from polychemmd.forcefields import assign_compact_periodic_alkanes


ROOT = Path(__file__).resolve().parents[1]


class LargerCellCliTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        source = json.loads((ROOT/'evidence/geometry_v03/cell.json').read_text())['payload']
        cls.geometry = source['molecules'][0]['source_geometry']
        cls.recipe = json.loads((ROOT/'examples/tiled_decane.json').read_text())
        cls.cell = tile_molecular_geometry(MolecularGeometry.from_mapping(cls.geometry), **cls.recipe).to_dict()
        cls.compact = assign_compact_periodic_alkanes(cls.cell).to_dict()

    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(args)
        return status, out.getvalue(), err.getvalue()

    def test_tiling_cli_preserves_geometry_and_recipe_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry, recipe = Path(directory)/'geometry.json', Path(directory)/'recipe.json'
            geometry.write_text(json.dumps(self.geometry)); recipe.write_text(json.dumps(self.recipe))
            status, text, err = self.invoke(['tile-cell', str(geometry), '--recipe', str(recipe)])
        self.assertEqual(status, 0, err)
        result = json.loads(text)
        self.assertEqual(result['input_sha256'], hashlib.sha256(json.dumps(self.geometry).encode()).hexdigest())
        self.assertEqual(result['payload']['recipe_input_sha256'], hashlib.sha256(json.dumps(self.recipe).encode()).hexdigest())
        self.assertEqual(result['payload']['atom_count'], 2304)
        self.assertAlmostEqual(result['payload']['construction_density_g_cm3'], .6987592397, places=9)
        self.assertEqual(result['payload']['molecules'][0]['source_geometry'], self.geometry)

    def test_large_source_requires_explicit_compact_choice(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory)/'cell.json'; source.write_text(json.dumps(self.cell))
            status, text, err = self.invoke(['parameterize-cell', str(source)])
            self.assertEqual(status, 2)
            self.assertEqual(text, '')
            self.assertIn('max_atoms', err)
            status, text, err = self.invoke(['parameterize-cell', str(source), '--compact'])
        self.assertEqual(status, 0, err)
        data = json.loads(text)['payload']
        self.assertEqual(data['schema_version'], '2')
        self.assertNotIn('nonbonded_pairs', data)
        self.assertEqual(len(data['atoms']), 2304)
        self.assertLess(len(data['special_pairs']), 2304*6)

    def test_engine_only_plan_has_no_reference_or_execution_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory); source = path/'compact.json'; run_dir = path/'absent'
            source.write_text(json.dumps(self.compact))
            with patch('polychemmd.engines.periodic_alkane.run_lammps') as run:
                status, text, err = self.invoke(['periodic-singlepoint', str(source), '--run-dir', str(run_dir), '--plan'])
                run.assert_not_called()
            self.assertEqual(status, 0, err)
            self.assertFalse(run_dir.exists())
        result = json.loads(text)
        self.assertEqual(result['record_kind'], 'periodic_ewald_engine_singlepoint_plan')
        self.assertEqual(result['payload']['independent_reference_status'], 'not_performed')
        self.assertFalse(result['payload']['engine_plan']['simulation_performed'])
        self.assertIn('2304 atoms', result['payload']['staged_inputs']['system.data'])
        self.assertIn('singlepoint.data', result['payload']['engine_plan']['expected_artifacts'])

    def test_large_reference_fails_before_creating_a_numeric_result(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory)/'compact.json', Path(directory)/'result.json'
            source.write_text(json.dumps(self.compact))
            status, _, err = self.invoke(['periodic-energy', str(source), '--output', str(output)])
            self.assertEqual(status, 2)
            self.assertIn('512', err)
            self.assertFalse(output.exists())

    def test_tiling_recipe_alias_or_unknown_field_cannot_overwrite_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            geometry, recipe = Path(directory)/'geometry.json', Path(directory)/'recipe.json'
            geometry.write_text(json.dumps(self.geometry)); recipe.write_text(json.dumps(self.recipe))
            before = recipe.read_bytes()
            status, _, err = self.invoke(['tile-cell', str(geometry), '--recipe', str(recipe), '--output', str(recipe)])
            self.assertEqual(status, 2)
            self.assertIn('overwrite', err)
            self.assertEqual(recipe.read_bytes(), before)
            invalid = dict(self.recipe, phase='liquid'); recipe.write_text(json.dumps(invalid))
            output = Path(directory)/'not-created.json'
            status, _, err = self.invoke(['tile-cell', str(geometry), '--recipe', str(recipe), '--output', str(output)])
            self.assertEqual(status, 2)
            self.assertIn('exactly', err)
            self.assertFalse(output.exists())

    def test_new_modules_remain_dependency_free_on_import(self):
        code = ('import sys;sys.path.insert(0,sys.argv[1]);'
                'import polychemmd.builders,polychemmd.forcefields,polychemmd.engines;'
                'assert not any(n in sys.modules for n in ("rdkit","numpy","lammps"))')
        result = subprocess.run([sys.executable, '-S', '-c', code, str(ROOT/'src')],
                                cwd=tempfile.gettempdir(), capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == '__main__':
    unittest.main()
