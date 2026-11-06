import contextlib
import hashlib
import importlib.util
import io
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.cli import main
from polychemmd.geometry import GeometryAtom, MolecularGeometry

ROOT = Path(__file__).resolve().parents[1]
HAS_RDKIT = importlib.util.find_spec('rdkit') is not None


class GeometryCliTests(unittest.TestCase):
    def invoke(self, args):
        out, err = io.StringIO(), io.StringIO()
        with contextlib.redirect_stdout(out), contextlib.redirect_stderr(err):
            status = main(args)
        return status, out.getvalue(), err.getvalue()

    @unittest.skipUnless(HAS_RDKIT, 'optional chemistry dependency absent')
    def test_conformer_pipeline_has_exact_recipe_hash_and_finite_geometry(self):
        source = ROOT/'examples'/'conformer_pe.json'
        status, text, error = self.invoke(['build-conformer', str(source)])
        self.assertEqual(status, 0, error)
        result = json.loads(text)
        self.assertEqual(result['input_sha256'], hashlib.sha256(source.read_bytes()).hexdigest())
        self.assertEqual(result['record_kind'], 'molecular_conformer')
        geometry = MolecularGeometry.from_mapping(result['payload'])
        self.assertEqual(len(geometry.atoms), 32)
        self.assertEqual(geometry.provenance['source_graph_identity']['molecular_formula'], 'C10H22')
        self.assertFalse(geometry.provenance['simulation_performed'])
        self.assertFalse(geometry.provenance['force_field_assigned'])

    @unittest.skipUnless(HAS_RDKIT, 'optional chemistry dependency absent')
    def test_actual_pack_and_export_preserve_cell_structure(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/'cell.json'
            xyz = Path(directory)/'cell.xyz'
            status, _, error = self.invoke(['pack-cell',str(ROOT/'examples'/'packed_pe.json'),'--output',str(path)])
            self.assertEqual(status, 0, error)
            data = json.loads(path.read_text())['payload']
            self.assertEqual((data['molecule_count'],data['atom_count']), (4,128))
            self.assertEqual(data['preparation']['conformer_seeds'], [73021,73022,73023,73024])
            self.assertEqual(data['equilibrium_status'], 'not_assessed')
            self.assertEqual(data['material_validation'], 'not_performed')
            self.assertGreaterEqual(data['packing']['minimum_inter_molecular_distance_angstrom'], 1.5)
            status, _, error = self.invoke(['export-xyz',str(path),'--output',str(xyz)])
            self.assertEqual(status, 0, error)
            lines = xyz.read_text().splitlines()
            self.assertEqual(int(lines[0]),128)
            self.assertEqual(len(lines),130)
            self.assertIn('Lattice=',lines[1])
            self.assertIn('pbc="T T T"',lines[1])

    def test_invalid_seeds_fail_before_embedding_and_write_nothing(self):
        for seeds in ([True], [1,1], [-1], [2**31], []):
            with self.subTest(seeds=seeds), tempfile.TemporaryDirectory() as directory:
                recipe = json.loads((ROOT/'examples'/'packed_pe.json').read_text())
                recipe['conformer_seeds'] = seeds
                source,target = Path(directory)/'input.json',Path(directory)/'output.json'
                source.write_text(json.dumps(recipe))
                with patch('polychemmd.builders.conformers.generate_conformer') as generate:
                    status, text, error = self.invoke(['pack-cell',str(source),'--output',str(target)])
                    generate.assert_not_called()
                self.assertEqual(status,2,error)
                self.assertEqual(text,'')
                self.assertFalse(target.exists())

    def test_invalid_packing_budget_and_units_fail_before_embedding(self):
        for key,value in [('max_distance_checks',True),('max_attempts_per_molecule',0),
                          ('box_lengths_angstrom',[1e-200]*3),
                          ('min_distance_angstrom',-1),('packing_seed',False)]:
            with self.subTest(key=key), tempfile.TemporaryDirectory() as directory:
                recipe=json.loads((ROOT/'examples'/'packed_pe.json').read_text())
                recipe[key]=value
                source=Path(directory)/'input.json'
                source.write_text(json.dumps(recipe))
                with patch('polychemmd.builders.conformers.generate_conformer') as generate:
                    status,_,_=self.invoke(['pack-cell',str(source)])
                    generate.assert_not_called()
                self.assertEqual(status,2)

    def test_export_requires_supported_kind_and_does_not_modify_source(self):
        with tempfile.TemporaryDirectory() as directory:
            source=Path(directory)/'input.json'
            source.write_text('{"kind":"invented_geometry"}')
            status,text,error=self.invoke(['export-xyz',str(source)])
            self.assertEqual(status,2)
            self.assertEqual(text,'')
            self.assertIn('unsupported structure',error)
            atom=GeometryAtom('analytic:He',1,'He',0,0,4.002602,'{}')
            geometry=MolecularGeometry((atom,),(),((0,0,0),),'0'*64,'{"kind":"analytic_fixture"}')
            source.write_text(json.dumps(geometry.to_dict()))
            before=source.read_bytes()
            status,text,error=self.invoke(['export-xyz',str(source),'--output',str(source)])
            self.assertEqual(status,2)
            self.assertEqual(source.read_bytes(),before)
            status,text,error=self.invoke(['export-xyz',str(source),'--wrapped'])
            self.assertEqual(status,2)
            self.assertIn('only to periodic cells',error)

    def test_geometry_commands_protect_input_alias_before_optional_calculation(self):
        for command,example in [('build-conformer','conformer_pe.json'),('pack-cell','packed_pe.json')]:
            with self.subTest(command=command), tempfile.TemporaryDirectory() as directory:
                source=Path(directory)/'input.json'
                source.write_bytes((ROOT/'examples'/example).read_bytes())
                before=source.read_bytes()
                status,text,error=self.invoke([command,str(source),'--output',str(source)])
                self.assertEqual(status,2)
                self.assertIn('must not overwrite',error)
                self.assertEqual(source.read_bytes(),before)

    def test_import_all_geometry_modules_without_optional_packages(self):
        code='import sys; sys.path.insert(0,sys.argv[1]); import polychemmd.geometry,polychemmd.builders,polychemmd.io; assert "rdkit" not in sys.modules; assert "numpy" not in sys.modules'
        result=subprocess.run([sys.executable,'-S','-c',code,str(ROOT/'src')],
                              cwd=tempfile.gettempdir(),text=True,capture_output=True)
        self.assertEqual(result.returncode,0,result.stderr)


if __name__ == '__main__':
    unittest.main()
