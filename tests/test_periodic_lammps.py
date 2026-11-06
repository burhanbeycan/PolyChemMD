"""Export/parser and failure contracts; actual Ewald runs are separate evidence."""

import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.engines.alkane import AlkaneParityError
from polychemmd.engines.periodic_alkane import (
    _COMPONENTS, _THERMO, build_periodic_alkane_input, build_periodic_alkane_setup,
    parse_periodic_alkane_force_dump, prepare_periodic_alkane_singlepoint,
    verify_periodic_alkane_singlepoint,
)
from polychemmd.engines.lammps import LammpsRunResult
from polychemmd.forcefields.periodic import assign_periodic_alkanes


ROOT=Path(__file__).resolve().parents[1]


def dump_text(atoms, box=(24.,24.,24.), forces=None):
    forces=forces or [(0.,0.,0.)]*len(atoms)
    lines=['ITEM: TIMESTEP','0','ITEM: NUMBER OF ATOMS',str(len(atoms)),
           'ITEM: BOX BOUNDS pp pp pp',*(f'0 {v:.17g}' for v in box),
           'ITEM: ATOMS id mol x y z ix iy iz fx fy fz']
    for i,(atom,force) in enumerate(zip(atoms,forces),1):
        lines.append(f'{i} {atom["molecule_index"]+1} '
                     +' '.join(format(v,'.17g') for v in atom['wrapped_coordinates_angstrom'])+' '
                     +' '.join(str(v) for v in atom['image_flags'])+' '
                     +' '.join(format(v,'.17g') for v in force))
    return '\n'.join(lines)+'\n'


class PeriodicDumpTests(unittest.TestCase):
    def setUp(self):
        self.atoms=[{'molecule_index':0,'wrapped_coordinates_angstrom':[.1,.2,.3],'image_flags':[-1,0,2]},
                    {'molecule_index':1,'wrapped_coordinates_angstrom':[23.9,22.,1.],'image_flags':[0,-1,0]}]

    def test_preserves_order_molecule_images_and_finite_forces(self):
        forces=((1.25,-2.5,3.75),(-.125,.25,-.5))
        lines=dump_text(self.atoms,forces=forces).splitlines()
        lines[-2:]=reversed(lines[-2:])
        self.assertEqual(parse_periodic_alkane_force_dump('\n'.join(lines),self.atoms,(24,24,24)),forces)

    def test_frame_identity_boundaries_box_and_columns_are_strict(self):
        good=dump_text(self.atoms)
        variants=[good+good,good.rsplit('\n',2)[0],good.replace('TIMESTEP\n0','TIMESTEP\n1'),
                  good.replace('pp pp pp','ff pp pp'),good.replace('0 24','0 25',1),
                  good.replace('id mol x y z','id mol xs ys zs')]
        for text in variants:
            with self.subTest(last=text[-50:]),self.assertRaises(ValueError):
                parse_periodic_alkane_force_dump(text,self.atoms,(24,24,24))

    def test_wrong_molecule_image_coordinate_duplicate_and_nonfinite_rejected(self):
        good=dump_text(self.atoms)
        atom_lines=good.splitlines()
        for index,value in ((0,'2'),(1,'3'),(2,'0.11'),(5,'0'),(8,'nan'),(10,'inf')):
            with self.subTest(index=index):
                lines=atom_lines.copy();row=lines[-2].split();row[index]=value;lines[-2]=' '.join(row)
                with self.assertRaises(ValueError):
                    parse_periodic_alkane_force_dump('\n'.join(lines),self.atoms,(24,24,24))


class PeriodicExportTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cell=json.loads((ROOT/'evidence'/'geometry_v03'/'cell.json').read_text())['payload']
        cls.assigned=assign_periodic_alkanes(cell)
        cls.record=cls.assigned.to_dict()

    def test_explicit_periodic_ewald_setup_and_complete_topology(self):
        setup,data=build_periodic_alkane_setup(self.assigned)
        self.assertIn('boundary p p p',setup)
        self.assertIn('pair_style lj/cut/coul/long 10',setup)
        self.assertIn('table 0',setup)
        self.assertIn('kspace_style ewald 1e-8',setup)
        self.assertIn('kmax/ewald 16 16 16',setup)
        self.assertNotIn('run 0',setup)
        self.assertNotIn('thermo_style',setup)
        self.assertIn('128 atoms\n124 bonds\n240 angles\n324 dihedrals',data)
        atoms=data.split('Atoms # full\n\n',1)[1].split('\n\n',1)[0].splitlines()
        self.assertEqual(len(atoms),128)
        first=atoms[0].split()
        self.assertEqual(len(first),10)
        self.assertEqual([int(v) for v in first[-3:]],self.record['atoms'][0]['image_flags'])
        self.assertEqual(int(first[1]),self.record['atoms'][0]['molecule_index']+1)
        self.assertEqual([float(v) for v in first[4:7]],self.record['atoms'][0]['wrapped_coordinates_angstrom'])
        self.assertIn('0 24 xlo xhi',data)

    def test_singlepoint_composes_stable_setup_and_runs_nothing_when_planned(self):
        setup,data=build_periodic_alkane_setup(self.assigned)
        script,other=build_periodic_alkane_input(self.assigned)
        self.assertTrue(script.startswith(setup))
        self.assertEqual(data,other)
        self.assertIn('ecoul elong pe',script)
        self.assertIn('norm no format float %.17g',script)
        self.assertIn('id mol x y z ix iy iz fx fy fz',script)
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'new'
            spec=prepare_periodic_alkane_singlepoint(self.assigned,run_directory=path)
            self.assertFalse(path.exists())
            self.assertEqual(spec.input_script,script)
            self.assertEqual(dict(spec.input_files)['system.data'],data)
            self.assertEqual(spec.units_style,'real')

    def test_invalid_cutoff_and_ewald_settings_fail_before_execution(self):
        for settings in ({'cutoff_angstrom':12.},{'cutoff_angstrom':True},
                         {'cutoff_angstrom':float('inf')},{'cutoff_angstrom':3.},
                         {'ewald_alpha_per_angstrom':0.},{'reciprocal_nmax':True},
                         {'reciprocal_nmax':0},{'reciprocal_nmax':21}):
            with self.subTest(settings=settings),self.assertRaises(ValueError):
                build_periodic_alkane_setup(self.assigned,**settings)

    def test_force_mismatch_is_failed_and_persisted_without_overflow(self):
        atoms=self.record['atoms'];zero_forces=[[0.,0.,0.] for _ in atoms]
        reference={'energy_components_kj_mol':{k:0. for k in _COMPONENTS},
                   'forces_kj_mol_angstrom':zero_forces}
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'run'
            def fixture_run(spec):
                spec.run_directory.mkdir()
                (spec.run_directory/'forces.dump').write_text(dump_text(atoms,forces=[(1e200,0.,0.) for _ in atoms]))
                # Synthetic process fixture: this is not an actual LAMMPS result.
                return LammpsRunResult({'thermo_blocks':[{'columns':list(_THERMO),'rows':[[0.]*len(_THERMO)]}],
                                        'engine_version_banner':'synthetic fixture, not LAMMPS'})
            with patch('polychemmd.forcefields.periodic_reference.evaluate_periodic_alkanes',return_value=reference),\
                 patch('polychemmd.engines.periodic_alkane.run_lammps',side_effect=fixture_run):
                with self.assertRaises(AlkaneParityError) as raised:
                    verify_periodic_alkane_singlepoint(self.assigned,run_directory=path)
            failed=json.loads((path/'periodic_parity_record.json').read_text())
            self.assertEqual(failed['status'],'failed')
            self.assertEqual(raised.exception.record['status'],'failed')
            self.assertGreater(failed['max_force_component_abs_error_kj_mol_angstrom'],1e200)
            self.assertGreater(failed['rms_force_component_error_kj_mol_angstrom'],1e200)

    def test_invalid_tolerance_fails_before_reference_or_engine(self):
        with patch('polychemmd.forcefields.periodic_reference.evaluate_periodic_alkanes') as reference,\
             patch('polychemmd.engines.periodic_alkane.run_lammps') as engine:
            for tolerance in (True,0.,-1.,float('nan')):
                with self.subTest(tolerance=tolerance),self.assertRaises(ValueError):
                    verify_periodic_alkane_singlepoint(self.assigned,run_directory='not-created',energy_abs_tol_kj_mol=tolerance)
            reference.assert_not_called();engine.assert_not_called()


if __name__=='__main__':
    unittest.main()
