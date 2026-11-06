import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.engines.alkane import (
    AlkaneParityError, _opls_coefficients, _opls_backend_diagnostics, build_alkane_input,
    parse_alkane_force_dump, prepare_alkane_singlepoint, verify_alkane_singlepoint,
)
from polychemmd.engines.lammps import LammpsError, LammpsRunResult
from polychemmd.forcefields.alkanes import assign_alkane
from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry


def analytic_ethane():
    """Specified tetrahedral ethane coordinates, not a computed/observed result."""
    atoms = [GeometryAtom('C1',1,'C',0,0,12.011,'{}'),
             GeometryAtom('C2',2,'C',0,0,12.011,'{}')]
    coordinates = [(-.7645,0,0),(.7645,0,0)]
    bonds = [GeometryBond(('C1','C2'),1)]
    radial = 1.09*math.sqrt(8)/3
    for carbon, side, phase in (('C1',-1,0.),('C2',1,math.pi/3)):
        center = -.7645 if side < 0 else .7645
        for j in range(3):
            identity=f'{carbon}-H{j}'
            atoms.append(GeometryAtom(identity,len(atoms)+1,'H',0,0,1.008,'{}'))
            angle=phase+j*2*math.pi/3
            coordinates.append((center+side*1.09/3,radial*math.cos(angle),radial*math.sin(angle)))
            bonds.append(GeometryBond((carbon,identity),1))
    return MolecularGeometry(tuple(atoms),tuple(bonds),tuple(coordinates),'0'*64,
                             '{"source_kind":"analytic_fixture"}')


def dump_text(coordinates, forces=None):
    if forces is None:
        forces=[(0.,0.,0.)]*len(coordinates)
    lines=['ITEM: TIMESTEP','0','ITEM: NUMBER OF ATOMS',str(len(coordinates)),
           'ITEM: BOX BOUNDS ff ff ff','-20 20','-20 20','-20 20',
           'ITEM: ATOMS id x y z fx fy fz']
    for i,(point,force) in enumerate(zip(coordinates,forces),1):
        lines.append(str(i)+' '+' '.join(format(v,'.17g') for v in (*point,*force)))
    return '\n'.join(lines)+'\n'


class ForceDumpTests(unittest.TestCase):
    def test_exact_atom_ids_coordinates_and_force_order(self):
        coordinates=((0.,1.,2.),(3.,4.,5.))
        forces=((.25,-.5,.75),(1.5,2.5,-3.5))
        text=dump_text(coordinates,forces)
        lines=text.splitlines()
        lines[-2:]=reversed(lines[-2:])
        self.assertEqual(parse_alkane_force_dump('\n'.join(lines),coordinates),forces)

    def test_truncation_extra_frame_duplicate_id_and_changed_coordinates_fail(self):
        coordinates=((0,1,2),(3,4,5))
        good=dump_text(coordinates)
        variants=[good.rsplit('\n',2)[0],good+good,
                  good.replace('2 3 4 5 0 0 0','1 3 4 5 0 0 0'),
                  good.replace('2 3 4 5 0 0 0','2 3.001 4 5 0 0 0'),
                  good.replace('ITEM: TIMESTEP\n0','ITEM: TIMESTEP\n1'),
                  good.replace('ff ff ff','pp pp pp')]
        for text in variants:
            with self.subTest(text=text[-50:]),self.assertRaises(ValueError):
                parse_alkane_force_dump(text,coordinates)

    def test_nonfinite_and_malformed_dump_rows_fail(self):
        coordinates=((0,1,2),)
        good=dump_text(coordinates,((.25,-.5,.75),))
        for text in [good.replace('.25','nan'),good.replace('1 0 1 2','0 0 1 2'),
                     good.replace('0.75','inf'),good.replace('id x y z fx fy fz','id xs ys zs fx fy fz'),
                     good.replace('-20 20','20 -20',1),good.replace('1 0 1 2','1 0 1')]:
            with self.subTest(text=text[-70:]),self.assertRaises(ValueError):
                parse_alkane_force_dump(text,coordinates)


class BackendGeometryDiagnosticTests(unittest.TestCase):
    def screen(self, phi):
        # Specified four-point geometry, not a complete alkane assignment.
        atoms=tuple(GeometryAtom(f'C{i}',i,'C',0,0,12.011,'{}') for i in range(1,5))
        points=((0.,1.,0.),(0.,0.,0.),(1.,0.,0.),(1.,math.cos(phi),math.sin(phi)))
        bonds=tuple(GeometryBond((f'C{i}',f'C{i+1}'),1) for i in range(1,4))
        geometry=MolecularGeometry(atoms,bonds,points,'0'*64,'{"kind":"analytic_geometric_fixture"}')
        assignment={'interactions':{'dihedrals':[
            {'id':'d1','type_id':'torsion_CCCC','atom_ids':['C1','C2','C3','C4']}]}}
        return _opls_backend_diagnostics(assignment,geometry,'test fixture, not an engine')

    def test_trans_and_near_trans_below_pinned_threshold_are_flagged(self):
        for phi in (math.pi, math.pi-5e-6, -math.pi+5e-6):
            with self.subTest(phi=phi):
                report=self.screen(phi)
                self.assertEqual(len(report['near_planar_dihedrals']),1)
                flagged=report['near_planar_dihedrals'][0]
                self.assertEqual(flagged['interaction_id'],'d1')
                self.assertEqual(flagged['atom_ids'],['C1','C2','C3','C4'])
                self.assertAlmostEqual(flagged['abs_sin_phi'],abs(math.sin(phi)),places=14)
                self.assertFalse(report['reference_or_potential_modified'])
                self.assertFalse(report['parity_tolerances_modified'])

    def test_well_angled_or_above_threshold_torsions_are_not_flagged(self):
        for phi in (math.pi/3, math.pi/2, math.pi-2e-5):
            with self.subTest(phi=phi):
                report=self.screen(phi)
                self.assertEqual(report['near_planar_dihedrals'],[])
                self.assertEqual(report['observed_source_sine_threshold'],1e-5)
                self.assertIn('not guaranteed',report['threshold_scope'])


class ExportTests(unittest.TestCase):
    def setUp(self):
        self.assigned=assign_alkane(analytic_ethane())

    def test_real_unit_conversion_and_complete_ethane_topology(self):
        script,data=build_alkane_input(self.assigned,cutoff_angstrom=11.5)
        self.assertIn('units real\n',script)
        self.assertIn('boundary f f f',script)
        self.assertIn('special_bonds lj/coul 0.0 0.0 0.5',script)
        self.assertIn('mix geometric shift no tail no',script)
        self.assertIn('thermo_modify norm no format float %.17g',script)
        self.assertIn('run 0\n',script)
        self.assertIn('8 atoms\n7 bonds\n12 angles\n9 dihedrals',data)
        sections={name:data.split(name+'\n\n',1)[1].split('\n\n',1)[0]
                  for name in ('Bond Coeffs # harmonic','Angle Coeffs # harmonic','Dihedral Coeffs # opls')}
        bonds={line.split('#')[1].strip():[float(v) for v in line.split('#')[0].split()[1:]]
               for line in sections['Bond Coeffs # harmonic'].splitlines()}
        self.assertAlmostEqual(bonds['bond_CC'][0],268.)
        self.assertAlmostEqual(bonds['bond_CC'][1],1.529)
        self.assertAlmostEqual(bonds['bond_CH'][0],340.)
        self.assertAlmostEqual(bonds['bond_CH'][1],1.09)
        torsions={line.split('#')[1].strip():[float(v) for v in line.split('#')[0].split()[1:]]
                  for line in sections['Dihedral Coeffs # opls'].splitlines()}
        for actual,expected in zip(torsions['torsion_CCCC'],(1.3,-.05,.2,0.)):
            self.assertAlmostEqual(actual,expected)
        self.assertAlmostEqual(torsions['torsion_HCCH'][2],.3)

    def test_arbitrary_phase_and_unrepresentable_periodicity_fail(self):
        for term in ({'k_kj_mol':1.,'periodicity':2,'phase_rad':0.},
                     {'k_kj_mol':1.,'periodicity':5,'phase_rad':0.}):
            with self.subTest(term=term),self.assertRaises(ValueError):
                _opls_coefficients({'terms':[term]})
        self.assertEqual(_opls_coefficients({'terms':[{'k_kj_mol':2.092,'periodicity':2,'phase_rad':math.pi}]}),
                         (0.,1.,0.,0.))

    def test_plan_stages_new_data_without_creating_run_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory)/'not-created'
            spec=prepare_alkane_singlepoint(self.assigned,run_directory=run)
            self.assertFalse(run.exists())
            self.assertEqual(spec.units_style,'real')
            self.assertEqual(spec.expected_artifacts,('log.lammps','forces.dump'))
            self.assertEqual(spec.input_files[0][0],'system.data')

    def test_bad_cutoff_and_assignment_are_explicit_errors(self):
        for cutoff in (True,0.,-1.,float('nan'),float('inf'),10**400):
            with self.subTest(cutoff_type=type(cutoff)),self.assertRaises(ValueError):
                build_alkane_input(self.assigned,cutoff_angstrom=cutoff)
        with self.assertRaises(ValueError):
            build_alkane_input(self.assigned.to_dict())


class MockParityPipelineTests(unittest.TestCase):
    """Synthetic engine/reference values test gates, not engine or FF accuracy."""
    def setUp(self):
        self.geometry=analytic_ethane()
        self.assigned=assign_alkane(self.geometry)
        self.reference={'energy_components_kj_mol':dict.fromkeys(('bond','angle','dihedral','lj','coulomb','total'),0.),
                        'forces_kj_mol_angstrom':[[0.,0.,0.] for _ in self.geometry.atoms]}

    def fake_run(self,spec,force=.0):
        spec.run_directory.mkdir()
        rows=[(force,0.,0.) for _ in self.geometry.atoms]
        (spec.run_directory/'forces.dump').write_text(dump_text(self.geometry.coordinates_angstrom,rows))
        return LammpsRunResult({'status':'completed','thermo_blocks':[
            {'columns':['Step','Temp','E_bond','E_angle','E_dihed','E_vdwl','E_coul','PotEng'],
             'rows':[[0.,0.,0.,0.,0.,0.,0.,0.]]}]})

    def test_mock_parity_success_has_zero_steps_and_no_material_claim(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch('polychemmd.forcefields.reference.evaluate_alkane',return_value=self.reference),\
                 patch('polychemmd.engines.alkane.run_lammps',side_effect=self.fake_run):
                result=verify_alkane_singlepoint(self.assigned,run_directory=Path(directory)/'run')
        self.assertEqual(result['status'],'passed')
        self.assertEqual(result['md_integration_steps'],0)
        self.assertFalse(result['trajectory_production_performed'])
        self.assertEqual(result['force_field_physical_accuracy'],'not_assessed')

    def test_actual_force_mismatch_is_not_converted_into_success(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory)/'run'
            with patch('polychemmd.forcefields.reference.evaluate_alkane',return_value=self.reference),\
                 patch('polychemmd.engines.alkane.run_lammps',side_effect=lambda spec:self.fake_run(spec,force=1.)):
                with self.assertRaises(AlkaneParityError) as caught:
                    verify_alkane_singlepoint(self.assigned,run_directory=run)
            saved=json.loads((run/'parity_record.json').read_text())
        self.assertEqual(saved['status'],'failed')
        self.assertGreater(caught.exception.record['max_force_component_abs_error_kj_mol_angstrom'],4.)

    def test_invalid_tolerance_stops_before_reference_or_engine(self):
        with patch('polychemmd.engines.alkane.run_lammps') as engine:
            with self.assertRaises(ValueError):
                verify_alkane_singlepoint(self.assigned,run_directory='unused',energy_abs_tol_kj_mol=float('nan'))
            engine.assert_not_called()

    def test_extreme_finite_force_mismatch_keeps_finite_failure_report(self):
        with tempfile.TemporaryDirectory() as directory:
            run=Path(directory)/'run'
            with patch('polychemmd.forcefields.reference.evaluate_alkane',return_value=self.reference),\
                 patch('polychemmd.engines.alkane.run_lammps',side_effect=lambda spec:self.fake_run(spec,force=1e200)):
                with self.assertRaises(AlkaneParityError):
                    verify_alkane_singlepoint(self.assigned,run_directory=run)
            saved=json.loads((run/'parity_record.json').read_text())
        self.assertEqual(saved['status'],'failed')
        self.assertTrue(math.isfinite(saved['rms_force_component_error_kj_mol_angstrom']))
        self.assertGreater(saved['rms_force_component_error_kj_mol_angstrom'],1e200)

    def test_engine_failure_is_not_substituted(self):
        with tempfile.TemporaryDirectory() as directory:
            with patch('polychemmd.forcefields.reference.evaluate_alkane',return_value=self.reference),\
                 patch('polychemmd.engines.alkane.run_lammps',side_effect=LammpsError('unavailable fixture')):
                with self.assertRaisesRegex(LammpsError,'unavailable fixture'):
                    verify_alkane_singlepoint(self.assigned,run_directory=Path(directory)/'run')


if __name__=='__main__':
    unittest.main()
