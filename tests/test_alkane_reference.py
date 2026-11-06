import math
import unittest

from polychemmd.forcefields.reference import (
    COULOMB_KJ_MOL_ANGSTROM_E2, _angle_energy, _bond_energy,
    _coulomb_energy, _differentiate, _dihedral_energy, _lj_energy,
    evaluate_alkane,
)
from polychemmd.geometry import GeometryAtom, GeometryBond, MolecularGeometry


def butane_geometry():
    # Hand-specified nondegenerate coordinates. They are a derivative-check
    # fixture, not a minimized geometry or a molecular-property observation.
    carbons = ((0.,0.,0.),(1.53,0.,0.),(2.10,1.4,0.),(3.5,1.7,.6))
    offsets = (((-.4,.8,.6),(-.4,-.8,.6),(-.4,0.,-1.)),
               ((.1,-.6,.9),(.1,-.6,-.9)),
               ((-.4,.4,.9),(-.4,.4,-.9)),
               ((.6,.7,.5),(.6,-.7,.5),(.6,0.,-.9)))
    atoms = [GeometryAtom(f'C{i}',i+1,'C',0,0,12.011,'{}') for i in range(4)]
    bonds = [GeometryBond((f'C{i}',f'C{i+1}'),1) for i in range(3)]
    positions = list(carbons)
    for carbon,(position,hydrogens) in enumerate(zip(carbons,offsets)):
        for offset in hydrogens:
            name=f'H{len(atoms)-4}'
            atoms.append(GeometryAtom(name,len(atoms)+1,'H',0,0,1.008,'{}'))
            bonds.append(GeometryBond((f'C{carbon}',name),1))
            positions.append(tuple(a+b for a,b in zip(position,offset)))
    return MolecularGeometry(tuple(atoms),tuple(bonds),tuple(positions),'0'*64,
                             '{"source_kind":"hand_specified_derivative_fixture"}')


class LocalDerivativeTests(unittest.TestCase):
    def test_harmonic_bond_exact_energy_and_opposing_gradient(self):
        energy,gradient=_differentiate(((0.,0.,0.),(2.,0.,0.)),
                                      lambda p:_bond_energy(p,{'k_kj_mol_angstrom2':10.,'r0_angstrom':1.}))
        self.assertEqual(energy,10.)
        self.assertEqual(gradient,((-20.,0.,0.),(20.,0.,0.)))

    def test_right_angle_gradient_matches_independent_coordinate_derivative(self):
        points=((1.,0.,0.),(0.,0.,0.),(0.,1.,0.))
        energy,gradient=_differentiate(points,lambda p:_angle_energy(p,{'k_kj_mol_rad2':2.,'theta0_rad':math.pi/3}))
        self.assertAlmostEqual(energy,2*(math.pi/6)**2)
        self.assertAlmostEqual(gradient[0][1],-2*math.pi/3)
        self.assertAlmostEqual(gradient[2][0],-2*math.pi/3)
        for axis in range(3):
            self.assertAlmostEqual(math.fsum(row[axis] for row in gradient),0.,places=12)

    def test_lj_zero_crossing_minimum_and_half_scale(self):
        energy,gradient=_differentiate(((0.,0.,0.),(2.,0.,0.)),
                                      lambda p:_lj_energy(p,.5,2.,.5))
        self.assertEqual(energy,0.)
        self.assertAlmostEqual(gradient[1][0],-3.)
        r=2**(1/6)*2
        energy,gradient=_differentiate(((0.,0.,0.),(r,0.,0.)),
                                      lambda p:_lj_energy(p,.5,2.,.5))
        self.assertAlmostEqual(energy,-.25)
        self.assertAlmostEqual(gradient[1][0],0.,places=12)

    def test_coulomb_exact_energy_and_charge_force_sign(self):
        charge_product=-.02
        energy,gradient=_differentiate(((0.,0.,0.),(2.,0.,0.)),
                                      lambda p:_coulomb_energy(p,charge_product,.5))
        self.assertAlmostEqual(energy,COULOMB_KJ_MOL_ANGSTROM_E2*charge_product*.5/2)
        self.assertAlmostEqual(gradient[1][0],-COULOMB_KJ_MOL_ANGSTROM_E2*charge_product*.5/4)
        self.assertGreater(gradient[1][0],0.)  # force is negative: attractive.
        self.assertAlmostEqual(gradient[0][0],-gradient[1][0])

    @staticmethod
    def points(phi):
        return ((0.,1.,0.),(0.,0.,0.),(1.,0.,0.),(1.,math.cos(phi),math.sin(phi)))

    def test_complete_opls_profile_uses_trans_pi_and_half_prefactors(self):
        terms=[{'k_kj_mol':v*4.184/2,'periodicity':n,'phase_rad':0. if n%2 else math.pi}
               for n,v in enumerate((1.3,-.05,.2,0.),1)]
        expected=((0.,1.5*4.184),(math.pi,0.),(math.pi/3,.9375*4.184))
        for phi,value in expected:
            with self.subTest(phi=phi):
                energy,_=_differentiate(self.points(phi),lambda p:_dihedral_energy(p,{'terms':terms}))
                self.assertAlmostEqual(energy,value,places=12)

    def test_signed_arbitrary_phase_torsion_derivative(self):
        for phi in (.7,-.7):
            with self.subTest(phi=phi):
                values={'terms':[{'k_kj_mol':2.,'periodicity':2,'phase_rad':.4}]}
                energy,gradient=_differentiate(self.points(phi),lambda p:_dihedral_energy(p,values))
                derivative=-4*math.sin(2*phi-.4)
                self.assertAlmostEqual(energy,2*(1+math.cos(2*phi-.4)),places=12)
                self.assertAlmostEqual(gradient[3][1],-math.sin(phi)*derivative,places=12)
                self.assertAlmostEqual(gradient[3][2],math.cos(phi)*derivative,places=12)

    def test_singular_geometries_reject_undefined_forces(self):
        with self.assertRaises(ValueError):
            _differentiate(((0.,0.,0.),(0.,0.,0.)),lambda p:_bond_energy(p,{'k_kj_mol_angstrom2':1.,'r0_angstrom':1.}))
        with self.assertRaises(ValueError):
            _differentiate(((0.,0.,0.),(1.,0.,0.),(2.,0.,0.)),lambda p:_angle_energy(p,{'k_kj_mol_rad2':1.,'theta0_rad':1.}))
        with self.assertRaises(ValueError):
            _differentiate(((0.,0.,0.),(1.,0.,0.),(2.,0.,0.),(2.,1.,0.)),
                           lambda p:_dihedral_energy(p,{'terms':[{'k_kj_mol':1.,'periodicity':1,'phase_rad':0.}]}))


class AssignedReferenceTests(unittest.TestCase):
    def setUp(self):
        from polychemmd.forcefields.alkanes import assign_alkane
        self.geometry=butane_geometry()
        self.assigned=assign_alkane(self.geometry)

    def test_all_forces_match_central_finite_difference_energy_derivatives(self):
        reference=evaluate_alkane(self.assigned)
        step=2e-6
        positions=[list(p) for p in self.geometry.coordinates_angstrom]
        for i in range(len(positions)):
            for axis in range(3):
                plus,minus=([list(p) for p in positions] for _ in range(2))
                plus[i][axis]+=step
                minus[i][axis]-=step
                ep=evaluate_alkane(self.assigned,coordinates_angstrom=plus)['energy_components_kj_mol']['total']
                em=evaluate_alkane(self.assigned,coordinates_angstrom=minus)['energy_components_kj_mol']['total']
                numerical=-(ep-em)/(2*step)
                force=reference['forces_kj_mol_angstrom'][i][axis]
                self.assertTrue(math.isclose(force,numerical,rel_tol=2e-7,abs_tol=2e-6),
                                (i,axis,force,numerical))
        self.assertFalse(reference['finite_difference_forces'])

    def test_force_torque_and_rigid_rotation_covariance(self):
        reference=evaluate_alkane(self.assigned)
        forces=reference['forces_kj_mol_angstrom']
        for axis in range(3):
            self.assertAlmostEqual(math.fsum(f[axis] for f in forces),0.,places=9)
        torque=[0.,0.,0.]
        for (x,y,z),(fx,fy,fz) in zip(self.geometry.coordinates_angstrom,forces):
            for axis,value in enumerate((y*fz-z*fy,z*fx-x*fz,x*fy-y*fx)):
                torque[axis]+=value
        for value in torque:
            self.assertAlmostEqual(value,0.,places=8)
        moved=[(-p[1]+3,p[0]-2,p[2]+4) for p in self.geometry.coordinates_angstrom]
        rotated=evaluate_alkane(self.assigned,coordinates_angstrom=moved)
        self.assertAlmostEqual(rotated['energy_components_kj_mol']['total'],reference['energy_components_kj_mol']['total'],places=9)
        for old,new in zip(forces,rotated['forces_kj_mol_angstrom']):
            for actual,expected in zip(new,(-old[1],old[0],old[2])):
                self.assertAlmostEqual(actual,expected,places=8)

    def test_cutoff_only_changes_nonbonded_terms_and_ambiguity_is_rejected(self):
        complete=evaluate_alkane(self.assigned)
        tiny=evaluate_alkane(self.assigned,cutoff_angstrom=.01)
        self.assertEqual(tiny['energy_components_kj_mol']['lj'],0.)
        self.assertEqual(tiny['energy_components_kj_mol']['coulomb'],0.)
        for name in ('bond','angle','dihedral'):
            self.assertEqual(tiny['energy_components_kj_mol'][name],complete['energy_components_kj_mol'][name])
        data=self.assigned.to_dict()
        pair=next(p for p in data['nonbonded_pairs'] if p['lj_scale'])
        ids={a.atom_id:i for i,a in enumerate(self.geometry.atoms)}
        cutoff=math.dist(*(self.geometry.coordinates_angstrom[ids[name]] for name in pair['atom_ids']))
        with self.assertRaisesRegex(ValueError,'discontinuous hard cutoff'):
            evaluate_alkane(self.assigned,cutoff_angstrom=cutoff)

    def test_invalid_positions_cutoffs_and_type_fail_without_result(self):
        for cutoff in (True,0.,-1.,float('nan'),10**1000,'12'):
            with self.subTest(cutoff_type=type(cutoff)),self.assertRaises(ValueError):
                evaluate_alkane(self.assigned,cutoff_angstrom=cutoff)
        for positions in ([],[[0,0,0]],[[0,0]]*14,[[False,0,0]]*14,[[float('inf'),0,0]]*14):
            with self.subTest(positions_type=type(positions)),self.assertRaises(ValueError):
                evaluate_alkane(self.assigned,coordinates_angstrom=positions)
        with self.assertRaises(ValueError):
            evaluate_alkane(self.assigned.to_dict())

    def test_result_keeps_model_state_and_method_boundaries(self):
        result=evaluate_alkane(self.assigned)
        self.assertEqual(result['atom_ids'],[a.atom_id for a in self.geometry.atoms])
        self.assertEqual(result['boundary_conditions'],'isolated_nonperiodic')
        self.assertEqual(result['material_validation'],'not_performed')
        self.assertEqual(result['equilibrium_status'],'not_assessed')
        self.assertFalse(result['trajectory_integration_performed'])
        self.assertEqual(set(result['energy_components_kj_mol']),{'bond','angle','dihedral','lj','coulomb','total'})
        components=result['energy_components_kj_mol']
        self.assertAlmostEqual(components['total'],sum(v for k,v in components.items() if k!='total'))


if __name__=='__main__':
    unittest.main()
