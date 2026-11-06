import math
import unittest
from unittest.mock import patch

from polychemmd.builders.packing import pack_cell
from polychemmd.forcefields.periodic import assign_periodic_alkanes
from polychemmd.forcefields.periodic_reference import (
    COULOMB_CONSTANT, _minimum_image, _pair_terms, _reciprocal_ewald,
    evaluate_periodic_alkanes,
)
from test_alkane_lammps import analytic_ethane


class EwaldKernelTests(unittest.TestCase):
    def test_minimum_image_across_faces(self):
        actual = _minimum_image((23.7, -23.6, 47.5), (24.,24.,24.))
        for x,y in zip(actual,(-.3,.4,-.5)):
            self.assertAlmostEqual(x,y,places=13)

    def test_screened_real_force_and_exclusion_match_independent_derivative(self):
        r,q,a=2.3,-.04,.31
        _,real,correction,force=_pair_terms((r,0.,0.),q,0.,1.,0.,.5,a)
        self.assertAlmostEqual(real,COULOMB_CONSTANT*q*math.erfc(a*r)/r,places=13)
        self.assertAlmostEqual(correction,-.5*COULOMB_CONSTANT*q/r,places=13)
        expected=(COULOMB_CONSTANT*q/r**2)*(
            math.erfc(a*r)+2*a*r/math.sqrt(math.pi)*math.exp(-(a*r)**2)-.5)
        self.assertAlmostEqual(force[0],expected,places=13)
        step=1e-6
        def energy(x):
            return COULOMB_CONSTANT*q*(math.erfc(a*x)-.5)/x
        self.assertAlmostEqual(force[0],-(energy(r+step)-energy(r-step))/(2*step),places=8)

    def test_lj_known_zero_crossing_and_geometric_force_sign(self):
        # Zero charges isolate a known half-scaled LJ derivative.
        lj,real,correction,force=_pair_terms((2.,0.,0.),0.,.5,2.,.5,1.,.3)
        self.assertEqual((lj,real,correction),(0.,0.,0.))
        self.assertAlmostEqual(force[0],3.)
        r=2*2**(1/6)
        lj,_,_,force=_pair_terms((r,0.,0.),0.,.5,2.,.5,1.,.3)
        self.assertAlmostEqual(lj,-.25)
        self.assertAlmostEqual(force[0],0.,places=12)

    def test_one_mode_reciprocal_force_matches_independent_energy_gradient(self):
        points=((1.,1.2,2.),(3.,2.4,1.))
        charges=(.2,-.2);box=(10.,11.,12.)
        energy,forces,metadata=_reciprocal_ewald(points,charges,box,.4,1)
        # Independently enumerate all 26 k-vectors with the full-sum prefactor.
        expected=0.
        for nx in (-1,0,1):
            for ny in (-1,0,1):
                for nz in (-1,0,1):
                    if (nx,ny,nz)==(0,0,0):
                        continue
                    k=tuple(2*math.pi*n/l for n,l in zip((nx,ny,nz),box))
                    k2=sum(v*v for v in k)
                    phases=[sum(v*x for v,x in zip(k,p)) for p in points]
                    sr=sum(q*math.cos(p) for q,p in zip(charges,phases))
                    si=sum(q*math.sin(p) for q,p in zip(charges,phases))
                    expected+=2*math.pi*COULOMB_CONSTANT/math.prod(box)*math.exp(-k2/(4*.4**2))/k2*(sr*sr+si*si)
        self.assertAlmostEqual(energy,expected,places=12)
        self.assertEqual(metadata['full_vector_count'],26)
        step=1e-6
        for atom,axis in ((0,0),(0,2),(1,1)):
            plus,minus=([list(p) for p in points] for _ in range(2))
            plus[atom][axis]+=step;minus[atom][axis]-=step
            ep=_reciprocal_ewald(plus,charges,box,.4,1)[0]
            em=_reciprocal_ewald(minus,charges,box,.4,1)[0]
            self.assertAlmostEqual(forces[atom][axis],-(ep-em)/(2*step),places=8)

    def test_finite_reciprocal_convergence_and_split_parameter_identity(self):
        # Analytic neutral point charges, not fitted alkane/material data.
        points=((1.,1.,1.),(3.,1.,1.));charges=(.2,-.2);box=(10.,10.,10.)
        values=[_reciprocal_ewald(points,charges,box,.6,n)[0] for n in (6,10,14)]
        self.assertGreater(abs(values[0]-values[-1]),1e-6)
        self.assertLess(abs(values[1]-values[-1]),1e-8)
        totals=[]
        for alpha in (.6,.8):
            reciprocal=_reciprocal_ewald(points,charges,box,alpha,16)[0]
            real=COULOMB_CONSTANT*charges[0]*charges[1]*math.erfc(2*alpha)/2
            self_energy=-COULOMB_CONSTANT*alpha/math.sqrt(math.pi)*sum(q*q for q in charges)
            totals.append(real+reciprocal+self_energy)
        self.assertAlmostEqual(*totals,places=7)

    def test_finite_cube_sphere_difference_is_bounded(self):
        points=((1.,1.,1.),(3.,2.,1.));charges=(.2,-.2);box=(10.,10.,10.)
        _,_,metadata=_reciprocal_ewald(points,charges,box,.6,2)
        bound=metadata['cited_backend_finite_list_difference']
        self.assertGreater(bound['energy_triangle_bound_kj_mol'],0.)
        self.assertGreater(bound['max_force_component_triangle_bound_kj_mol_angstrom'],0.)
        self.assertIn('not an infinite-sum',bound['scope'])

    def test_reciprocal_work_and_invalid_inputs_fail_early(self):
        for n in (True,0,25,1.2):
            with self.subTest(n=n),self.assertRaises(ValueError):
                _reciprocal_ewald(((0,0,0),(1,0,0)),(.2,-.2),(24,24,24),.3,n)
        with self.assertRaisesRegex(ValueError,'budget'):
            _reciprocal_ewald(((0,0,0),(1,0,0)),(.2,-.2),(24,24,24),.3,4,100)
        with self.assertRaises(ValueError):
            _pair_terms((0.,0.,0.),-.01,.1,2.,1.,1.,.3)
        for box in ((24.,24.,1e300),(1e-300,1e150,1e150)):
            with self.subTest(box=box),self.assertRaises(ValueError):
                _reciprocal_ewald(((0,0,0),(1,0,0)),(.2,-.2),box,.3,1)


class PeriodicAssignedReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cell=pack_cell((analytic_ethane(),)*2,box_lengths_angstrom=(24.,24.,24.),
                           seed=671,min_distance_angstrom=1.5).to_dict()
        cls.assigned=assign_periodic_alkanes(cls.cell)
        cls.positions=[list(a['coordinates_angstrom']) for a in cls.assigned.to_dict()['atoms']]

    def evaluate(self,**kwargs):
        return evaluate_periodic_alkanes(self.assigned,reciprocal_nmax=4,**kwargs)

    def test_components_identity_and_scope(self):
        result=self.evaluate()
        components=result['energy_components_kj_mol']
        self.assertEqual(set(components),{'bond','angle','dihedral','lj','coulomb_realspace','coulomb_kspace','total'})
        self.assertAlmostEqual(components['total'],math.fsum(v for k,v in components.items() if k!='total'),places=12)
        self.assertEqual(result['atom_ids'],[a['cell_atom_id'] for a in self.cell['atoms']])
        self.assertEqual(result['md_integration_steps'],0)
        self.assertEqual(result['equilibrium_status'],'not_assessed')
        self.assertEqual(result['material_property_validation'],'not_performed')
        self.assertFalse(result['finite_difference_forces'])

    def test_selected_total_force_components_match_finite_differences(self):
        result=self.evaluate();step=2e-6
        for atom,axis in ((0,0),(0,2),(3,1),(8,0),(12,2)):
            plus,minus=([list(p) for p in self.positions] for _ in range(2))
            plus[atom][axis]+=step;minus[atom][axis]-=step
            ep=self.evaluate(coordinates_angstrom=plus)['energy_components_kj_mol']['total']
            em=self.evaluate(coordinates_angstrom=minus)['energy_components_kj_mol']['total']
            expected=-(ep-em)/(2*step)
            actual=result['forces_kj_mol_angstrom'][atom][axis]
            self.assertTrue(math.isclose(actual,expected,abs_tol=2e-6,rel_tol=2e-7),(atom,axis,actual,expected))

    def test_global_translation_and_whole_molecule_image_shift_invariance(self):
        baseline=self.evaluate()
        for move in ((3.7,-8.2,24.1),(24.,-48.,72.)):
            shifted=[[x+d for x,d in zip(p,move)] for p in self.positions]
            translated=self.evaluate(coordinates_angstrom=shifted)
            for name,value in baseline['energy_components_kj_mol'].items():
                self.assertAlmostEqual(value,translated['energy_components_kj_mol'][name],places=10)
            for first,second in zip(baseline['forces_kj_mol_angstrom'],translated['forces_kj_mol_angstrom']):
                for x,y in zip(first,second):
                    self.assertAlmostEqual(x,y,places=10)
        shifted=[list(p) for p in self.positions]
        for i in range(8):
            shifted[i][0]+=24.
        result=self.evaluate(coordinates_angstrom=shifted)
        self.assertAlmostEqual(result['energy_components_kj_mol']['total'],baseline['energy_components_kj_mol']['total'],places=10)

    def test_net_force_and_special_correction_are_explicit(self):
        result=self.evaluate()
        for axis in range(3):
            self.assertAlmostEqual(math.fsum(f[axis] for f in result['forces_kj_mol_angstrom']),0.,places=10)
        parts=result['coulomb_decomposition_kj_mol']
        self.assertNotEqual(parts['special_pair_exclusion'],0.)
        self.assertEqual(parts['surface'],0.)
        self.assertEqual(parts['neutralizing_background'],0.)
        self.assertAlmostEqual(parts['reciprocal']+parts['self'],result['energy_components_kj_mol']['coulomb_kspace'],places=12)

    def test_invalid_cutoffs_geometry_and_budget_fail(self):
        for cutoff in (0.,True,float('nan'),12.,20.):
            with self.subTest(cutoff=cutoff),self.assertRaises(ValueError):
                self.evaluate(cutoff_angstrom=cutoff)
        bad=[list(p) for p in self.positions];bad[0][0]+=24
        with self.assertRaisesRegex(ValueError,'bonded unwrapped'):
            self.evaluate(coordinates_angstrom=bad)
        bad=[list(p) for p in self.positions];bad[0][0]=float('inf')
        with self.assertRaises(ValueError):
            self.evaluate(coordinates_angstrom=bad)
        with self.assertRaisesRegex(ValueError,'budget'):
            with patch('polychemmd.forcefields.periodic_reference._differentiate') as differentiate:
                try:
                    self.evaluate(max_reciprocal_atom_evaluations=1)
                finally:
                    differentiate.assert_not_called()


if __name__=='__main__':
    unittest.main()
