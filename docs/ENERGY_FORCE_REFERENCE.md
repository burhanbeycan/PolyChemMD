# Isolated OPLS energy and analytical force reference

`forcefields/reference.py` independently evaluates the explicitly assigned
GROMACS v2025.3 OPLS-AA n-alkane numerical subset. It is a small verification
kernel, not an MD integrator, equilibrium sampler or bulk electrostatics method.
Its source does not import or copy an external force-field evaluator.

```python
from polychemmd.forcefields.alkanes import assign_alkane
from polychemmd.forcefields.reference import evaluate_alkane

assigned = assign_alkane(geometry)
result = evaluate_alkane(assigned, cutoff_angstrom=12.0)
```

The output retains atom order/IDs, coordinates actually evaluated, separate
bond/angle/dihedral/LJ/Coulomb energy contributions and total forces. Molecular
energies are kJ/mol; forces are kJ/(mol angstrom). An optional
`coordinates_angstrom` override retains the assigned atom ordering and changes
only the positions used for this evaluation. It does not alter the assignment.

## Explicit potential and derivatives

Bond and angle energies use the project's no-half convention:

\[
E_b=K_b(r-r_0)^2,\qquad E_a=K_a(\theta-\theta_0)^2.
\]

All canonical Fourier terms are retained:

\[
E_d=\sum_n K_n[1+\cos(n\phi-\delta_n)].
\]

The signed dihedral follows adjacent oriented plane normals, with cis equal to
zero and trans equal to plus/minus pi. The local implementation uses `atan2`
and therefore retains explicit phases rather than assuming an unsigned angle.
The pinned OPLS subset has phases zero/pi, including its negative second
coefficient. The OPLS one-half prefactor has already been converted into each
canonical Fourier amplitude by the parameter-data layer.

For eligible distinct atom pairs:

\[
E_{LJ}=s_{LJ}\,4\epsilon_{ij}
[(\sigma_{ij}/r)^{12}-(\sigma_{ij}/r)^6],\qquad
E_q=s_q\,Cq_iq_j/r.
\]

Both sigma and epsilon mix geometrically. Explicit graph distances one/two are
excluded, graph distance three scales both terms by one half, and larger graph
distances have unit scaling. Pair records come from the complete assigned graph
and are counted once. All hydrogen-ended proper torsions remain present.

The reference is isolated/nonperiodic, with relative dielectric constant one.
Unshifted LJ and Coulomb terms are included only for `r < cutoff_angstrom`;
there is no periodic image, reciprocal-space sum, tail correction, switching,
reaction field or energy shift. The Coulomb conversion is exactly
`332.06371 * 4.184` kJ angstrom/(mol e^2), matching the numerical convention in
the pinned LAMMPS `units real` implementation. This numerical engine convention
is distinguished from an exact modern SI-derived physical constant.

Forces are minus the analytical coordinate derivatives of these equations.
The original implementation carries each interaction's value and local
coordinate derivatives through forward-mode arithmetic. A bonded interaction
has at most twelve differentiation coordinates; gradients are accumulated into
the full atom array. No finite-difference step is used to calculate the
reported forces. Automatic differentiation is an established mathematical
method, not a claimed new scientific algorithm.

## Numerical gates and checks

- Coordinate arrays must match the atom count, contain exactly three finite
  numeric components per atom and reject booleans.
- Zero geometric norms, coincident eligible nonbonded points and collinear or
  nearly collinear angular/torsional geometries are rejected. The relative
  sine threshold is `1e-12`; no arbitrary force is substituted at singularities.
- Nonrepresentable energies or derivatives fail explicitly.
- Hard-cutoff energies are discontinuous. A nonbonded separation within
  `max(1e-10 angstrom, 1e-12 * cutoff)` of the cutoff is rejected instead of
  silently selecting an ambiguous derivative.

Executed tests independently specify harmonic derivatives, a right-angle
coordinate derivative, LJ zero-crossing/minimum and half scaling, Coulomb force
sign, complete trans/gauche OPLS torsion energies, and a signed arbitrary-phase
torsion derivative. Full assigned C4H10 force components are also checked
against central finite differences of energy on a hand-specified nondegenerate
coordinate fixture. That derivative check uses a `2e-6` angstrom step and
predeclared relative/absolute tolerances `2e-7` and `2e-6` respectively. It is a
test of analytical derivatives, not the method used to produce forces.
Rigid-rotation covariance, vanishing total force/torque and invalid/cutoff
inputs are checked separately.

## Scientific boundary

Engine agreement establishes convention, implementation and export parity for
the stated isolated verification model. It does not establish force-field
physical accuracy, condensed-phase long-range electrostatics, polymer density,
equilibration or transport. Historical OPLS-AA and later OPLS variants must not
be mixed silently. Bulk protocols need separate parameter-applicability,
electrostatics, cutoff, finite-size, relaxation and sampling assessment.

## Primary equation and convention references

- [LAMMPS harmonic bonds](https://docs.lammps.org/bond_harmonic.html)
- [LAMMPS harmonic angles](https://docs.lammps.org/angle_harmonic.html)
- [LAMMPS OPLS dihedrals](https://docs.lammps.org/dihedral_opls.html)
- [LAMMPS LJ/Coulomb cut model](https://docs.lammps.org/pair_lj_cut_coul.html)
- [Pinned LAMMPS real-unit conversion](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L177)
- [GROMACS 2025.3 bonded forms and RB/Fourier conventions](https://manual.gromacs.org/documentation/2025.3/reference-manual/functions/bonded-interactions.html)

These references define the equations and engine conventions; the project's
explicit implementation and its executed tests supply implementation evidence.
