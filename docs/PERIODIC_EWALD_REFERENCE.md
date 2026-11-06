# Explicit periodic Ewald energy and force reference

This independent verification kernel consumes the strict, frozen
`AssignedPeriodicAlkaneCell` result. It supports a neutral orthorhombic cell with
three periodic axes and the selected OPLS-AA n-alkane parameters. It does not
integrate trajectories, infer equilibrium, establish bulk parameter accuracy or
certify convergence of a selected Ewald setup.

```python
from polychemmd.forcefields.periodic_reference import evaluate_periodic_alkanes

result = evaluate_periodic_alkanes(
    assigned,
    cutoff_angstrom=10.0,
    ewald_alpha_per_angstrom=0.30,
    reciprocal_nmax=16,
)
```

Molecular energies are kJ/mol and forces kJ/(mol angstrom). The numerical Coulomb
conversion is `332.06371 * 4.184`, explicitly matching the pinned LAMMPS
`units real` convention rather than claiming an exact SI-derived constant.
Bond, angle and complete Fourier-torsion terms reuse the project's analytical
local-differentiation kernels. Real/reciprocal nonbonded forces are independent
analytical coordinate derivatives of the equations below. Finite differences
are used only for separate derivative tests.

## Periodic geometry and real-space terms

Bonded terms use the complete molecular **unwrapped** coordinates. Every bond
component must be strictly shorter than half its corresponding cell axis.
Real-space nonbonded terms use orthorhombic minimum-image displacements and an
explicit hard cutoff strictly below half the shortest box side. The model uses
unshifted LJ with geometric sigma/epsilon mixing and no dispersion tail term.
Its fixed graph-based pair scales remain zero for 1–2/1–3, one half for 1–4 and
one for other/intermolecular pairs.

For Ewald real space, the screened term is present for every atom pair inside
the cutoff, including excluded/scaled bonded pairs:

\[
E_{\rm screened}=C\sum_{i<j,r_{ij}<r_c}q_iq_j\,
\frac{\operatorname{erfc}(\alpha r_{ij})}{r_{ij}}.
\]

Reciprocal space already contains all charges, so exclusions require a separate
full-Coulomb correction rather than multiplication of the screened term by the
special scale:

\[
E_{\rm exclusion}=-C\sum_{\rm special\ pairs}(1-s_{ij})\frac{q_iq_j}{r_{ij}}.
\]

Every such short pair must use the same molecular and minimum image, and must
lie inside the real cutoff. Ambiguous images are rejected instead of applying a
scale to the wrong periodic pair. The main reported `coulomb_realspace`
component is the sum of screened and exclusion terms. Python `math.erfc` is
used; the reference does not imitate a backend-specific erfc approximation.
Consequently an engine's polynomial approximation can produce a small real-space
residual even when the mathematical model and reciprocal calculation agree.

## Reciprocal and self terms

Let `k = 2*pi*(nx/Lx, ny/Ly, nz/Lz)` and
`S(k) = sum_i q_i exp(i k dot r_i)`. The finite reference sum is

\[
E_{\rm reciprocal}=\frac{2\pi C}{V}
\sum_{\substack{n_x,n_y,n_z=-N\\\mathbf n\ne0}}^{N}
\frac{\exp[-k^2/(4\alpha^2)]}{k^2}|S(\mathbf k)|^2,
\qquad E_{\rm self}=-\frac{C\alpha}{\sqrt\pi}\sum_iq_i^2.
\]

The code evaluates one member of each `k/-k` pair with the exactly doubled
prefactors. Reciprocal forces use the analytical derivative of the same
structure factors; the self term has zero positional derivative. Conducting
(tin-foil) macroscopic boundary conditions mean no surface term. A neutralizing
background is not silently introduced: nonneutral cells are unsupported.
`coulomb_kspace` is reciprocal plus self energy, while detailed screened,
exclusion, reciprocal and self values are retained separately to prevent
double counting in the main energy-component sum.

`reciprocal_nmax` defines a **cube of integer vectors**. The cited pinned LAMMPS
Ewald implementation additionally screens vectors using a spherical squared-
wavevector limit. Those finite lists are not generally identical, particularly
at small N. The report includes floating-point triangle expressions bounding
the energy/force difference caused by the cube's additional vectors, using the
actual charge magnitudes and the cited spherical limit. These bounds concern
only those two finite vector lists. They do **not** bound the omitted infinite
reciprocal tail, real-space cutoff error, numerical backend approximations or
physical-model accuracy. The default N=16/24-angstrom-cube example makes this
particular list difference very small; different cells/settings need their own
assessment.

## Inputs, budgets and tests

The reference is limited to 512 atoms, integer N from 1 through 24 and at most
20 million atom-mode evaluations by default. The reciprocal budget is checked
before bonded/pair work. The optional explicit budget may be changed to another
positive integer. Undefined norms, nonrepresentable reciprocal scales, squared
wavevector underflow/overflow and nonfinite energies/forces fail explicitly.
Pairs on the hard cutoff are rejected within the documented numerical band
`max(1e-10 angstrom, 1e-12*cutoff)`; no arbitrary force is substituted there.

An optional `coordinates_angstrom` argument evaluates new **unwrapped** molecular
coordinates in the same atom order. It does not rewrite the assignment, source
geometry or images. Whole-cell translations and whole-molecule shifts by lattice
vectors preserve energy/forces; independently wrapping atoms without their
bonded images is not an equivalent unwrapped input and can fail the image gate.

Tests include independently enumerated reciprocal vectors/prefactors, known LJ
and screened-Coulomb derivatives, finite reciprocal convergence, Ewald splitting
parameter consistency for an analytic neutral point-charge fixture, periodic
translations and molecule-image shifts, net force, complete component sums and
selected total-force components against central energy differences. Analytic
point-charge fixtures are mathematical checks, not material data. The periodic
system tests use explicitly labeled specified ethane geometry, not a measured
or relaxed thermodynamic state.

## Primary equation and backend sources

- [GROMACS Ewald equations and convergence context](https://manual.gromacs.org/documentation/current/reference-manual/functions/long-range-electrostatics.html#ewald-summation)
- [LAMMPS kspace styles](https://docs.lammps.org/kspace_style.html)
- [LAMMPS explicit gewald/kmax settings](https://docs.lammps.org/kspace_modify.html)
- [Pinned engine finite-vector screen](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/KSPACE/ewald.cpp#L291-L310)
- [LAMMPS LJ and long-range Coulomb pair equations](https://docs.lammps.org/pair_lj_cut_coul.html)

The implementation here is original project code for established equations.
Its numerical checks are distinct from production sampling and independent
experimental/property validation.
