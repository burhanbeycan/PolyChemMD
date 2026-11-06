# Periodic alkane Ewald export and numerical parity

`engines/periodic_alkane.py` consumes a validated
`AssignedPeriodicAlkaneCell`. It exports the selected OPLS-AA numerical subset
into a fully periodic orthorhombic LAMMPS cell and compares an actual run-zero
single point with the independent periodic reference. It does not infer a
material phase or validate a production simulation.

```python
from polychemmd.engines.periodic_alkane import (
    build_periodic_alkane_setup, build_periodic_alkane_input,
    prepare_periodic_alkane_singlepoint, verify_periodic_alkane_singlepoint,
)

report = verify_periodic_alkane_singlepoint(
    assigned,
    run_directory="artifacts/new-periodic-check",
    cutoff_angstrom=10.0,
    ewald_alpha_per_angstrom=0.30,
    reciprocal_nmax=16,
    energy_abs_tol_kj_mol=1e-4,
    force_abs_tol_kj_mol_angstrom=1e-4,
)
```

## Scientific setup and identity

The model explicitly selects `units real`, `atom_style full`, three periodic
axes, harmonic bonds/angles, complete OPLS proper dihedrals, geometric sigma
and epsilon mixing, and `special_bonds lj/coul 0 0 0.5`. LJ is unshifted with
no tail correction. Real-space Coulomb uses `lj/cut/coul/long`, dielectric one,
with the pair table disabled. Reciprocal Coulomb uses `kspace_style ewald 1e-8`
plus explicit `kspace_modify gewald alpha kmax/ewald n n n`. This is a separate
periodic model, rather than silently reusing the isolated truncated Coulomb
setup from v0.4.

Atom rows contain wrapped positions and all three integer image flags, with
one-based molecule IDs derived from the source molecule index. All global and
source IDs are retained in comparison records. Data files include complete
mass, pair, bond, angle and dihedral coefficients in converted real units.

The cutoff must be strictly less than half the shortest box dimension. For
every excluded/scaled pair at graph distance one through three, the supplied
unwrapped displacement must be its minimum-image displacement and its distance
must lie inside the real-space cutoff. This restricts ambiguous bonded images
and guarantees the selected Ewald exclusion correction is represented by the
same short pair in the engine and reference.

`build_periodic_alkane_setup()` returns `(setup_script, data)` ending after the
scientific model/read-data/neighbor setup. It contains no thermo, integration,
minimization or output commands. Bounded protocol builders reuse that helper.
`build_periodic_alkane_input()` appends single-point thermo, `run 0`, a force
dump and completion marker. Preparation functions create no files or processes.

## Output and numerical comparison

The parser requires exactly one step-zero dump with the periodic box and fields
`id mol x y z ix iy iz fx fy fz`. Missing/duplicate atom IDs, changed molecule
IDs/images, altered box/positions, extra frames and nonfinite data fail
explicitly. Forces are mapped back into the source atom sequence.

Thermo energies use absolute system normalization (`norm no`) and 17-digit
formatting. Bond, angle, dihedral, LJ, real-space Coulomb, reciprocal/self
Coulomb and total energies are compared separately. `E_coul` maps to the
reference real-space screened term plus special-pair correction; `E_long` maps
to reciprocal plus self energy. All Cartesian forces are compared after
kcal-to-kJ conversion. Tolerances are validated before running and are never
adjusted after observing an output. Failed comparisons persist a finite JSON
failure record and raise `AlkaneParityError`.

Both implementations receive the same alpha and integer nmax settings. The
independent reference enumerates an integer reciprocal cube, whereas the pinned
LAMMPS Ewald implementation screens a sphere. The finite lists are therefore
not asserted to be identical for arbitrary settings. At the stated 24 angstrom
box, alpha 0.30/angstrom and nmax16 verification settings, the omitted mode
weights are negligible; lower nmax or larger boxes still need the stated
reference bound and actual parity gate. A solver accuracy keyword by itself is
not evidence of convergence of the complete physical model.

## Pinned backend numerical effects

The report retains the v0.4 near-planar OPLS sine diagnostic. Its reference
force is not changed to imitate the pinned engine derivative safeguard.

The pinned `lj/cut/coul/long` implementation also evaluates a polynomial erfc
approximation and a rounded `2/sqrt(pi)` value in its real-space kernel, even
with pair tables disabled. The reference uses `math.erfc` and the direct
mathematical constant. `realspace_backend_numerical_context` records those
sources and their version scope; it does not change the potential, reference,
geometry or acceptance tolerances.

An actual 128-atom starting-cell comparison showed real-space energy difference
about 6.97e-5 kJ/mol, reciprocal/self difference about 4e-13 kJ/mol and maximum
force-component difference about 2.59e-5 kJ/(mol angstrom), inside the declared
1e-4 limits. The real-space energy residual was quantitatively explained by
the pinned erfc approximation. The force residual was dominated by one nearly
trans CCCC torsion flagged by the OPLS diagnostic. These are numerical-model
checks; they are not roundoff agreement for every observable or evidence of
accurate polymer density.

## Primary sources

- [LAMMPS kspace settings](https://docs.lammps.org/kspace_modify.html)
- [LAMMPS Ewald style](https://docs.lammps.org/kspace_style.html)
- [Wrapped positions and image flags](https://docs.lammps.org/read_data.html)
- [LJ/Coulomb long-range pair style](https://docs.lammps.org/pair_lj_cut_coul.html)
- [Pinned real-space erfc kernel](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/KSPACE/pair_lj_cut_coul_long.cpp#L136-L141)
- [Pinned erfc constants](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/ewald_const.h#L19-L25)
- [Pinned Ewald reciprocal mode screen](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/KSPACE/ewald.cpp)

The actual source snapshot, preparation and parameter provenance accompany the
assignment and execution records. No trajectory integration, equilibration,
bulk property, force-field transferability or experimental accuracy is
validated by this single-point adapter.
