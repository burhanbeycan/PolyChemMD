# Isolated OPLS-AA n-alkane single points

TASK-006B connects a supported chemical assignment to an actual LAMMPS
single-point energy/force calculation. The exporter, parser and comparison
logic are independently written. LAMMPS is an external numerical engine; no
other polymer workflow package is imported or copied.

```python
from polychemmd.forcefields.alkanes import assign_alkane
from polychemmd.engines.alkane import (
    build_alkane_input, prepare_alkane_singlepoint,
    verify_alkane_singlepoint,
)

assigned = assign_alkane(geometry)
script, data = build_alkane_input(assigned, cutoff_angstrom=12.0)
spec = prepare_alkane_singlepoint(
    assigned, run_directory="artifacts/new-alkane-run", executable="lmp",
)
# The two calls above create no files and execute no engine.
report = verify_alkane_singlepoint(
    assigned, run_directory="artifacts/another-new-run", executable="lmp",
    energy_abs_tol_kj_mol=1e-6,
    force_abs_tol_kj_mol_angstrom=1e-5,
)
```

The run directory must be new. The existing bounded process adapter identifies
the executable, stages `system.data`, records input hashes and rejects failed,
stale, empty or incomplete output. A separately prepared environment may need
`UCX_TLS=self OMP_NUM_THREADS=1` for the pinned serial-container smoke setup;
this does not validate MPI/GPU execution on another host.

## Exported model and unit conversion

The only supported model is the explicitly assigned n-alkane subset documented
in [ALKANE_ASSIGNMENT.md](ALKANE_ASSIGNMENT.md). The export uses:

- `units real`, `atom_style full`, fixed nonperiodic boundaries, one molecule,
  and a box extending beyond every atom by more than the specified cutoff.
- Harmonic bonds and angles with canonical `K * displacement^2`. The canonical
  kJ/mol coefficients are divided by 4.184 for LAMMPS kcal/mol; an extra factor
  of one half is **not** applied. Angle reference values become degrees, while
  the coefficient remains per radian squared.
- Explicit `dihedral_style opls` coefficients. A canonical Fourier term becomes
  `K_n = 2*k_n/4.184`; its phase must be zero for odd periodicity and pi for even
  periodicity. Unsupported phases/periodicities fail instead of silently changing
  the potential. Every assigned proper dihedral is exported.
- Unshifted `lj/cut/coul/cut`, relative dielectric one, geometric sigma/epsilon
  mixing, no tail correction, and `special_bonds lj/coul 0 0 0.5`. The latter
  excludes 1-2/1-3 interactions and applies one-half scaling to both 1-4 terms.

This is a deliberately isolated, finite-cutoff comparison model. It is not a
periodic bulk electrostatics, Ewald/PPPM, density or production sampling protocol.
Using a larger cutoff does not by itself validate condensed-phase accuracy.

## Actual output and parity checks

`run 0` evaluates the current coordinates without integration or minimization.
Thermo normalization is disabled and float precision is requested explicitly.
The parser requires one completed step-zero thermo row and one force dump with
exactly the requested atoms. It checks finite values, unique sequential engine
IDs, nonperiodic box fields and coordinate identity. Output rows are mapped back
to the original atom order; engine/source atom identities are retained in the
report. Lost atoms, duplicate IDs, extra frames, changed coordinates or malformed
fields produce explicit failures.

Bond, angle, dihedral, LJ, Coulomb and total potential energies are compared
separately with the independently written analytical-gradient reference.
Every Cartesian force component is also compared. LAMMPS kcal/mol and
kcal/(mol angstrom) are multiplied by 4.184 before comparison. The reference
uses the same documented numerical Coulomb conversion as the pinned engine:
`332.06371 kcal angstrom/(mol e^2)`. This numerical convention is explicit,
not an undocumented change in the physical model.

Absolute tolerances are checked before execution and retained in the report.
They are not adjusted after seeing the result. `AlkaneParityError` retains a
failed comparison record, which is also written as `parity_record.json` in the
run directory. Even a very large finite mismatch is reported with finite,
scaled RMS arithmetic. Unavailable or failed engines never produce a substituted
reference-as-engine success.

A passing report establishes implementation/convention agreement for the stated
coordinates, model and numerical tolerances. It does **not** establish parameter
transferability, a relaxed geometry, equilibrium, a polymer material property,
an experimental prediction or a production-ready simulation. The report keeps
zero MD integration steps and explicit unassessed material/equilibrium fields.

## Pinned-backend near-planar diagnostic

The pinned LAMMPS 22 Jul 2025 Update 4 `dihedral_style opls` source floors the
sine used in torsional derivative ratios when its magnitude is below `1e-5`.
Near a planar torsion, this numerical treatment can produce force differences
larger than floating-point roundoff while the potential energies still agree.
For example, the fixed decane verification geometry contains a nearly trans
CCCC torsion and exhibits a force-component difference of approximately
`5.1e-6 kJ/(mol angstrom)`, within the predeclared `1e-5` tolerance. That result
must not be described as machine-precision force agreement.

`backend_numerical_diagnostics` records a geometric screen for assigned
dihedrals with `abs_sin_phi < 1e-5`, listing their interaction IDs, atom IDs and
computed absolute sine. The screen uses project-written unit-plane geometry;
it does not copy or substitute the engine force algorithm. Its threshold is
identified as an observed detail of the cited pinned backend, not a guarantee
about other engine versions/builds. The actual executable banner is retained
separately. Source:
https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/MOLECULE/dihedral_opls.cpp#L171-L180

The diagnostic is informational. It does not change the potential, independent
reference forces, input seeds or acceptance tolerances, and it never turns a
failed comparison into a pass. Exact and near-trans analytical geometry cases
test the screen independently of engine execution.

## Verification and primary sources

Unit tests cover independently known conversion values, complete ethane export,
strict dump identities, malformed/poisoned outputs and failure gates. Synthetic
engine/reference fixtures are labeled as process-comparison tests; they are not
evidence of actual LAMMPS execution. Actual integration records are separate.

- LAMMPS harmonic bond convention:
  https://docs.lammps.org/bond_harmonic.html
- LAMMPS harmonic angle convention:
  https://docs.lammps.org/angle_harmonic.html
- LAMMPS OPLS Fourier convention:
  https://docs.lammps.org/dihedral_opls.html
- LAMMPS LJ/Coulomb cutoff model:
  https://docs.lammps.org/pair_lj_cut_coul.html
- LAMMPS special bonded-pair scaling:
  https://docs.lammps.org/special_bonds.html
- Exact pinned-engine real-unit numerical constant (`qqr2e`):
  https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L177

The parameter source snapshot and its exact numerical subset are recorded by
the assignment layer. These references support methods and conventions, not
PolyChemMD physical-property accuracy.
