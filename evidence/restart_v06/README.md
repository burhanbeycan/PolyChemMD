# Executed native NPT restart audit — v0.6

`report.json` records an actual bounded LAMMPS execution, including a same-cell
40-step continuation versus two 20-step continuations. It also reviews the
actual continuous branch with the density sampling diagnostics. Restart parity
passed; sampling remains **insufficient**. No equilibrium density or liquid
phase is validated.

Run from an editable source checkout with the pinned working engine:

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_restart_v06.py --output-dir artifacts/new-restart-audit
```

The output directory must be new. This command regenerates native checkpoints;
the committed evidence deliberately contains their hashes rather than binary
files. Native restart compatibility is restricted to the matching engine build,
platform and recorded configuration. Linked system libraries are not
exhaustively fingerprinted.

## Declared system and finite resources

The source is `evidence/geometry_v03/cell.json`: four explicit C10H22 molecules,
128 atoms, in a 24 angstrom cube. The frozen OPLS-AA subset uses 10 angstrom
real-space cutoff, Ewald alpha 0.30 inverse angstrom and reciprocal index 16.
The requested controls are 300 K and 1 bar, with 0.25 fs timesteps. These are
control settings, not evidence that the instantaneous state matches them.

The auditor writes its declared input, tolerances and budget before invoking
the engine: four calculation processes, at most 60 s requested per process
and 240 s aggregate audit wall budget. The initial CG minimization permits
200 iterations and 1000 evaluations, followed by 100 NVT and 100 NPT steps.
Continuations preserve the active NPT fix state, atom identities, velocities,
images, box and parameter identity. They do not reseed, minimize or reset time.
An initial run-zero in each continuation audits the restored state; its frame
is excluded from the sampling series.

The final audit took **8.886093162 s** on the executing Linux x86-64 host
(Intel Xeon Platinum 8573C; one MPI rank, one OpenMP thread, no GPU requested).
This is measured container runtime, not a benchmark of another computer.
LAMMPS identified itself as **22 Jul 2025 — Update 4**. Exact engine identity,
native package hashes, source hashes, input hashes and artifact hashes are
retained in the records.

| Actual stage | MD steps executed | Final global step | Measured stage wall time (s) |
| --- | ---: | ---: | ---: |
| Preparation | 200 | 200 | 5.436452 |
| Continuous branch | 40 | 240 | 1.142131 |
| First split fragment | 20 | 220 | 1.099150 |
| Second split fragment | 20 | 240 | 1.120971 |

These executions total 280 steps because alternative branches repeat the same
physical path. The longest individual path is 240 steps, **60 fs including
preparation**. Minimizer iterations add no MD duration. There is **one** prepared
cell; the branches, restart fragments and frames are not independent replicates.

## Predeclared endpoint comparison

| Quantity | Absolute tolerance | Observed maximum absolute difference |
| --- | ---: | ---: |
| Position component (angstrom) | 1e-9 | 1.0658141036401503e-14 |
| Velocity component (angstrom/fs) | 1e-10 | 1.93421667571414e-15 |
| Force component (kJ/mol/angstrom) | 1e-6 | 2.5501378786429996e-11 |
| Box component (angstrom) | 1e-9 | 0 |
| Potential/total energy (kJ/mol) | 1e-7 | 5.684341886080801e-13 |

Atom, molecule, fixed-charge and image identities match exactly. Both wrapped
and unwrapped positions are checked; minimum-image equivalence does not hide
a discontinuity. Actual initial/final text snapshots and data read-backs make
these checks reviewable. Native binary hashes and checkpoint parent identities
retain the restart lineage.

## Sampling outcome and scope

The continuous branch has five thermo frames at steps 200, 210, 220, 230 and
240: **10 fs = 0.01 ps**, with a 2.5 fs sampling interval. Canonical density is
recomputed from supplied molar masses and actual thermo volume using the exact
SI Avogadro constant. Pressure is converted from native atmospheres to bar.
The native engine density remains a separate field.

The predeclared default review requires 100 ps and complete statistical blocks.
The actual window fails physical duration, complete-block and coarsened-block
gates, with too few frames for variability diagnostics. No standard error or
independent-cell uncertainty is fabricated. The five frames and correlated
restart branches cannot be promoted into independent observations.

The dilute four-molecule fixture also cannot establish a liquid under its
10 angstrom cutoff guard: keeping every box axis above 20 angstrom prevents
arbitrary compression. A separate size/cutoff-qualified condensed-phase system,
independent preparations, adequate sampling and state-matched external evidence
remain necessary before reporting a bulk density prediction.

## Retained failed attempts

`failed_attempt.json` and `failed_attempt_2.json` preserve the initial actual
failed outcomes, including their source/input hashes and raw-report hash.
Neither failure was erased or reclassified as a successful calculation.

The first failed during strict parameter read-back: native angle coefficients
are written by the pinned LAMMPS `AngleHarmonic::write_data` using C `%g`, which
defaults to six significant digits. For example, the native coefficient
58.350023900573618 is printed as 58.35. The corrected audit compares the exact
declared text serialization while separately checking the binary hash and
restored force behavior. It does not relax the scientific parity tolerances.
See the [pinned primary source](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/MOLECULE/angle_harmonic.cpp).

The second failed because `write_data` assigns new interaction-row numbers
when it emits bonds, angles and proper dihedrals. Those row numbers are not
global atom identities. The corrected audit compares the complete typed
topology multiset, preserving multiplicities and all atom endpoints under the
supported harmonic/OPLS full-reversal symmetry. Missing, duplicated or changed
interactions still fail. Global atom and molecule identities remain strict.

The successful rerun followed these audited read-back fixes with the original
predeclared endpoint tolerances. This evidence verifies restart implementation
for this bounded case; it does not establish force-field physical accuracy,
polymer relaxation, numerical bulk-property convergence or experimental agreement.
