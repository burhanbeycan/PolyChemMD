# Executed numerical and preparation review — v0.8

The final actual audit completed **eight engine requests in 73.651719068 s**.
Engine execution and state/read-back checks passed. Its independent outcomes
remain **`numerical_status: requires_review`** and
**`preparation_readiness.status: requires_preparation`**. These outcomes were
retained with the original declared thresholds. No equilibrium density, phase
identity or external material accuracy is validated.

Reproduce from this editable checkout with the pinned working engine:

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_qualification_v08.py --output-dir artifacts/new-qualification-v08
```

The output directory must be new. The runner writes its exact recipe, review
specifications, input identities and budgets before requesting any engine work.
It retains failed/inconclusive records rather than substituting a successful
number. The runtime source hashes in `report.json` match the final implementation.

## Declared geometry and resource scope

The input is one stored C10H22 geometry from
`evidence/large_cell_v07/common/source_geometry.json`, repeated on an explicitly
ordered grid using `examples/qualification_cell.json`:

- 3 x 7 x 7 copies with spacings 13 x 5.1 x 5.1 angstrom.
- **147 molecules / 4704 explicit atoms**, box **39 x 35.7 x 35.7 angstrom**.
- Construction density **0.6987592396767457 g/cm^3**, from specified mass/volume.
- Independent motif-image contact minimum **1.609940515 angstrom**, above the
  unchanged 1.5 angstrom geometric floor. A finite neighboring-offset calculation
  and a bound for farther offsets establish this geometry result independently
  of the builder's linked bins.
- The largest 17 angstrom cutoff is strictly less than half every box axis.

This is one ordered repeated-motif configuration, with an unassigned phase.
Neither the 147 copies, six numerical settings nor continuation fragments are
independent material-cell replicates. Construction density is not a prediction.

Compact assignment stores **25,284 special pairs**, with **11,061,456 global
unordered pairs** represented by explicit default rules rather than a global
pair array. The full canonical assignment is 19,107,924 bytes locally and is
retained by hash/reproduction recipe instead of being duplicated in Git.

The predeclared limits are eight requests, 60 s per request and 240 s aggregate
audit wall time. The numerical subprotocol has its own 120 s aggregate limit;
the parent refuses to start it without that entire allowance plus overhead.
Actual execution used one direct engine process and one OpenMP thread on a
Linux x86-64 Intel Xeon Platinum 8573C host. No GPU was requested. These measured
container timings are not a benchmark of the project owner's computer.
LAMMPS identified itself as 22 Jul 2025 — Update 4.

| Actual stage | Engine elapsed time (s) |
| --- | ---: |
| Cutoff 10 angstrom, reciprocal bound 16 | 0.954534 |
| Cutoff 12 angstrom, reciprocal bound 16 | 0.977264 |
| Cutoff 14 angstrom, reciprocal bound 16 | 1.073103 |
| Cutoff 16 angstrom, reciprocal bound 16 | 1.236769 |
| Cutoff 17 angstrom, reciprocal bound 16 | 1.215704 |
| Cutoff 17 angstrom, reciprocal bound 20 | 1.482111 |
| Bounded preparation | 26.786631 |
| Native continuation | 8.608059 |

The numerical phase took **19.031732843 s** including its orchestration. Complete
preparation/continuation calls took 28.102404721 and 11.966016206 s. Per-point
wrapper wall times were not separately measured; the report uses null rather
than relabeling engine elapsed times. The complete audit also includes source,
construction, assignment, planning and postprocessing work.

## Fixed-geometry numerical evidence

All six evaluations retain the same source geometry, assignment, atom order,
parameter model and producer identity. They use zero velocities and an explicit
**static virial pressure** compute. Scalar/tensor pressure is in bar with
positive compression; it is not an equilibrium pressure average or mechanical
stress tensor. The Ewald splitting parameter is 0.30 inverse angstrom, with
geometric LJ mixing, no LJ tail correction, no Coulomb lookup table and the same
frozen OPLS-AA numerical subset.

The cutoff criteria were declared before execution: energy change per atom at
most **0.02 kJ/mol per atom**, maximum force-component change at most
**0.05 kJ/(mol angstrom)**, and scalar **and every tensor pressure-component**
change at most **20 bar**. **Two consecutive finest refinement pairs** must
pass. These are development review settings, not universal physical error bounds.

| Cutoff change (angstrom) | Absolute energy change / atom (kJ/mol) | Maximum force change (kJ/mol/angstrom) | Static scalar pressure change (bar) | Maximum tensor change (bar) | Review |
| --- | ---: | ---: | ---: | ---: | --- |
| 10 to 12 | 0.0410599783 | 0.00730841671 | 128.816535 | 138.454167 | requires review |
| 12 to 14 | 0.0191982380 | 0.00552827110 | 60.327102 | 65.759950 | requires review |
| 14 to 16 | 0.0112543475 | 0.00125250042 | 35.368892 | 36.201920 | requires review |
| 16 to 17 | 0.00394312592 | 0.000433551858 | 12.392564 | 12.772921 | diagnostics passed |

Only the last refinement pair passes every gate. In particular, the preceding
14-to-16 pair passes energy/force limits but fails pressure limits. The required
two consecutive pairs therefore do not pass. **No settings are automatically
selected** (`selected_settings` is null), and no threshold was widened.
Refining a truncated cutoff changes its Hamiltonian; these differences are not
errors against an exact infinite-cutoff or bulk reference.

At fixed cutoff 17 angstrom, reciprocal bound 16 to 20 produced:

- Energy change per atom **7.39954371397e-11 kJ/mol**.
- Maximum force-component change **6.18456397206e-10 kJ/(mol angstrom)**.
- Scalar static pressure change **1.74962769961e-6 bar**.
- Maximum tensor-component change **2.30985460803e-6 bar**.

The reciprocal diagnostic passes its separately predeclared limits. Comparing
two finite bounds does not establish the infinite Ewald error or convergence of
a material property. Bonded-component invariants also passed the comparisons.

## Preparation and native continuation

The 17 angstrom / reciprocal-20 preparation probe was declared before the sweep
for technical integration/readiness evidence. It is not automatic production
execution under numerically qualified settings. Its physical controls were
300 K, 1 bar, 0.25 fs, CG max10/eval100, 10 NVT and 10 NPT steps, followed by ten
native continuation steps.

CG stopped at **max iterations**, after ten iterations and 18 evaluations.
Its force tolerance was not met. Native continuation preserved the saved NPT
fix state, positions, velocities, images and box; full topology, parameter and
unit read-back checks passed. There were **30 MD steps / 7.5 fs** in total.
Minimizer iterations do not add MD duration.

The preparation readiness result retains these failed gates:

- Minimizer force-convergence policy.
- Liquid-candidate orientation screen.
- NPT density sampling adequacy.
- Starting-configuration numerical qualification.

The final configuration also has no new numerical-sensitivity audit. A review of
the initial immutable geometry cannot be transferred to a thermally deformed
state. Readiness always remains weaker than a phase/equilibrium certificate.

Preparation-snapshot maximum and RMS atom-force norms were **361.201817829** and
**90.3336634833 kJ/(mol angstrom)**. They are instantaneous thermal forces, not
minimization residuals. Optional gross-force thresholds are null in the declared
review, so these values are reported without inventing a universal cutoff.
The actual minimizer's own stopping criterion supplies the minimization gate.

The continuation has two thermo frames spanning **0.0025 ps** and returns
`insufficient_sampling`. Final instantaneous values were 130.012647968 K,
7219.11241323 bar and canonical density 0.6987032421592682 g/cm^3. These are
execution outputs, not ensemble averages or validated predictions. The short
physical duration and declared gates leave preparation/state sampling open.

## Static chain-structure outputs

The new structure calculation preserves its chemical definitions: all-atom
mass-weighted Rg, carbon-backbone Rg, terminal-carbon end-to-end Re, gyration
tensors/eigenvalues/relative shape anisotropy, and an equal-molecule apolar
orientation tensor. The largest orientation eigenvalue is the nematic S.

| Snapshot | Mean all-atom Rg (angstrom) | Mean backbone Rg (angstrom) | Mean terminal-carbon Re (angstrom) | Nematic S |
| --- | ---: | ---: | ---: | ---: |
| Constructed | 3.068135950 | 2.972163706 | 8.574463350 | 1.000000000 |
| Prepared | 3.191134812 | 3.098150643 | 9.049731954 | 0.999970104 |
| Continued | 3.191867841 | 3.098711925 | 9.050717256 | 0.999940608 |

The nearly parallel axes retain the ordered preparation bias; S is not a phase
classifier. These are statistics of individual static configurations, not
converged ensemble observables. No uncertainty is estimated from them.
`continued_chain_structure.json` contains the complete 147-chain continued
snapshot dataset. Earlier duplicated per-chain rows are hash-identified in
summaries and reproducible with the runner; aggregate values, definitions and
units remain visible. See [chain structure](../../docs/CHAIN_STRUCTURE.md) and
[preparation review](../../docs/PREPARATION_READINESS.md).

## Retained initial failure and resource-only correction

`failed_attempt.json` retains the actual first run, which took 24.276107122 s.
Four static points completed, then the 17 angstrom point failed before
preparation with **LAMMPS neighbor-list overflow**. The default allowance was
2000 neighbors per atom; the list distance included the unchanged 2 angstrom
skin. The periodic geometric cutoff bound remained valid. Partial input/log
text for all five attempted points is retained in `failed_inputs/` with hashes.

The corrected implementation explicitly allocates neighbor storage for larger
cells: `one=max(2000,2*N)` and `page=max(100000,10*one)` for N above 512. For this
case, those values are **9408 and 100000**. The same control is restored during
native continuation. Small-cell input bytes remain unchanged. The force field,
physical cutoffs, skin, neighbor rebuilding policy and numerical review limits
were not modified. The final audit was a new declared run with fresh artifacts.
See the primary [LAMMPS neighbor storage documentation](https://docs.lammps.org/neigh_modify.html).

## Public raw evidence and reproducibility

The committed subset contains the fixed geometry/recipe, one shared
`common/system.data`, every stage input/log/summary, the finest reciprocal
single-point force frame, final velocity-bearing checkpoint frame, and the
complete continued static-chain dataset. Public numerical input scripts expect
`system.data` in their working directory; use the runner to stage it correctly,
or copy the identical shared file for manual inspection/execution.

Other full frames, native read-backs, the complete assignment and binary
checkpoints are retained locally by explicit SHA256/size and can be regenerated
with the runner. They are not silently presented as committed downloads.
Native checkpoints require the matching engine build/platform/configuration.
Fifteen source-file hashes, 25 public raw-artifact references and the full
continued-structure dataset identity were verified against the recorded bytes.
The earlier failed attempt is historical evidence with its own source hashes.

The next scientific requirements remain numerical pressure/cutoff qualification,
appropriate alternative preparation histories, independently prepared cells,
adequate phase/state sampling, finite-size assessment and a state-matched
external density comparison.
