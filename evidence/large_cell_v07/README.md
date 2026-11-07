# Executed larger-cell integration and sensitivity audit — v0.7

The bounded integration completed with **passed_engine_integration**. The separate
numerical-sensitivity review remains **requires_review**.
This is an actual 2304-atom construction/parameter/engine/restart check. It is
not a liquid-phase validation or an equilibrium density prediction.

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_large_cell_v07.py --output-dir artifacts/new-large-cell-audit
```

Use an editable source checkout, the pinned working engine and a new output
directory. The runner declares inputs, budgets and review limits before invoking
any engine. Runtime source hashes and actual artifact hashes remain in
`report.json`. The final run took **19.433876717 s** on
its recorded single-rank CPU host; this is not a GPU or other-computer benchmark.

## One dense geometric construction

The starting motif is the first explicit C10H22 geometry from
`evidence/geometry_v03/cell.json`. The unchanged motif is tiled 2×6×6 at
13×5.1×5.1 angstrom spacings. The cell contains **72 molecules / 2304 atoms**
in a **26×30.6×30.6 angstrom** box. All identities, source lineage, charges,
parameters, bonded topology, wrapped/unwrapped positions and images remain
explicit.

Its construction density is **0.6987592396767455 g/cm³**, computed from supplied
molar masses and requested volume. That value is a construction input/result,
not a density predicted at 300 K and 1 bar. The ordered repeated motif is not
identified as a crystal, liquid or amorphous polymer phase. Seventy-two molecules
and all numerical variants represent **one prepared cell**, not 72 independent
cell replicates. Repeating one conformation does not sample conformational
entropy or establish an equilibrium molecular ensemble.

An independent contact oracle checks 26 neighboring motif translations and
all 32×32 atom pairs for each translation: 26624 distance evaluations. The
minimum periodic intermolecular distance is **1.6099405152823842 angstrom**.
Every farther lattice offset has a coordinate-span lower bound of
**6.063975729010847 angstrom**, so it cannot lower that minimum. This check is
independent of the builder's linked spatial bins. The 1.5 angstrom contact floor
is geometric, not a van der Waals or force-equilibrium criterion.

The compact assignment stores **12384 special-pair records**, with **2653056**
all-unordered pairs represented by explicit default rules plus exceptions.
There is no global `nonbonded_pairs` array. The exact canonical compact record
is 9326465 bytes and remains with local raw artifacts; this is a measured
serialization size, not a speedup or asymptotic-memory benchmark. Per-molecule
assignment work remains explicitly bounded. The same frozen OPLS-AA parameter
model and short-topology rules are retained.

## Actual engine stages and fixed numerical choices

| Stage | Actual status | Measured wall time (s) |
| --- | --- | ---: |
| baseline | completed_singlepoint | 2.050967 |
| cutoff_12 | completed_singlepoint | 1.601355 |
| reciprocal_20 | completed_singlepoint | 1.211193 |
| preparation | completed_preparation_pilot | 8.598165 |
| continuation | completed_npt_continuation | 4.356109 |

The three singlepoints use the identical geometry and assignment:

- Baseline: 10 angstrom cutoff, Ewald alpha 0.30 inverse angstrom, reciprocal index 16.
- Cutoff sensitivity: only the cutoff increases to 12 angstrom.
- Reciprocal sensitivity: only the reciprocal index increases to 20.

Every box axis is strictly longer than twice the largest cutoff. The model uses
geometric LJ mixing, the fixed OPLS exclusions/1–4 scaling, periodic tin-foil
Ewald and no LJ tail correction. Large-cell singlepoint execution does not
include an independent analytical large-cell reference comparison. Separate
small-cell tests check compact/legacy equivalence, generated inputs and the
bounded independent Ewald implementation.

Preparation used the explicitly resource-calibrated CG budget of 10 iterations /
100 evaluations, followed by 10 NVT and 10 NPT steps at 0.25 fs. CG stopped at
**10 iterations / 18 evaluations: budget_exhausted**, without force convergence.
A further 10-step native NPT continuation preserved the active fix state,
coordinates, velocities, images and step. Native topology, coefficients and
units were independently read back by the continuation audit. There are **30 MD
steps / 7.5 fs** in the entire path; minimizer iterations add no MD duration.

## Numerical sensitivity remains a review outcome

| Comparison | Absolute total-energy difference / atom (kJ/mol) | Maximum force-component difference (kJ/mol/angstrom) | Review |
| --- | ---: | ---: | --- |
| Cutoff 10 → 12 angstrom | 0.04105997831218941 | 0.0073084167087529295 | requires_review |
| Reciprocal index 16 → 20 | 0 | 1.4210854715202004e-14 | diagnostics_passed |

The signed cutoff changes are -94.6441132497821 kJ/mol in LJ energy and
+0.04192321849131986 kJ/mol in real-space Coulomb energy, giving
-94.60219003128441 kJ/mol total. Bonded terms and unchanged reciprocal terms
remain within their predeclared invariance checks.

The cutoff energy-per-atom difference exceeds the **predeclared 0.02 kJ/mol**
development review margin; the force difference is below its 0.05 margin.
The reciprocal margins were 1e-5 kJ/mol per atom and 1e-4 kJ/mol/angstrom.
These thresholds were not widened after execution. Per-atom energy means the
cell energy difference divided by atom count; it is not an experimental molar
material property. The cutoff comparison changes the truncated Hamiltonian,
so this is sensitivity rather than error against an exact model. Agreement at
two finite reciprocal choices does not prove convergence of an infinite Ewald
sum. Bulk-property and cutoff qualification remain open before production.

## Unrelaxed final state and insufficient sampling

The final instantaneous state at global step 30 has temperature
**127.63842142870227 K**, pressure **7097.816097594815 bar** and canonical density
**0.6987018285061131 g/cm³**. The raw native density is
0.6987017670143656 g/cm³; the engine and exact-SI conversion conventions remain
separate. These values describe this short unrelaxed state, not the requested
300 K / 1 bar equilibrium averages.

The continuation contains two frames spanning **0.0025 ps**. Default physical
duration and complete/coarsened block gates fail, returning
`insufficient_sampling`, `production_ready=false` and
`density_prediction_validation=not_performed`. No density uncertainty or
independent-cell error bar is manufactured. Independent cells, observed phase,
adequate relaxation/sampling, cutoff/tail and size checks, and state-matched
external evidence remain necessary.

## Retained resource failure and public artifact policy

`failed_attempt.json`, `failed_preparation.in` and `failed_preparation.log`
preserve the first actual attempt and its hashes. The original CG request of
200 iterations / 1000 evaluations hit its 60-second preparation timeout after
the log reached iteration 30, before any MD. The full failed audit took
65.794599895 s, including its completed singlepoints. It remains failed.

A **new recipe was declared before the next engine invocation**: CG maximum
10 iterations / 100 evaluations, with the same 10 NVT + 10 NPT + 10 continuation
steps, physical model, numerical choices, comparison margins, five-call count,
60-second per-call cap and 240-second aggregate audit cap. This resource
calibration enables a bounded technical probe; it does not establish minimizer
convergence or repair the failed cutoff review.

Public files include the motif geometry and fixed recipe, one shared
`common/system.data`, each stage's input/log, three singlepoint force dumps and
the final continuation's velocity-bearing `checkpoint.dump`. Repeated large
arrays are omitted from summary records with canonical field hashes. Native
binary checkpoints, the full compact record and other large native read-backs
remain local, with explicit raw hashes/sizes and this regeneration recipe.
The public set is approximately 2.014 MB; it contains single snapshots rather
than a large trajectory. Source and artifact hashes were verified after the
actual run. The broad independent polymer platform remains under development.
