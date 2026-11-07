# Preparation readiness review

`analysis/readiness.py` evaluates one full, completed preparation record against
predeclared development criteria. It launches no engine and writes no files.
Its scope is the supported neutral explicit n-alkane model and an explicitly
unverified `liquid_candidate` target. A successful process is different from
adequate preparation, a sampled phase or a validated density prediction.

```python
from polychemmd.analysis.readiness import (
    PreparationReadinessSpec, analyze_preparation_readiness,
)

spec = PreparationReadinessSpec.from_mapping(specification_mapping)
review = analyze_preparation_readiness(
    assigned,
    full_preparation_record,
    spec,
    numerical_review=full_starting_geometry_numerical_review,
)
```

`examples/preparation_readiness.json` records the explicit defaults. Select and
justify a study's policy before inspecting its results. This function does not
tune thresholds or change a protocol to produce a passing review.

## Record and unit consistency

The original assignment must be a validated legacy or compact periodic object.
Its canonical SHA256 and model must match the preparation. The complete protocol
specification must retain the requested temperature, pressure and timestep and
match the readiness declaration. The nested `sampling_spec` uses those same
physical controls; differing sampling units or targets cannot pass silently.

Completed records require actual process completion, the independent native
`write_data` header confirming real units, all three minimization/NVT/NPT raw
thermo blocks, and full NVT/NPT stage records. The analysis checks counts,
initial/final steps, physical durations and stage-boundary consistency. Raw
temperatures, native-atmosphere pressure and volume are compared with the
projected canonical rows. Canonical density is independently recomputed from
the supplied atomic molar masses and exact-SI Avogadro conversion; the raw
engine density uses the pinned real-unit constant separately. Raw box volume
and minimum-image cutoff bounds are checked too.

The final state must retain all atoms, source/molecule identities, fixed charges,
images and coordinates. The chain-structure analyzer revalidates these fields
and bonded image continuity. The final step, box volume, canonical density,
temperature and pressure must agree with the preparation's final NPT row.
Reduced public reports that intentionally omit atom arrays or raw blocks cannot
establish readiness. Use the full local execution record or reproduce it with
the corresponding auditor.

These are consistency checks on supplied records, not digital signatures or a
cryptographic authentication of external engine output. Preserve original file
hashes and execution provenance alongside the review.

## Minimizer policy and thermal forces are distinct

The default policy requires the actual minimizer's `force tolerance` stopping
criterion and consistent `force_tolerance_met` flag. An iteration/time budget or
energy-tolerance stop is not relabeled as force convergence. The caller may
prespecify `require_minimizer_force_convergence=false` for a scientifically
justified alternative protocol; doing so after inspecting a failure invalidates
a predeclared comparison. Even true force convergence does not establish a
minimum rather than another critical point, let alone thermal equilibration.
See the [primary LAMMPS minimize contract](https://docs.lammps.org/minimize.html).

After finite-temperature MD, instantaneous forces fluctuate. They must **not**
be required to approach zero as if they were a minimized-configuration residual.
The report therefore gives the maximum atom-vector norm and

\[
F_{\mathrm{RMS,atom}}=\sqrt{\frac1N\sum_i(F_{ix}^2+F_{iy}^2+F_{iz}^2)}.
\]

This is an RMS over atom-vector magnitudes, not over individual Cartesian
components. Finite-value checks are mandatory. Optional `max_force_norm_...`
and `max_rms_force_...` bounds default to null. If declared, they are labeled
gross-instability development screens, not universal physical limits or an
equilibrium test. The existing protocol does not retain minimizer force-norm
statistics separately, so their report fields remain explicitly unknown rather
than borrowing final thermal forces. NVT/NPT dynamics and their target controls
are described in the [primary engine documentation](https://docs.lammps.org/fix_nh.html).

## Structure and sampling screens

The chain analyzer reports static radius of gyration, end-to-end distance and
apolar orientation. For the liquid candidate, the default nematic-order screen
requires `S <= 0.20`. This is an illustrative preparation-bias diagnostic. It
does not classify a phase: low order alone cannot distinguish a liquid from
every other disordered state. A repeated ordered motif generally fails this
screen; copied molecules are one cell, not independent material samples.
Degenerate endpoints leave orientation inconclusive rather than receiving a
fabricated direction.

The actual uninterrupted NPT stage is passed to the existing density sampling
review. Its defaults require at least 100 ps plus complete/coarsened blocks,
autocorrelation closure, precision/drift and sample-mean state checks. A bounded
short preparation cannot satisfy that default duration. Later production-window
and continuation analysis must be qualified separately; this review does not
merge boundaries, discard rows or count restarts as independent cells. Its
uncertainty estimates remain conditional, as documented in
[density sampling](DENSITY_SAMPLING.md).

## A numerical review belongs to its source configuration

An optional numerical qualification must retain its full specification and
stage records, including the force arrays. The readiness analyzer reruns the
pure numerical comparison and requires the nested summary and top-level
status, identities, selected settings and producer identity to agree with that
derived result. Missing cases, rewritten passing labels, different atom/geometry
hashes, changed cutoff/alpha/reciprocal settings and a different producer engine
cannot establish the initial numerical gate. Public force-omitted summaries are
insufficient for this adoption step.

A passing fixed-geometry review applies only to the **original preparation
starting geometry**. Its cutoff, alpha and refined reciprocal bound must equal
the actual preparation controls. Evolved final coordinates and a thermal
trajectory require their own numerical assessment. The report always exposes
this missing transfer evidence; it cannot turn a starting-state numerical pass
into a production certificate.

## Outcomes

All gates and failures remain visible. `requires_preparation` indicates an
unmet minimizer policy, requested force screen, orientation screen or verified
starting numerical review. If those pass but sampling does not, the status is
`insufficient_sampling`. If every implemented local gate passes, the status is
only `diagnostics_passed_but_unqualified`.

Every outcome retains `production_ready=false`, `equilibrium_status=not_assessed`,
`phase_identification=not_performed` and no validated density prediction. Final
configuration/trajectory numerical sensitivity, independently prepared cells,
preparation-history sensitivity, phase and slow-mode/chain relaxation, finite
size and a state-matched external comparison remain unverified requirements.
A properly identified failed engine record is retained as a failed preparation
review with no invented force or sampling statistics.

Tests use actual archived short-run records for consistency/rejection checks and
explicit synthetic numerical and time-series fixtures for controlled local-gate
cases. Synthetic fixtures are not material simulations or physical evidence.
