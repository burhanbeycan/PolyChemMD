# Fixed-geometry numerical qualification

`protocols/numerics.py` implements an independently written, predeclared
cutoff/reciprocal refinement review for the supported periodic n-alkane model.
It consumes actual static thermodynamic single points from the new virial
adapter. Passing this review concerns **one fixed source geometry** and explicit
development thresholds. It never identifies a phase, qualifies a final thermal
configuration or validates an equilibrium material property.

```python
from polychemmd.protocols.numerics import (
    NumericalQualificationSpec,
    plan_numerical_qualification,
    run_numerical_qualification,
)

spec = NumericalQualificationSpec.from_mapping(explicit_specification)
plan = plan_numerical_qualification(
    assigned, spec, run_directory="runs/numerics",
)
# Only this call discovers/executes an external solver.
result = run_numerical_qualification(
    assigned, spec, run_directory="runs/numerics",
)
```

The run directory must be new. Planning creates no directory, reads no previous
run and discovers no executable. It generates every scientific script, input
hash and request before the first engine call. Common `system.data` is retained
once, and each planned stage includes its exact script and SHA256 values. All
cutoffs must satisfy the source box and short-topology minimum-image constraints
before execution can begin.

## Requests, identities and evidence

The default `examples/numerical_qualification.json` contains five ascending
cutoffs, 10, 12, 14, 16 and 17 angstrom, at Ewald alpha 0.30 inverse angstrom
and reciprocal index 16. A sixth evaluation uses the same finest 17 angstrom
cutoff and alpha with reciprocal index 20. No source coordinates, atom order,
charges, molecular identities, box, potential parameters or other numerical
controls change between these evaluations.

`compare_numerical_records(records, spec)` supplies the pure review kernel. It
requires records in the complete predeclared order, including coarse evaluations
that fail the thresholds. Missing or failed cases remain incomplete; they cannot
be discarded to produce a selected setting. Every provided record, including
partial failures, must contain finite JSON data and supported status/reason types.

Successful records must agree in their canonical assignment and geometry hashes,
model, ordered atom IDs, counts, density and volume. Recorded engine build and
configuration identities must be identical. Fixed controls are explicitly the
supported tin-foil Ewald model, requested relative accuracy 1e-8, no real-space
Coulomb table and no LJ tail correction. Canonical units, zero temperature and
kinetic energy, static zero-velocity output, identity audit, completed engine
record and an actual real-unit native-header proof must be retained. Numerical
review checks serialized consistency; these records are not cryptographic
signatures authenticating untrusted external calculations.

The geometry SHA256 projects the complete assignment `box` and ordered atoms'
`cell_atom_id`, `molecule_index`, `charge_e`, `coordinates_angstrom`,
`wrapped_coordinates_angstrom` and `image_flags`, serialized as finite UTF-8
JSON with sorted keys and compact separators. It describes the starting source;
its numerical result cannot silently transfer to later MD coordinates.

## Predeclared comparison measures

For each adjacent cutoff pair, the implementation records every signed component
energy change and evaluates `abs(delta total energy)/atom_count`, maximum
absolute force-component change, force-component RMS change, scalar static
pressure change and maximum pressure-tensor component change. Tensor checks keep
an anisotropic change from disappearing in the scalar trace. Static pressure
uses the virial-only adapter, not a thermal pressure average. See
[the static thermodynamic adapter](PERIODIC_THERMO.md) and the primary
[LAMMPS pressure-compute contract](https://docs.lammps.org/compute_pressure.html).

| Development gate | Cutoff refinement | Reciprocal refinement |
| --- | ---: | ---: |
| Absolute total energy change per atom, kJ/mol per atom | 0.02 | 1e-5 |
| Maximum force-component change, kJ/(mol angstrom) | 0.05 | 1e-4 |
| Scalar static pressure change, bar | 20 | 0.01 |
| Maximum static pressure-tensor component change, bar | 20 | 0.01 |

Bond, angle and dihedral energies must remain unchanged within the separately
predeclared absolute 1e-7 kJ/mol tolerance. Reciprocal energy is additionally
invariant during cutoff variation; LJ and real-space Coulomb energies are
invariant during reciprocal variation. Violating those component invariants
prevents selection even if a small total change happens to cancel errors.

At least **two consecutive finest cutoff pairs** must pass, and the reciprocal
comparison must pass. All coarse outcomes remain visible. The report also
retains whether absolute adjacent total-energy changes are monotone; it does
not erase a nonmonotone sequence or promote that finite trend into a rigorous
error bound. The default selection, when all required gates pass, is the finest
actually evaluated cutoff (17 angstrom) and refined reciprocal index (20), with
the original alpha.

These thresholds are transparent development review settings, not universal
published convergence constants. They must be justified and declared for the
study before its outcomes are examined. They must not be widened after a failed
comparison. Changing a real-space cutoff changes the truncated Hamiltonian;
the differences are sensitivities, not errors against a known exact bulk model.
Two finite reciprocal settings likewise do not bound the infinite Ewald error.
LJ tails, finite-size effects, force-field accuracy and observable-specific
thermal sampling remain separate requirements.

## Finite resources and retained failures

The default request allows six engine requests, at most 60 seconds per stage
and 180 seconds of accounted protocol wall time. Hard supported ceilings are
120 seconds per stage and 240 aggregate seconds. Callers may explicitly choose
smaller budgets; the full requested stage count must fit the declared maximum.
Planning time is included in aggregate accounting. Before each engine call,
the remaining budget reduces its requested timeout, with a two-second reserve.

The underlying process adapter bounds engine identification and the calculation
subprocess. This is not a real-time scheduler that can atomically preempt every
Python/filesystem operation. Aggregate time is checked after stages, comparison
postprocessing and final record serialization. An overrun cannot return a
successful protocol result; its failed record and completed/partial evaluations
are preserved. A stale existing root is never overwritten.

`engine_requests_started` counts adapter requests, including an unavailable
executable. `engine_execution_attempts_reported` separately counts stages whose
actual process record reports `simulation_performed=true`; unavailable engines
are not counted as executed calculations. Successful six-case sweeps normally
have six in both fields. Process completion itself is not equilibrium evidence.

Before the first engine request, `numerical_qualification.json` contains exact
settings, scripts, hashes and limits. The manifest is updated after every stage.
An engine error, wrong input hash, mismatched source, malformed numeric output,
resource exhaustion or serialization overrun leaves a failed record and raises
`NumericalQualificationError`. No numeric placeholder substitutes for a missing
calculation. Full engine results remain in each stage's `record`; publication
runners can store selected static frames and hashes while retaining a reproduction
recipe for larger raw outputs.

## Result boundary

A completed orchestration has `status="completed_numerical_review"`. Its separate
`numerical_status` is `fixed_geometry_diagnostics_passed` or `requires_review`.
Only a passing review contains `selected_settings`; otherwise that field is null.
An execution failure also clears selection. Every outcome retains
`production_ready=false`, `equilibrium_status="not_assessed"`,
`bulk_property_convergence="not_established"` and
`density_prediction_validation="not_performed"`.

The accompanying tests use explicitly synthetic engine records and hand-set
energy/pressure/force sequences to exercise exact metrics and failures. Their
numbers are not molecular simulation results. Actual numerical evidence must
come from the separately executed reproducibility auditor.
