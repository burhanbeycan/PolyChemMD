# Bounded periodic alkane preparation pilot

This protocol connects an assigned periodic n-alkane cell to actual LAMMPS
minimization, NVT and NPT stages. It is a finite preparation workflow, not an
automatic equilibrator or a validated polymer-density calculation.

```python
from polychemmd.protocols import (
    PeriodicAlkaneProtocolSpec, plan_periodic_alkane_protocol,
    run_periodic_alkane_protocol,
)

spec = PeriodicAlkaneProtocolSpec(
    temperature_k=300.0, pressure_bar=1.0, timestep_fs=0.25,
    minimize_max_iterations=200, minimize_max_evaluations=1000,
    nvt_steps=100, npt_steps=100, velocity_seed=49217,
    timeout_seconds=60.0, cutoff_angstrom=10.0,
    ewald_alpha_per_angstrom=0.30, reciprocal_nmax=16,
)
plan = plan_periodic_alkane_protocol(
    assigned, spec, run_directory="artifacts/new-protocol-run")
result = run_periodic_alkane_protocol(
    assigned, spec, run_directory="artifacts/new-protocol-run")
```

`plan` discovers/runs no engine and creates no files. It retains the exact setup,
protocol commands, staged data, hashes, model identity, seed, requested resources
and unit conversions. Execution requires an installed working engine and a new
run directory. Reuse of an existing directory is rejected; an unavailable or
failed engine never receives substituted numerical outputs.

## Scientific model and stage controls

The scientific setup comes directly from
`engines.periodic_alkane.build_periodic_alkane_setup`: the frozen attributed
OPLS-AA alkane parameters, full explicit-atom topology, fixed charges,
orthorhombic periodic boundaries, geometric LJ mixing, 0/0/0.5 special-pair
scaling and the declared periodic Ewald model. The protocol does not slice or
silently replace an isolated-model input script. Cutoff, Ewald alpha and
reciprocal truncation are retained in its specification.

1. CG minimization has a maximum per-iteration atom move of 0.1 angstrom.
   Its energy tolerance is **1e-6, dimensionless**. Its native LAMMPS real-unit
   force tolerance is **1e-6 kcal/(mol angstrom)**, equivalent to
   **4.184e-6 kJ/(mol angstrom)**. Iteration and evaluation budgets are separate.
2. Velocities are drawn with the explicit integer seed, Gaussian distribution,
   and linear/angular momentum removal. The MD timestep is in femtoseconds.
3. NVT uses the target kelvin temperature and a fixed 100 fs thermostat damping.
4. NPT uses the same temperature control and isotropic target pressure, with a
   fixed 1000 fs barostat damping. Canonical bar input is converted to LAMMPS
   real-unit atmospheres by **P_atm = P_bar / 1.01325**. Thermo pressures are
   converted back to bar in the public protocol record.

These are transparent preparation choices, not optimized or physically
validated sampling parameters. The defaults integrate only **25 fs NVT plus
25 fs NPT**. Minimizer iterations are not MD timesteps and are not added to
that physical trajectory duration.

The minimizer's actual stopping criterion and reported iteration/evaluation
counts are parsed separately from process success. Only an explicit engine
`force tolerance` stop is recorded as meeting that criterion. An energy-only
stop retains force-unverified status. A maximum-iteration/evaluation stop is
`budget_exhausted`; other or missing stop information is `inconclusive`.
The fixed policy permits the short MD pilot after a budget-exhausted minimizer,
while retaining that limitation. It never turns a completed process into a
minimizer-convergence certificate.

## Periodic cutoff and actual-output gates

The initial setup requires cutoff less than half the shortest box length. NPT
can change that box. A `fix halt` checks explicit per-axis Boolean conditions
on every minimizer and MD step; any axis reaching `2*cutoff` is a **hard engine error**. Reported thermo
boxes and the actual final box are independently checked again. Thus the pilot
cannot silently claim successful execution outside its prescribed unique
minimum-image real-space regime.

A second every-step hard guard requires the maximum bond length to remain
below `cutoff/3`. By the triangle inequality, every supported 1–4 endpoint
separation is then less than cutoff and half the minimum box length, preserving
the selected bonded-image/exclusion regime. The final unwrapped bond lengths
are checked independently. This sufficient condition is intentionally
conservative; it is not a bond-breaking or reactive-model criterion.

The halt uses an equal-style variable. The generic adapter correctly retains
its conservative recipe-unit inference as `unresolved`, while the run
specification separately declares `real`. The protocol independently inspects
the actual `write_data` header and requires `units = real` before accepting its
dimensional output. A caller declaration is not used to overwrite ambiguous
recipe inference.

Actual minimization/NVT/NPT markers and three complete thermo blocks are
required. NVT and NPT first/last steps must match their requested counts exactly.
Missing completion, premature stopping, malformed/nonfinite output or wrong
sampled/final cutoff regime yields `PeriodicProtocolError` with a failed record.

The final single-frame dump contains fixed atom/molecule IDs and charges,
wrapped positions, image integers and forces. Atom loss, duplicate IDs, changed
charge/molecule identity, wrong step, nonperiodic/triclinic bounds or positions
outside the half-open box are rejected. Unwrapped positions are reconstructed
as `wrapped + image*L`, preserving the actual possibly shifted NPT box origin.
The final mass/volume density is checked against the actual thermo density using
the explicitly pinned engine convention. LAMMPS real units use a molar-mass
conversion factor `1/0.602214129` in the tested source, whereas the exact SI
Avogadro constant gives `1/0.602214076`. The small numerical difference is not
hidden by relaxing a tolerance: `density_g_cm3` uses exact-SI mass/volume,
`engine_density_g_cm3` preserves actual thermo output, and the native-factor
audit/source are retained separately. These are observed instantaneous pilot
states, not property predictions.

## Resources, state and reproducibility

The 12 specification fields serialize exactly and can be loaded through
`PeriodicAlkaneProtocolSpec.from_mapping`. Booleans and nonfinite numerical
inputs are rejected. Limits are explicit: at most 10000 total integration
steps, 5000 minimizer iterations, 50000 evaluations, 120 s engine timeout,
1 fs unconstrained explicit-H timestep and reciprocal truncation index 20.
These API limits are separate from the smaller default pilot budget; raising
a budget does not establish equilibrium or force-field applicability.

Successful output remains `completed_preparation_pilot`, with
`equilibrium_status: insufficient_sampling`, `production_ready: false` and
`density_prediction_validation: not_performed`. It includes exact stage counts,
source-assignment and input hashes, scientific controls, actual engine identity
and artifacts. `final_state` is distinct from the initial rigid source packing:
real minimization/MD positions are not misrepresented as a new rigid packing.
It does not establish amorphous phase, chain relaxation, density convergence,
transport, uncertainty or experimental agreement. A later observable-specific
campaign needs independent preparations and justified sampling gates.

The expected scientific artifacts are `input.lammps`, `system.data`,
`log.lammps`, `prepared.data` and `final.dump`; process and protocol records
remain separate. Tests with synthetic output frames check these contracts only;
actual engine integration evidence must be reported separately.

## Primary method references

- [LAMMPS NVT/NPT and damping conventions](https://docs.lammps.org/fix_nh.html)
- [Minimizer tolerances and stopping conditions](https://docs.lammps.org/minimize.html)
- [Runtime halt conditions for dynamics/minimization](https://docs.lammps.org/fix_halt.html)
- [LAMMPS real units](https://docs.lammps.org/units.html)
- [Seeded velocity initialization](https://docs.lammps.org/velocity.html)
- [Pinned seed range in RanMars](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/random_mars.cpp#L38)
- [Pinned real-unit mass-density conversion](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L175)

These references support controls and numerical conventions, not the claim that
the short pilot equilibrates a polymer material.
