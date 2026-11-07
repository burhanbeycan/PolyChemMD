# Declared density campaigns and capacity gates

`protocols/campaign.py` implements a strict, dependency-free planning contract.
It launches no engine, creates no run directory, and returns no density estimate.
`DensityCampaignSpec.from_mapping` freezes the declaration;
`plan_density_campaign` reports requested work and missing evidence. This is the
TASK-006C2 planning foundation, not completion of equilibration qualification.

## State declaration

The implemented plan domain is neutral explicit n-alkanes using the frozen
`OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET` model. Each finite molecule has an explicit
carbon count, formula C_n H_(2n+2), and ordinary-isotope masses 12.011 and 1.008
g/mol. The plan retains all chain lengths and molecule counts. It does not
relabel short oligomers as a representative bulk polymer ensemble.

The declared phase is **liquid_candidate**. Temperature, pressure and a requested
phase do not establish the observed phase. A supplied density range is explicitly
a planning hypothesis, not an experimental reference or a fitted output. A null
range leaves the target-capacity gate blocked. Actual phase assessment, suitable
molecular relaxation and an external state-matched reference remain necessary.

Legacy initialization uses `dilute_rigid_then_npt`. The compact declaration also
allows `seeded_rigid_then_npt`, whose name makes no dilute-state assumption.
Both require construction mass/volume density, minimum geometric contact
distance, bounded packing attempts and independent conformer/packing/velocity
seeds. The initial cubic box is computed from exact-SI Avogadro conversion. It is
a requested geometry, not a successful pack, a melt, or proof that bounded NPT
can prepare the target phase. Actual conformer-extent/packing and runtime
box/bond guards can fail, particularly when requesting dense rigid placement.

## A state-capacity limitation exposed before simulation

For an orthorhombic box with each axis strictly larger than twice the largest
requested cutoff, the strict mass-density upper bound is

\[
\rho < \frac{M_\mathrm{cell}}{N_A (2r_c)^3\,10^{-24}},
\]

where `M_cell` is summed atomic molar mass in g/mol and lengths are in angstrom.
This is a geometric upper bound. Passing it is insufficient for packing,
liquid-phase preparation or sampling adequacy.

The example requests four C10H22 molecules, initial construction density
0.06837 g/cm^3, baseline cutoff 10 angstrom, and a paired 11 angstrom cutoff
sensitivity case. Its strict upper bound is approximately **0.0887571 g/cm^3**
when all three settings must be supported. Its declared candidate target range
0.65–0.80 g/cm^3 therefore fails capacity assessment before any engine run.
That range is an illustrative planning hypothesis; this document does not assert
that it is the experimental range at the requested state.

Even sixteen decanes (512 explicit atoms, the preserved legacy/default
periodic-assignment capacity) only allow a strict bound of approximately
0.472543 g/cm^3 at 10 angstrom, or 0.355028 g/cm^3 at 11 angstrom. Thus extending
a dilute small-cell pilot cannot by itself establish a dense-liquid benchmark
at the illustrated range. A larger, separately verified representation is
required. The compact declaration below increases capacity while retaining
numerical-sensitivity requirements; it does not silently reduce the cutoff or
treat short execution as a liquid. The old four-decane example and its
incompatibility remain unchanged.

## Explicit compact capacity and a larger candidate

The top-level optional field selects the assignment contract:

| `assignment_representation` | Total explicit atoms | Per-chain explicit atoms | Additional audit budget |
| --- | --- | --- | --- |
| Omitted, or `explicit_all_pairs_v1` | 512 | 512 | Legacy explicit global-pair representation |
| `compact_sparse_pairs_v2` | 8192 | 512 | At most 2,000,000 summed intramolecular pairs |

The legacy example needs no migration. Omitting the field preserves its exact
v0.6 plan payload and canonical hash; a regression check compares it with the
archived evidence. Unknown selectors are rejected. Compact plans explicitly
report the selected capacity and schema 2 assignment interface. The actual
compact assignment's `nonbonded_rules.representation` is
`sparse_topological_exceptions`; the planner selector is a separate API label.
It retains the same frozen parameter model and topology rules. It does not
select a new potential or change 1–4 scaling, Ewald settings or cutoff.

The intramolecular audit budget is
`sum(n_m * (n_m - 1) / 2)` over molecules, where `n_m = 3*C_m + 2` for the
supported n-alkanes. It limits per-molecule graph/geometry audit work and is
not the number of global atom pairs. For example, sixteen 512-atom chains fit
8192 total atoms but exceed the two-million intramolecular-pair budget and are
rejected. Both global size and per-chain/audit bounds must pass.

[The larger example](../examples/density_campaign_large.json) explicitly selects
compact assignment and requests 72 C10H22 molecules, or **2304 atoms**. Its
construction density input is 0.70 g/cm^3; the corresponding requested cubic
box is approximately **28.96556 angstrom** per axis. It retains a 10 angstrom
baseline cutoff and uses a 12 angstrom cutoff variant, with reciprocal bounds
16 and 20 in the baseline and Ewald variant. The strict all-settings density
upper bound is approximately **1.23058 g/cm^3**. Its illustrative target range
0.65–0.80 g/cm^3 therefore passes this box/cutoff bound. The range remains a
planning hypothesis, and 0.70 g/cm^3 is a construction input rather than a
predicted or observed equilibrium density.

The new file declares three separate cell recipes, each with 72 distinct
conformer seeds and its own packing/velocity seeds. These are independently
seeded **planned rigid-packing requests**, with zero independently prepared
cells executed by the planner. A separately executed bounded tiled-cell audit
is one constructed cell; tiling that cell or translating it does not realize
these three independent preparations. No tiling evidence is silently counted
as an independent replicate or as successful dense random packing.

The compact plan remains `execution_supported: false`: there is no automatic
campaign dispatcher, validated dense initialization, phase certificate, or
finite-size production campaign. The geometric pass removes one size barrier;
it does not close state qualification, sampling or external comparison.

## Cells, segments and numerical variants

At least three independently constructed cell recipes are declared. Each chain
in each primary cell gets a distinct conformer seed; packing and velocity seeds
are also distinct across the primary recipes. Different seeds are a reproducible
construction recipe, not a proof of statistical independence. Realized source
hashes and independent preparation histories still need to be audited.

Three numerical settings are required:

- The baseline cutoff, Ewald splitting parameter and reciprocal bound.
- A larger cutoff with the other two settings fixed.
- A larger reciprocal bound with cutoff and splitting parameter fixed.

These comparisons change one numerical control at a time. Their measured
sensitivity can guide further refinement, but two parameter choices alone do not
prove convergence of the infinite Ewald sum or LJ truncation. The numerical
variants reuse each primary cell's recipe for paired comparisons. Three cells
with three numerical settings means three primary cell recipes and nine paired
work requests, not nine independent samples. Restart segments and trajectory
frames are correlated observations of their parent cell and are not new cells.

Both examples request 5000 NVT plus 5000 NPT preparation steps at 0.25 fs, then
40000 production steps split into four 10000-step continuations. The requested
production duration is 10 ps per cell and numerical setting, with 2.5 ps of
preparation. These are declared durations, not an assertion that 10 ps suffices.
The request contains 45 engine runs and 450000 total integration steps. None of
this work is launched by the planner. The legacy example is capacity-incompatible;
the compact example passes only the declared box/cutoff size bound.

## Resource accounting

The declaration specifies hardware identity, threads, maximum aggregate wall
budget, maximum run count and per-run timeout. The plan counts requested engine
runs and steps separately from minimizer iteration/evaluation budgets. An exceeded
run-count budget is retained explicitly. The sum of timeouts is a worst-case
allowance if all calls independently time out, not a measured runtime or a total
budget guarantee. There is currently no campaign dispatcher enforcing an
aggregate wall clock, so `execution_supported` remains false.

Optional calibration fields are:

```json
{
  "hardware_label": "host_identifier",
  "engine_executable_sha256": "64 lowercase hexadecimal characters",
  "atom_count": 128,
  "integration_steps": 1000,
  "wall_seconds": 2.0,
  "cutoff_angstrom": 10.0,
  "reciprocal_nmax": 16
}
```

The numbers above illustrate the schema, not an actual calibration. Both example
files use null calibration. Calibration atom count is checked against the
selected representation's capacity and must match the actual planned cell
(128 atoms in the legacy example or 2304 in the compact example). The planner requires matching declared hardware,
atom count, cutoff and reciprocal bound, and a supplied executable hash. It
cannot authenticate a user-supplied measurement. It returns only a labeled
linear scheduling projection for baseline integration. Startup, checkpoint I/O
changes, minimization, different numerical settings and finite-size comparisons
are excluded. This projection cannot verify the full campaign budget and must
not be reported as measured future work. All dimensional and derived values
must remain finite and positive; underflow/overflow gives an explicit error.

## Qualification and external comparison

The declaration retains preselected duration, density drift and relative
interval-width limits, temperature/pressure mean deviations, numerical-sensitivity
limits and minimum independent-cell count. These are review thresholds for this
study, not universal validated convergence criteria. The planner does not
calculate or pass any of them. Their sampling methods and relevant correlation
and state assumptions need actual evidence.

The following gates remain `not_evaluated`: realized cell/topology/parameter
identity; actual same-state NPT sampling and restart lineage; temperature and
pressure sampling; observed phase; density stationarity, correlation-aware
uncertainty and block-size sensitivity; agreement across independently prepared
cells; cutoff/Ewald sensitivity; finite-size comparison; and recorded resource
use with retained failures. Finite-size assessment cannot be disabled in this
contract and is not implemented by the planner. Passing a duration threshold or
obtaining stable density does not establish chain relaxation or the convergence
of a transport property.

External comparison is deferred until these state and sampling gates pass. Its
reference must match chemistry, chain-length distribution, temperature, pressure,
phase and density definition. A failed/inconclusive state remains visible;
planning never upgrades it to an equilibrium density or validated prediction.
