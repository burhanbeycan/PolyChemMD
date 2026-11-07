# Initial bounded development tasks

These are repository task specifications, not claims that GitHub issues have already been opened. Keep one task per reviewed change where feasible. Read AGENTS.md and the scientific contracts before implementation.

## Milestone status — 7 November 2026

| Task | Current state | Evidence boundary |
| --- | --- | --- |
| 001 | Implemented for mapped linear repeat units | Explicit graph/attachment contract and atomic lineage; broad reaction chemistry remains future work |
| 002 | Implemented for supported single-bond head-to-tail graphs | H/methyl caps, graph tests, explicit unsupported stereochemical cases; no coordinates/cells |
| 003 | Implemented | Analytic energy checks for LAMMPS Fourier/OPLS conversions, not material parameter validation |
| 004 | Implemented contract | Assigned-type coverage only; no chemical atom typing or production parameter set |
| 005 | Implemented process adapter | Real reduced-unit LJ smoke, plus independent failure/parser tests |
| 006 | Geometry, restricted parameters, TASK-006C1 plus TASK-006C2A/B/C1 restart, compact cells and qualification-analysis foundations implemented | No validated polymer-density result exists; sampling and external comparison remain open |
| 007 | Planned | No representation advantage or chemical-family benchmark claimed |
| 008 | Development evidence audit implemented for v0.2–v0.8 | Machine-readable geometry/parameter/periodic/engine checks; complete scientific benchmarks remain open |

## TASK-001 — Chemical contract and atom-lineage schema

**Objective:** specify and implement attachment roles, atom mapping, repeat-unit identity and termini without changing existing validated behavior.

**Inputs:** the current material-record schema; independently designed examples for a linear repeat unit, two-unit linkage and capped chain.

**Acceptance:** reject duplicate maps, incompatible attachments and missing unit/state information; preserve atom lineage through round trips; distinguish repeat-unit assembly from monomer transformation; document migration if the schema changes.

**Excluded:** coordinate generation, force-field assignment, polymerization kinetics and MD property output.

## TASK-002 — Independent mapped repeat-unit linear-chain assembly

**Prerequisite:** TASK-001 contract is stable.

**Objective:** implement a finite linear chain builder with explicit bond and cap rules.

**Acceptance:** hand-checked one-, two- and several-unit graphs; exact atom-count balance; deterministic output; no duplicate atom identities; valence/formal-charge and stereochemical checks; clear errors for unsupported patterns. Optional RDKit validation must be isolated from the lightweight core.

**Excluded:** treating a free monomer as an already transformed repeat unit, generating a physical reaction trajectory, or reporting density from a graph.

## TASK-003 — Torsion import and descriptor audit

**Objective:** audit the existing generic Fourier routine, add an explicitly named external-convention importer and test conversion with independent analytical values.

**Acceptance:** every term, phase, periodicity, prefactor and unit is preserved; numerical grids are labeled; multi-term examples distinguish profiles with the same largest coefficient; round-trip and convention tests use independently calculated references.

**Excluded:** claiming that generic toy coefficients are validated parameters for a named polymer.

## TASK-004 — Parameter-set contract and coverage report

**Objective:** represent functional forms, parameter provenance, exclusions, mixing rules and 1–4 scaling. Produce a coverage report from actual assigned records.

**Acceptance:** missing parameters stop production preparation; compatibility failures are explicit; parameter assumptions are logged; no guessed or silently substituted production parameters.

## TASK-005 — First external-engine adapter

**Objective:** implement engine discovery, planned-input creation, a bounded run and parsing of current-run artifacts.

**Acceptance:** missing binary, timeout, nonzero exit, nonconvergence and stale files are tested independently; dry runs cannot generate computed observables; actual integration tests record executable identity and source input hashes.

**Excluded:** fabricating a successful physical result when an executable is unavailable.

## TASK-006 — First end-to-end density protocol

**Prerequisites:** supported chemistry, chain/cell construction, parameter coverage and a working engine adapter.

**Objective:** run and analyze a documented amorphous-polymer density example, with independently prepared repeats and sampling diagnostics.

**Acceptance:** actual input and trajectory artifacts, physical units, thermodynamic state, convergence status, uncertainty and a suitable external comparison. Dataset count must distinguish independently prepared cells from correlated frames.

The analysis kernel and the bounded TASK-006A preparation layer are implemented.
The full milestone remains open:

- **TASK-006A:** independent conformer/coordinate generation and cell packing,
  with topology/lineage preservation, overlap checks and deterministic seeds.
  **Implemented in v0.3:** optional RDKit ETKDGv3 conformers with coordinate-derived
  specified-stereo checks; independent rigid placement in orthorhombic cells;
  minimum-image intermolecular checks, conservative self-image exclusion and
  wrapped/unwrapped coordinates with image flags. Molecules must fit the
  conservative extent bound. No melt packing, relaxation, equilibrium sampling,
  conformational search or physical density prediction is implied.
- **TASK-006B:** chemically assigned production parameters and charges from
  attributed sources or independently fitted references; demonstrate convention,
  energy/force and coverage checks. The analytic parameter-contract demo is not
  a production force field and must not be used as one.
  **Implemented in v0.4 for neutral explicit n-alkanes:** exact sourced OPLS-AA
  CT/HC numerical subset, CH3/CH2/H typing, fixed charges, all proper torsions,
  exclusions/1–4 pairs, strict frozen-assignment identity, analytic energy/force
  reference and actual isolated LAMMPS single-point parity. This is parameter
  assignment and convention/implementation verification, not physical accuracy
  validation or a complete periodic production protocol.
- **TASK-006C:** actual polymer equilibration/production campaign, independent
  cells, observable-specific diagnostics and phase/state-matched comparison.

A reduced-unit LJ smoke test establishes engine execution, not any of these
polymer-specific prerequisites. Density diagnostics accept only explicitly
converted g/cm^3 and do not certify equilibrium.

### TASK-006C1 — implemented within the restricted periodic model

Assemble the supported per-molecule assignments into an actual periodic cell
without losing global/source identities, images or pair conventions. Define and
test periodic bonded images, nonbonded cutoff/long-range treatment, charge
neutrality and inter-molecular terms against an independent reference or engine
parity case. Record the selected bulk model rather than silently reusing the
isolated truncated verification setup. Add an explicit relaxation/equilibration
protocol with finite resources and honest failure/inconclusive states. A long
density campaign needs a declared hardware/time budget and sampling gates;
short energy/parity checks cannot close TASK-006C. Transferability beyond the
supported n-alkane domain remains separate work.

**Implemented in v0.5:** strict global/source identity and image-preserving cell
assembly; complete intermolecular pairs; explicit 3D Ewald model and independent
analytical force reference; actual periodic engine parity including a translated
boundary-crossing representation; finite minimization/NVT/NPT protocol with
resource/box guards and final-state/unit audits. Both a source-invariant geometry
record and a subsequent deformed engine state have separate roles. Native and
canonical density values remain explicit. Short preparation does not close the
full density milestone.

### TASK-006C2 — state and sampling qualification (still open)

Qualify the preparation and sampling protocol for a precisely declared state and
chemistry before reporting a density prediction. Define target phase, chain/cell
sizes, initial-density strategy, long-range/cutoff sensitivity, replicate seeds,
physical duration and measured hardware/time budget. Add restartable production
states and observable-specific diagnostics; distinguish independent cells from
correlated frames. Retain failed/inconclusive cases and perform a state-matched
external comparison only after the sampling gates pass. The dilute four-oligomer
pilot and its 200 steps do not establish these conditions.

### TASK-006C2A — restart, diagnostics and campaign foundations

Implemented in v0.6 within the restricted periodic n-alkane model: native NPT
checkpoint manifests and bounded continuation, same-cell split/restart auditing,
explicit physical-time density diagnostics, independent-cell summary contracts,
and a state/resource-aware no-execution campaign planner. A diagnostic pass is
not an equilibrium certificate. Actual short runs remain insufficiently sampled;
no external density comparison or validated polymer density is reported.

### TASK-006C2B — admissible larger periodic cells

Remove the present periodic assignment's 512-atom/all-pair capacity bottleneck
with a separately tested representation while retaining frozen parameter,
identity, exclusion and image invariants. Verify that declared target density,
cell mass and cutoff permit a three-axis minimum-image cell. For the current
four-decane fixture, a 10 angstrom cutoff requires a volume strictly above
8000 angstrom cubed; the pilot cannot represent the planned liquid density.
Do not shrink the cutoff just to fit. Demonstrate a larger bounded cell and
numerical cutoff/long-range sensitivity before committing to a production
campaign. Measure cost on available hardware; long sampling, finite-size checks,
phase/chain relaxation and state-matched external evaluation remain later gates.

**Implemented in v0.7 within the restricted model:** compact sparse pair rules
with exact counts and strict frozen-source reconstruction; bounded ordered motif
construction and periodic linked-bin contacts; engine export/preparation/restart
support for larger cells; fixed-geometry numerical sensitivity evidence. The
legacy model and serialized schema remain unchanged. The larger geometric cell
removes the former capacity obstacle for the declared density range, but does
not identify a liquid or certify sampling. Actual outcomes are in VALIDATION.md.

### TASK-006C2C — state and numerical preparation qualification (open)

Establish an appropriate preparation history for the declared material/state
instead of treating the ordered repeated motif as a liquid or independent
replicate. Retain minimizer stopping/force evidence, assess alternative starting
structures and independent cells, and qualify cutoff/long-range settings with
observable-specific evidence. Select and justify thermalization/discard and
physical sampling windows under a measured finite compute budget. Continue to
retain failed/inconclusive states. Phase/chain relaxation, finite-size sampling
and state-matched external density comparison remain open until actual gates
pass; a short dense geometric pilot cannot close them.

### TASK-006C2C1 — numerical and preparation review layer

Implemented in v0.8: explicit static virial scalar/tensor pressure, a bounded
multi-point cutoff/reciprocal protocol with consecutive finest-refinement gates,
static Rg/end-to-end/gyration/nematic calculations, and a full-record preparation
readiness review. Source/geometry/engine identity and all numerical settings are
prespecified; initial-state numerical coverage does not transfer through MD.
Short ordered-cell execution evidence is distinct from state qualification.
Actual outcomes and limitations are recorded in VALIDATION.md.

### Next bounded task — TASK-006C2C2: preparation history and evolved-state coverage

Resolve the retained pressure-sensitive cutoff result with declared numerical
refinements and appropriate cell size before adopting a production model.
Establish suitable, independently prepared starting structures and an affordable
relaxation/thermalization protocol. Retain actual minimizer stopping evidence
without requiring finite-temperature forces to vanish. Check whether imposed
motif orientation persists; do not relabel perturbations of one parent cell as
independent material replicates. Extend numerical sensitivity to actual evolved
states or representative trajectory windows before adopting it for production.
Prespecify physical duration/discard and observable sampling gates, measure
resources, and retain failed/unqualified outcomes. Full phase/chain relaxation,
finite-size and state-matched external density evaluation remain open.

## TASK-007 — Chemically conditioned torsion descriptors

**Objective:** partition local environments by backbone, sidechain, termini and junctions; retain full torsion landscapes and explicitly defined sequence correlations.

**Acceptance:** matched-composition examples test whether the representation distinguishes topology/sequence changes; compare with simple chemistry and marginal FF baselines; record feature computation cost and avoid target leakage.

## TASK-008 — Scientific release/evidence audit

**Objective:** generate a human-readable matrix from implemented capabilities, tests, execution records and external comparisons.

**Acceptance:** each claimed feature links to actual code and evidence; planned functionality remains planned; the manuscript and README describe the same release state; no invented benchmark values or citations.

## Change report template

- Objective and scientific assumption.
- Files and public contracts changed.
- Independent verification performed, with commands and outcomes.
- External engine runs actually executed, or an explicit statement that none ran.
- Remaining limitations and next bounded task.

Test results must come from executed commands. A task checklist is not a test result.
