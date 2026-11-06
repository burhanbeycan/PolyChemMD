# Initial bounded development tasks

These are repository task specifications, not claims that GitHub issues have already been opened. Keep one task per reviewed change where feasible. Read AGENTS.md and the scientific contracts before implementation.

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
