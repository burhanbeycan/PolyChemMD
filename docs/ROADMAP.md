# Development roadmap

This is a staged research-software plan, not a statement that the complete platform exists. Version 0.4 adds a sourced n-alkane parameter/charge assignment and independent analytic-versus-LAMMPS single-point parity to the graph/geometry/packing/analysis core. Broader chemical coverage, periodic parameter assembly, bulk protocol and material-property validation remain open. Current task evidence is recorded in docs/TASKS.md and VALIDATION.md.

## P0 — Freeze the chemical and scientific contracts

Deliver a versioned schema for repeat units, attachment-site roles, termini, microstructure, finite realizations and material state. Retain an explicit distinction between measured, assumed and computed information.

Acceptance:

- Invalid or ambiguous attachment definitions are rejected with actionable errors.
- Formal-charge and valence requirements are specified independently of any force field.
- Atom-lineage and mapping rules survive JSON serialization.
- Distribution targets and generated statistics have separate fields.
- A deterministic specification fingerprint has tests for scientific-input changes; volatile timestamps are not scientific inputs.
- Dependency-free tests remain runnable without external solvers.

## P1 — Independent mapped linear-chain builder

Implement repeat-unit connectivity assembly first. A repeat unit with attachment markers is not a free monomer; monomer-to-repeat-unit reactions require a separate transformation module.

Acceptance:

- Deterministic one-, two- and multi-unit examples have hand-checked connectivity and atom counts.
- Atom mappings and chemical lineage are preserved without duplicate identifiers.
- End-group transformations and any atom loss are explicit.
- Bond orders, valence, stereochemistry and formal charge are checked.
- The builder rejects unsupported attachment patterns rather than guessing.
- Optional RDKit-based sanitization is tested when installed; its absence does not break core imports.
- No coordinates, density or physical-property accuracy is claimed from connectivity assembly alone.

## P2 — Parameter and execution contracts

Define complete functional forms, units, exclusions, scaling factors and source identities. Build engine discovery, input generation, execution and artifact parsing as separate interfaces.

Acceptance:

- A test fixture checks energy/force and torsion conventions independently.
- Unknown or unmapped chemistry stops production execution.
- No automatic silent mixing of incompatible parameter conventions.
- Missing executable, timeout, nonconvergence and stale output are distinct failures.
- A dry run produces only planned inputs and never reports a simulation result.
- Actual engine tests are separate from lightweight tests and record executable versions.

## P3 — First validated polymer workflow

Complete one supported chemistry/state/protocol from specification to actual MD artifacts and one basic observable, then expand to independent chemical families. Density is an appropriate early target only after chemical and parameter checks pass.

TASK-006A supplies the starting geometry layer. It deliberately does not assign
a thermodynamic phase or optimize the requested cell density. TASK-006B now
supplies sourced typing/charges and energy/force parity for restricted neutral
n-alkanes. That static finite-cutoff check is not a periodic production polymer
model; TASK-006C1 must assemble the periodic parameterized system and define its
bulk treatment before relaxation, sampling and an external density comparison.

Acceptance:

- Fresh preparation and independent random seeds are retained.
- Stationarity and autocorrelation-aware uncertainty are reported.
- Chain length, preparation and system-size sensitivity are assessed for the claim.
- Reference comparisons use matching temperature, phase and material information.
- Failed and inconclusive runs remain visible in the execution record.

## P4 — Polymer-specific descriptors and mechanistic analyses

Add full torsion-landscape, linkage-conditioned, sequence and conformational descriptors. Add persistence/Kuhn analysis, hydrogen-bond networks, structure factors and probe-defined free volume.

Acceptance:

- Definitions and units are explicit, including backbone selection and hydrogen-bond criteria.
- Descriptor ordering and schema versions are stable.
- Baselines and ablations assess which information supplies a benefit.
- Pre-simulation and post-simulation features are separated.
- Molecular-weight, sequence, tacticity and unseen-chemistry tests match the proposed claim.

## P5 — Property families

Develop each family through its own protocol and validation gate; a registry entry remains planned until that evidence exists.

| Family | Main additions | Gate |
| --- | --- | --- |
| Mechanical | Elastic tensors, stress–strain and strain-rate response | Finite-temperature and linear-response validation |
| Dynamic | Relaxation, creep, viscosity and oscillatory response | Sampling and terminal/linear-regime evidence |
| Transport | Species diffusion, sorption, permeability and ionic correlations | Correct loading, state and transport estimator |
| Interfacial | Adsorption PMF and interfacial thermal/structural response | Surface model and enhanced-sampling convergence |
| Electronic/spectroscopic | Ensemble response, IR/Raman/UV–vis and selected reaction references | Appropriate QM/polarization model and reference conditions |

## P6 — Research release

Prepare a tagged reproducible release with installation instructions, small verification fixtures, source and parameter attribution, an execution manifest, benchmark recipes, and an evidence matrix linking manuscript claims to actual outputs. Freeze a release only after its relevant checks pass. Do not manufacture a completed benchmark or submission-ready results from planned capabilities.

## Execution environments

Use a lightweight environment for routine development/CI. Use a separately recorded scientific environment for RDKit, LAMMPS and QM engines. GPU availability, CUDA support and engine acceleration must be measured on the selected host; the coding service is not assumed to provide a suitable scientific GPU.

Milestones are gated by evidence rather than a fixed calendar estimate. Production cost should be estimated after representative engine timing and convergence pilots.
