# PolyChemMD

Independent polymer chemistry simulation and property-analysis research project.

**Release state: 0.7.0.dev0 — compact periodic cells, ordered molecular construction and larger-cell engine sensitivity audits; equilibrium polymer-property validation remains open.**

The long-term scope is a general-purpose connection from polymer chemistry and
microstructure to quantum-informed parameterization, molecular dynamics,
descriptors and material properties. Application examples do not narrow the
platform to a particular filler, solvent, polymer family or device.

Start with [the Turkish Codex guide](CODEX_START_HERE_TR.md), then read
[AGENTS.md](AGENTS.md), [architecture](docs/ARCHITECTURE.md),
[scientific contracts](docs/SCIENTIFIC_SPEC.md) and [bounded tasks](docs/TASKS.md).

## Current implementation and evidence boundaries

| Capability | Status and boundary |
| --- | --- |
| Material metadata records | Implemented: schema, fraction, microstructure label and state validation; no molecular graph/valence validation |
| Provenance envelopes | Implemented: software/Python versions, UTC time and exact input-file SHA-256 |
| Repeat-unit chemical contract | Implemented for explicitly mapped, already transformed units; attachment and source-atom lineage checks |
| Finite linear-chain graph builder | Single-bond head-to-tail assembly with H/methyl caps; unsupported stereo cases fail explicitly; no coordinates/cell preparation |
| Molecular coordinate contract | Immutable atom/topology/coordinate records, explicit angstrom units, lineage and strict JSON round trips |
| Seeded 3D conformer generation | Optional RDKit ETKDGv3 backend, finite budgets, atom/graph preservation and coordinate-derived specified-stereo checks; no conformational search or energy ranking |
| Periodic cell preparation | Independent rigid-molecule rejection packing in orthorhombic cells, proper rotations, periodic contacts, conservative self-image bounds and image flags |
| Structure visualization export | Extended XYZ with atom order, lattice and periodicity; topology, isotopes and provenance remain in companion JSON |
| Complete Fourier torsion kernels | Implemented and analytically tested: all coefficients, phases and multiplicities retained |
| Torsion descriptors | Grid extrema/barrier/minima and isolated-model Boltzmann bin probabilities; illustrative coefficients only |
| LAMMPS torsion imports | Explicit real-unit Fourier and OPLS conversions with independent equation tests |
| Parameter contracts and coverage | Units, forms, provenance, exclusions and 1-4 scaling; assigned-type completeness, not chemical atom typing or force-field accuracy |
| Restricted chemical atom typing | Neutral explicit n-alkanes: CH3/CH2/HC types, complete bonds/angles/all proper four-atom paths and shortest-path pair scaling; unsupported chemistry fails |
| Sourced fixed-charge parameter assignment | Frozen GROMACS v2025.3 OPLS-AA CT/HC numerical subset with per-type source lines and hashes; empirical charges, no new QM charge fitting |
| Independent energy/force reference | Original local analytic differentiation of harmonic, Fourier, LJ and Coulomb equations; isolated, unshifted finite-cutoff model |
| LAMMPS energy/force verification | Generated full-atom input, run-zero single points, strict identity/parser checks and component/force comparison; no relaxation or property validation |
| Periodic parameter assembly | Legacy exhaustive schema 1 and compact schema 2 preserve global/source identities, images and frozen parameters; compact storage uses unit-scale defaults and short-topology exceptions |
| Larger periodic starting cells | Compact assignment supports at most 8192 atoms, with per-molecule/local-work limits; actual engine evidence uses 2304 atoms |
| Ordered motif construction | Independent deterministic molecular tiling with periodic linked-bin contacts and rigid/source checks; constructed density and phase claims remain separate |
| Periodic Ewald reference and engine parity | Explicit real/reciprocal/self/exclusion terms, analytic forces, minimum-image gates and declared finite-sum budgets; no automatic convergence claim |
| Bounded periodic preparation | Reviewable minimization, NVT and NPT pilot with box/resource guards, native/canonical units and final-state checks; short execution never certifies equilibrium |
| LAMMPS subprocess adapter | Bounded execution, exclusive run directory, output hashes and thermo checks; real reduced-unit LJ smoke executed |
| Native NPT checkpoint/continuation | Hash-bound restart with parameter/state/engine identity checks and retained thermostat/barostat state; restricted to the same supported engine build and platform |
| Density analysis | Physical-unit conversion, autocorrelation/block diagnostics and independent-cell summaries; diagnostics never certify equilibrium or external accuracy |
| Density campaign planning | Declared state, seed identities, physical durations, numerical sensitivities and resource/cell-capacity gates; no automatic campaign execution |
| Dependency inspection and CLI | Implemented; presence inspection does not execute or validate engines |
| Property catalogue | 19 planned material-property families; the catalogue is not an executable property engine |
| Broader chemistry, production sampling and state-matched bulk benchmarks | Planned |
| Validated polymer-property MD and QM protocols | Planned |
| Experimental validation and publication benchmarks | Planned; no performance or material-property result is claimed |

RadonPy and other polymer workflow packages are neither dependencies nor vendored
source. The lightweight runtime core uses only the Python standard library.
Optional RDKit chemistry and external solvers are clearly separated from this
project's independently written code and scientific contributions.

## Install and check

Python 3.11 or later, from the repository root:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
polychemmd doctor
polychemmd validate examples/pe_record.json
polychemmd catalogue
polychemmd torsion-describe examples/torsion_demo.json --n-grid 720 --temperature-k 300 --output artifacts/torsion_demo.json
```

Install optional chemistry support and assemble a molecular graph:

```bash
python -m pip install -e '.[chemistry]'
polychemmd build-linear examples/linear_pe.json --output artifacts/linear-pe.json
polychemmd torsion-import examples/torsion_lammps_fourier.json --source lammps-fourier
polychemmd parameter-check examples/parameter_contract_demo.json --require-complete
polychemmd build-conformer examples/conformer_pe.json --output artifacts/pe-conformer.json
polychemmd pack-cell examples/packed_pe.json --output artifacts/pe-cell.json
polychemmd export-xyz artifacts/pe-cell.json --output artifacts/pe-cell.xyz
```

The geometry examples prepare four C10H22 molecules in a 24 angstrom cube,
using explicit conformer and packing seeds and a 1.5 angstrom all-atom geometric
cutoff. The cutoff is not a van der Waals model. The reported construction
density is mass divided by the requested volume, not a density prediction.
Random rigid placement does not establish an equilibrated amorphous phase.
See [coordinate generation](docs/COORDINATE_GENERATION.md) and
[cell packing](docs/CELL_PACKING.md) for supported geometry and numerical limits.
XYZ is a viewing format; preserve the JSON for identity, isotope and topology
information. Default cell XYZ coordinates are unwrapped; `--wrapped` is available.

Execute the bounded examples and an independent lattice-image/rigidity audit:

```bash
python scripts/verify_geometry_v03.py --output-dir artifacts/geometry-v03
```

The output directory must be new. The script also generates finite PE, PVDF and
PS repeat-graph conformers; their tacticity/phase is not inferred. Actual example
records and their limitations are in [geometry evidence](evidence/geometry_v03/README.md).

Assign the supported force field and inspect an isolated single point:

```bash
polychemmd parameterize-alkane artifacts/pe-conformer.json --output artifacts/assigned-pe.json
polychemmd alkane-energy artifacts/pe-conformer.json --output artifacts/pe-energy.json
polychemmd alkane-lammps artifacts/pe-conformer.json --run-dir artifacts/pe-singlepoint-plan --plan --output artifacts/pe-plan.json
```

The plan writes no engine directory and runs nothing. Omit `--plan` to execute
the bounded run-zero comparison with an installed LAMMPS. The actual three-case
verification recipe (ethane, butane, decane) is:

```bash
python -m pip install -e '.[chemistry,md-smoke]'
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_alkane_v04.py --output-dir artifacts/alkane-v04
```

Use a new output directory. See [assignment](docs/ALKANE_ASSIGNMENT.md),
[reference equations](docs/ENERGY_FORCE_REFERENCE.md) and
[LAMMPS parity](docs/ALKANE_LAMMPS.md). This subset excludes branches, rings,
ions, heteroatoms and labelled isotopes. Its model is tied to an exact upstream
parameter snapshot; it is not an assertion that every OPLS-AA release is identical.
The finite-cutoff, nonperiodic verification model is separate from a bulk MD
protocol with independently chosen long-range treatment and sampling.

Assemble the supported periodic cell and inspect or run its explicit Ewald setup:

```bash
polychemmd parameterize-cell artifacts/pe-cell.json --output artifacts/assigned-cell.json
polychemmd periodic-lammps artifacts/assigned-cell.json --run-dir artifacts/periodic-plan --plan --output artifacts/periodic-plan.json
polychemmd periodic-prepare artifacts/assigned-cell.json --protocol examples/periodic_protocol.json --run-dir artifacts/preparation-plan --plan --output artifacts/preparation-plan.json
```

Both plans run nothing and create no engine directory. Omit `--plan` to execute
the separately bounded calculation in a new run directory. The first command
accepts the original packed cell or its result envelope; later commands also
accept the strictly validated assigned-cell record. `periodic-energy` evaluates
the independent Ewald reference without an engine.

The reproducibility runner checks the original 128-atom cell, a boundary-crossing
translation of the same state, and a short minimization/NVT/NPT preparation:

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_periodic_v05.py --output-dir artifacts/periodic-v05
```

Its translated cells are an invariance test, not independent physical samples.
The selected model is 3D periodic Ewald with explicit alpha/reciprocal bounds,
real-space cutoff, geometric LJ mixing, and no LJ tail correction. See
[periodic assignment](docs/PERIODIC_ASSIGNMENT.md),
[Ewald reference](docs/PERIODIC_EWALD_REFERENCE.md),
[engine verification](docs/PERIODIC_LAMMPS.md) and
[bounded protocol](docs/PERIODIC_PROTOCOL.md). The default 200 integration steps
are a preparation pilot, not evidence of equilibrated polymer density.

Review a campaign and resume a locally retained native checkpoint:

```bash
polychemmd density-plan examples/density_campaign.json --output artifacts/density-plan.json
polychemmd periodic-resume artifacts/assigned-cell.json --checkpoint artifacts/preparation-result.json --checkpoint-dir artifacts/preparation --protocol examples/periodic_continuation.json --run-dir artifacts/continuation-plan --plan --output artifacts/continuation-plan.json
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_restart_v06.py --output-dir artifacts/restart-v06
```

The restart auditor compares uninterrupted and split continuation of the **same**
cell. Fragments and frames are not independent material samples. Binary restart
files stay with local run artifacts; public evidence retains their hashes and
reproduction recipes. A native restart requires the same engine build/platform
and preserves the original NPT fix state without reseeding velocities.

`density-sampling INPUT --protocol examples/density_sampling.json` accepts
explicit physical-unit thermo rows and predeclared duration/statistical gates.
Passing those diagnostics does not establish phase or full chain relaxation.
See [restart](docs/RESTART_PROTOCOL.md), [sampling](docs/DENSITY_SAMPLING.md),
[campaign planning](docs/DENSITY_CAMPAIGN.md) and [actual restart evidence](evidence/restart_v06/README.md).
The earlier example plan retains the mismatch between the legacy 512-atom limit,
10 angstrom cutoff and its declared candidate liquid-density range. It does not
silently shorten the cutoff or report the dilute pilot as a liquid prediction.

Construct and parameterize the bounded larger-cell example:

```bash
python scripts/verify_large_cell_v07.py --output-dir artifacts/large-v07
polychemmd tile-cell evidence/large_cell_v07/common/source_geometry.json --recipe examples/tiled_decane.json --output artifacts/ordered-cell.json
polychemmd parameterize-cell artifacts/ordered-cell.json --compact --output artifacts/compact-cell.json
polychemmd periodic-singlepoint artifacts/compact-cell.json --run-dir artifacts/large-plan --plan --output artifacts/large-plan.json
polychemmd density-plan examples/density_campaign_large.json --output artifacts/large-campaign-plan.json
```

The fixed motif example has 72 decane molecules, 2304 atoms and a
26 by 30.6 by 30.6 angstrom box. Its approximately 0.69876 g/cm^3 density is a
construction input/output identity, not a prediction. Repeated motifs form one
ordered starting cell; neither a liquid nor a crystal is inferred. The 12
angstrom cutoff fits within half every box axis. Linked bins verify geometric
contacts without enumerating every global pair. Compact assignment stores only
excluded/scaled pairs and explicit unit defaults; local molecular validation
still has bounded quadratic work.

`periodic-singlepoint` provides actual engine output without claiming independent
reference parity for a large cell. `periodic-lammps` retains its independently
verified small-cell role; the pure-Python Ewald reference remains capped at 512
atoms. Both representations work with bounded preparation and native restart.
The larger-cell auditor compares fixed-geometry cutoff and reciprocal-bound
sensitivity before short preparation/continuation. It does not establish bulk
convergence, phase identity or equilibrium density. See
[compact assignment](docs/COMPACT_PERIODIC.md),
[molecular tiling](docs/MOLECULAR_TILING.md),
[engine boundary](docs/COMPACT_ENGINE.md) and
[actual large-cell evidence](evidence/large_cell_v07/README.md).

For the separately installed Linux/x86-64 engine smoke environment:

```bash
python -m pip install -e '.[md-smoke]'
polychemmd lammps-run examples/lammps_smoke.in --run-dir artifacts/new-smoke --expect smoke.data --output artifacts/engine-result.json
```

The smoke script uses **LJ reduced units**. Its density is not g/cm^3 and is not
a polymer-density measurement. The current hosted-container run uses single-rank
`UCX_TLS=self OMP_NUM_THREADS=1`; MPI/GPU execution on another host needs separate
verification. A summary output cannot overwrite files inside its engine run
directory. Scripts are trusted recipes, not sandboxed input programs.

In a prepared offline environment:

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

Alternatively, `PYTHONPATH=src python -m unittest discover -s tests -v` runs the
source tests on POSIX shells without installation. The `python -m polychemmd`
entry point is equivalent to `polychemmd` after installation.

The `validate` command remains metadata-only. `build-linear` invokes separate
molecular-graph validation and assembly. The torsion example is an analytic numerical example,
not a force field fitted to PE or any other polymer. Probabilities refer to an
isolated one-dimensional angular grid; discrete Shannon entropy is grid-dependent
and is not thermodynamic entropy. All angles are radians and molecular energies
are kJ/mol; imported force-field conventions require explicit conversion.

## Development sequence

1. Freeze the attachment/atom-lineage contract (TASK-001).
2. Implement and verify a mapped-repeat-unit linear-chain builder (TASK-002).
3. Define parameter functional forms and validate coverage/convention conversion.
4. Build a real external-engine adapter and test failures as well as success.
5. Validate one end-to-end physical observable before expanding its family.
6. Benchmark chemically conditioned torsion/sequence representations separately
   from ordinary implementation features and from post-simulation analysis.

Each step needs executed evidence. Graph verification and a real process smoke
test do not establish a polymer production MD protocol or experimental
predictivity. TASK-006A now supplies bounded coordinate/cell preparation.
TASK-006B adds attributed parameter assignment and engine energy/force parity for
the restricted n-alkane graph. TASK-006C1 now adds periodic model assembly,
independent Ewald/engine verification and a bounded preparation protocol.
TASK-006 remains open until an appropriate state/preparation, production
sampling, independent cells and external comparison are established. TASK-006C2A supplies restart and sampling/planning foundations. TASK-006C2B
adds a bounded compact representation and larger geometric cell verification.
The next step is TASK-006C2C: qualify preparation history, numerical settings,
independent cells and physical sampling for a declared state.

## Reproducibility and dependencies

The lightweight core has no runtime package dependencies. Separate CI jobs target
the core, optional RDKit chemistry, and a bounded LAMMPS integration smoke. See
VALIDATION.md for executed results and the GitHub checks for hosted status;
creating a workflow file alone does not mean it passed. Optional environment
versions are pinned to those exercised for this milestone. Never commit large
trajectories, access tokens or private laboratory data.

The public-release license is not selected in this development starter; decide
and add an appropriate license with the project owner before a public release.
Existing algorithm, parameter and software sources must always be attributed.
