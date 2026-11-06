# PolyChemMD

Independent polymer chemistry simulation and property-analysis research project.

**Release state: 0.5.0.dev0 — periodic n-alkane assignment, Ewald verification and bounded preparation pilots; equilibrium polymer-property validation remains open.**

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
| Periodic parameter assembly | Global/source atom and molecule identity, images, all intermolecular pairs and frozen per-molecule assignments are revalidated |
| Periodic Ewald reference and engine parity | Explicit real/reciprocal/self/exclusion terms, analytic forces, minimum-image gates and declared finite-sum budgets; no automatic convergence claim |
| Bounded periodic preparation | Reviewable minimization, NVT and NPT pilot with box/resource guards, native/canonical units and final-state checks; short execution never certifies equilibrium |
| LAMMPS subprocess adapter | Bounded execution, exclusive run directory, output hashes and thermo checks; real reduced-unit LJ smoke executed |
| Density analysis | Physical-unit conversion and conditional block diagnostics; no polymer equilibrium or external validation claim |
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
sampling, independent cells and external comparison are established. The next
step is TASK-006C2 sampling/protocol qualification under a declared compute budget.

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
