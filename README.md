# PolyChemMD

Independent polymer chemistry simulation and property-analysis research project.

**Release state: 0.2.0.dev0 — chemical graph and execution core, not a validated polymer-property platform.**

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
| Complete Fourier torsion kernels | Implemented and analytically tested: all coefficients, phases and multiplicities retained |
| Torsion descriptors | Grid extrema/barrier/minima and isolated-model Boltzmann bin probabilities; illustrative coefficients only |
| LAMMPS torsion imports | Explicit real-unit Fourier and OPLS conversions with independent equation tests |
| Parameter contracts and coverage | Units, forms, provenance, exclusions and 1-4 scaling; assigned-type completeness, not chemical atom typing or force-field accuracy |
| LAMMPS subprocess adapter | Bounded execution, exclusive run directory, output hashes and thermo checks; real reduced-unit LJ smoke executed |
| Density analysis | Physical-unit conversion and conditional block diagnostics; no polymer equilibrium or external validation claim |
| Dependency inspection and CLI | Implemented; presence inspection does not execute or validate engines |
| Property catalogue | 19 planned material-property families; the catalogue is not an executable property engine |
| Polymer coordinates, cell packing and production atom typing | Planned |
| End-to-end polymer MD and QM protocols | Planned |
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
```

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
predictivity. TASK-006 remains open until coordinates/cell preparation,
production parameterization, sampling and external comparison are complete.

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
