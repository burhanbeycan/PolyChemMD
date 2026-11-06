# PolyChemMD

Independent polymer chemistry simulation and property-analysis research project.

**Release state: 0.1.0.dev0 — development starter, not a validated polymer MD/QM platform.**

The long-term scope is a general-purpose connection from polymer chemistry and
microstructure to quantum-informed parameterization, molecular dynamics,
descriptors and material properties. Application examples do not narrow the
platform to a particular filler, solvent, polymer family or device.

Start with [the Turkish Codex guide](CODEX_START_HERE_TR.md), then read
[AGENTS.md](AGENTS.md), [architecture](docs/ARCHITECTURE.md),
[scientific contracts](docs/SCIENTIFIC_SPEC.md) and [bounded tasks](docs/TASKS.md).

## What actually works in this starter

| Capability | Status and boundary |
| --- | --- |
| Material metadata records | Implemented: schema, fraction, microstructure label and state validation; no molecular graph/valence validation |
| Provenance envelopes | Implemented: software/Python versions, UTC time and exact input-file SHA-256 |
| Complete Fourier torsion kernels | Implemented and analytically tested: all coefficients, phases and multiplicities retained |
| Torsion descriptors | Grid extrema/barrier/minima and isolated-model Boltzmann bin probabilities; illustrative coefficients only |
| Dependency inspection and CLI | Implemented; presence inspection does not execute or validate engines |
| Property catalogue | 19 planned material-property families; the catalogue is not an executable property engine |
| Polymer graph/coordinate construction | Planned |
| Force-field assignment and QM/MD execution | Planned |
| Experimental validation and publication benchmarks | Planned; no performance or material-property result is claimed |

RadonPy and other polymer workflow packages are neither dependencies nor vendored
source. The current runtime core uses only the Python standard library. Later
general libraries and external solvers must remain clearly separated from this
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

In a prepared offline environment:

```bash
python -m pip install --no-deps --no-build-isolation -e .
```

Alternatively, `PYTHONPATH=src python -m unittest discover -s tests -v` runs the
source tests on POSIX shells without installation. The `python -m polychemmd`
entry point is equivalent to `polychemmd` after installation.

The example PE record does not build polyethylene. Its mapped SMILES is stored,
not chemically validated. The torsion example is an analytic numerical example,
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

Each step needs executed evidence. The current analytical tests are not evidence
of a working polymer builder, a valid MD protocol, or experimental predictivity.

## Reproducibility and dependencies

The starter has no runtime package dependencies. The CI specification targets
Python 3.11/3.12 on Linux and does not require an MD/QM executable. Hosted CI has
not been executed merely by creating its workflow file. Later RDKit and solver
environments need their own version locks and chemistry-specific integration
tests. Never commit large trajectories, access tokens or private laboratory data.

The public-release license is not selected in this development starter; decide
and add an appropriate license with the project owner before a public release.
Existing algorithm, parameter and software sources must always be attributed.
