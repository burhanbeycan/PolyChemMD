# PolyChemMD project instructions

## Scientific scope
Develop an independent, general-purpose polymer chemistry simulation platform.
Support chemistry, microstructure, quantum-informed parameterization, molecular
dynamics, descriptors, and physical-property analysis. Preserve the broad scope;
application examples must not narrow the whole platform to one material system.

## Independent implementation
- Write project code and architecture independently. Do not copy, vendor, import,
  or wrap RadonPy, SwiftPol, HTPolyNet, or another polymer workflow package.
- Established scientific equations must be attributed. General libraries and
  external solvers (for example RDKit, NumPy, SciPy, LAMMPS, Psi4, PySCF, CP2K)
  may be dependencies, with their roles and versions documented.
- Keep literature parameters separate from original fitted parameters. Record
  origin, convention, units, applicability, and license where relevant.
- This repository is an independent implementation, not a claimed clean-room
  process: the comparative design discussion examined public upstream behavior.

## Evidence and numerical integrity
- Never invent calculations, measurements, trajectories, citations, or timings.
- Explicitly distinguish planned capability, implemented calculation, analytic
  demo, unit verification, engine integration, and experimental validation.
- Missing engines or unsupported chemistry must give an explicit unavailable or
  unsupported outcome. Never substitute plausible numeric outputs.
- Store units, phase, microstructure, preparation history, method, parameters,
  random seeds, software versions, input hashes, and uncertainty with results.
- Validate atom mapping, valence, charge, linkage chemistry, and end groups when
  building molecular graphs. Record validation does not imply chemical validity.
- A stable density is not proof that chain dynamics or a transport property has
  converged. Define sampling gates per observable.
- Separate pre-simulation descriptors from post-simulation mechanism analyses;
  prevent target and source leakage in machine-learning benchmarks.
- Torsion demo coefficients are analytic test inputs, not fitted polymer data.

## Development workflow
- Read README.md, docs/ARCHITECTURE.md, docs/SCIENTIFIC_SPEC.md, and docs/TASKS.md.
- Implement one scoped task at a time with acceptance criteria and meaningful
  scientific tests. Update the capability status and task status with changes.
- Keep domain models, numerical kernels, engine adapters, and CLI separate.
- Use explicit conventions and SI-derived units: kJ/mol for molecular energies,
  radians internally, kelvin, bar for pressure input; convert at engine boundaries.
- Use deterministic independent seeds for stochastic construction. Reject NaN,
  infinity, boolean numeric values, malformed and unsupported inputs.
- Do not launch long simulation campaigns as an incidental test. Integration
  tasks must define their system, resource budget, and expected evidence first.
- Do not add runtime LLM/API dependence to the scientific core. Codex is a
  development tool; the resulting package must run independently.

## Commands
Python 3.11 or later. Install: `python -m pip install -e .`
Offline/dev without build isolation: `python -m pip install --no-deps --no-build-isolation -e .`
Tests: `python -m unittest discover -s tests -v`
Smoke checks: `polychemmd doctor`, `polychemmd validate examples/pe_record.json`,
`polychemmd torsion-describe examples/torsion_demo.json --n-grid 720 --temperature-k 300`

## Finish a task
Report changed behavior, tests actually run, scientific limits, and remaining
work. Do not describe a runnable prototype as a validated production platform
or a manuscript as submission-ready without the required evidence.
