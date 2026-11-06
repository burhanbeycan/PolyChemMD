# Executed validation history

## Starter v0.1 — historical baseline

Date: 6 November 2026. Runtime: Linux, Python 3.12.14.

## Executed successfully

- `python -m pip install --no-deps --no-build-isolation -e .`
- `python -m unittest discover -s tests -v`: **31 tests passed**.
- `polychemmd doctor`
- `polychemmd validate examples/pe_record.json`
- `polychemmd torsion-describe examples/torsion_demo.json --n-grid 720 --temperature-k 300 --output artifacts/torsion_demo.json`
- Built a normal wheel with `python -m pip wheel --no-deps --no-build-isolation --wheel-dir artifacts/wheels .`.
- Imported that wheel from a separate temporary directory and verified its
  packaged 19-family catalogue; this avoids relying only on editable/source imports.

The tests cover analytic torsion periodicity, complete Fourier terms, known
barriers/minima, cyclic plateaus, stable normalized probabilities and invalid
inputs. Record/CLI checks cover fractions, schema/state consistency, exact input
hashes, evidence labeling, duplicate JSON fields, unrepresentable numeric values
and protection against overwriting the input through path/filesystem aliases.

A separate read-only code review identified the unrepresentable-integer and
input/output-alias issues; both were fixed and regression-tested before delivery.

## Evidence boundaries

- The example Fourier coefficients are analytic test inputs, not fitted material
  parameters. No QM, MD, experimental or ML benchmark result was generated.
- At this starter-validation stage the runtime did not expose RDKit,
  Psi4/PySCF, LAMMPS, GROMACS, xTB or CP2K. Later installations are recorded below.
- Molecular graph validation, polymer/cell building, solver adapters and physical
  property pipelines remain planned.
- The GitHub Actions workflow is provided but was not run on hosted GitHub.
- Python 3.11 is a configured CI target, not an independently executed runtime in
  this validation session.

## Chemical/execution core v0.2 — executed milestone

Executed 6 November 2026 UTC (7 November in Türkiye), Linux, Python 3.12.14.
The starter limitations above are historical; the following checks supersede
them for this development milestone.

- Offline editable installation succeeded using
  `PIP_NO_INDEX=1 python -m pip install --no-deps --no-build-isolation -e .`.
- Full unittest discovery: **132 tests passed, zero skipped**. The command-equivalent
  discovery/runner invocation and exact runtime versions are in
  [evidence/tests_v02.json](evidence/tests_v02.json).
- Optional packages actually exercised: RDKit 2026.9.1, LAMMPS wheel
  2025.7.22.4.0 and MPICH 5.0.2.
- A regular wheel was built offline and imported in a separate temporary
  directory under `python -S`, without editable/source paths or site-packages.
  Core/optional-module imports and the packaged catalogue were checked;
  [evidence/wheel_v02.json](evidence/wheel_v02.json) records the artifact hash.
- CLI graph assembly, strict parameter-coverage checking and LAMMPS Fourier
  conversion executed successfully using the checked-in examples.
- Five H-capped PE repeat graphs yielded C10H22 with complete atom lineage;
  [evidence/linear_graph_v02.json](evidence/linear_graph_v02.json) is a graph
  result, not an MD trajectory or prepared cell.
- Actual LAMMPS 22 Jul 2025 Update 4 executed the checked-in smoke recipe:
  108 LJ particles, 20 integration steps, return code zero and a completed thermo
  block. Input/output hashes, runtime controls and raw thermo values are in
  [evidence/engine_smoke_v02.json](evidence/engine_smoke_v02.json).
  Inspection of the actual `write_data` header independently confirmed LJ units.
- A separate dependency-absent review used `python -S` for the chemistry,
  builder and v0.2 CLI tests: 38 discovered, 29 skipped as intended, 9 passed;
  RDKit was not imported. This subset is not an additional 38 independent tests.
- Review findings fixed and regression-tested: explicit/implicit hydrogen role
  consistency, output-summary aliasing with engine artifacts, and `clear`
  resetting LAMMPS units. Dynamic recipe units remain unresolved rather than
  trusting a manual label.

The reduced density 0.8442 in the engine smoke is **not g/cm^3**. No polymer
equilibrium, production force field, polymer density, QM calculation or
experimental prediction was validated. Density conversion and blocking kernels
were checked against analytical/numerical cases only. TASK-006 remains open.

Hosted CI is recorded by GitHub checks, separately from this local report.
Python 3.11 and platform compatibility must be judged from the actual hosted
jobs, not inferred from the Python 3.12 local run.
