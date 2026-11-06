# Executed validation — starter v0.1

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
- This runtime did not expose RDKit, Psi4/PySCF, LAMMPS, GROMACS, xTB or CP2K.
- Molecular graph validation, polymer/cell building, solver adapters and physical
  property pipelines remain planned.
- The GitHub Actions workflow is provided but was not run on hosted GitHub.
- Python 3.11 is a configured CI target, not an independently executed runtime in
  this validation session.
