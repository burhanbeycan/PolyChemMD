# Bounded LAMMPS adapter

The independently written adapter executes an installed external LAMMPS binary;
it does not implement or vendor LAMMPS, a polymer force field, equilibration, or
experimental property validation.

## API

```python
from pathlib import Path
from polychemmd.engines import LammpsRunSpec, plan_lammps, run_lammps

spec = LammpsRunSpec(
    input_script=Path('examples/lammps_smoke.in').read_text(),
    run_directory=Path('artifacts/unique-smoke-run'),
    executable='lmp',  # or a known executable path; None searches normal names
    timeout_seconds=30,
    expected_artifacts=('log.lammps', 'smoke.data'),
)
plan = plan_lammps(spec)  # no files, process or computed observables
result = run_lammps(spec).to_dict()
```

An existing run directory is rejected. Use a new directory for each run. Small
data files may be supplied as immutable `(relative_name, text)` pairs through
`input_files`. Expected output names may not overlap those staged inputs.

## Execution gates

1. Resolve an executable, then actually invoke `-h`. File presence alone does not
   imply a working MPI/runtime installation.
2. Hash the executable bytes, exact executed input bytes, staged input files and
   captured help output. Save the help identity, command arguments, timing and
   selected runtime variables (`OMP_NUM_THREADS`, `UCX_TLS`) without dumping the
   whole environment. `engine_version_banner` searches for an actual LAMMPS
   banner; the first help line is retained separately and may be a host diagnostic.
3. Run with an argument list, no shell, and a finite timeout (default 30 seconds,
   explicitly configurable up to 3600). POSIX timeout handling kills the process
   group; no automatic MPI launcher is constructed.
4. Require zero process status and nonempty current-run artifacts. Symlinks and
   outputs outside the new run directory are rejected. Modification times older
   than launch by more than one second fail the freshness check; this tolerance
   accommodates coarse filesystem timestamps.
5. Require the standalone `POLYCHEMMD_RUN_COMPLETE` output marker after the final
   completed thermo block. Merely echoing its `print` command is insufficient.
6. Parse numeric, finite single-line thermo blocks beginning with `Step`; require
   configured columns (default Step/Temp), nondecreasing steps and a closing
   `Loop time of` line. YAML/multi-line thermo and `run ... post no` are outside
   this initial parser contract.

The completion marker is a recipe contract, not a security or scientific proof.
Scripts are trusted LAMMPS programs and may contain file/system operations; the
adapter is not a sandbox. Ordinary host diagnostics containing an `ERROR` word do
not replace process, artifact and thermo checks. A successful record explicitly
retains `equilibration_validation=not_performed` and
`material_validation=not_performed`.

Failures after creating a run directory save `run_record.json` with the actual
failure phase/reason. Missing executables and pre-existing directories fail
before any current-run artifacts are created. `LammpsError.record` exposes a
saved failure record when available. Dry-run plans cannot contain numeric thermo
results.

`declared_units_style` retains the optional manual `units_style` declaration.
`input_units_style` is a separate conservative recipe inference; a manual value
never substitutes for an unresolved inference. Straight-line analysis tracks
literal `units` and `clear` commands: `clear` restores LAMMPS defaults including
LJ units. The label refers to units active at `run`/`minimize`, so a later reset
cannot relabel earlier thermo rows. Different unit systems across runs leave the
scalar label null and preserve the list `recipe_run_units_styles`. A contradictory
explicit specification is rejected when the straight-line recipe is known.

Any include, jump, variable, conditional, Python/shell execution, restart load,
run-every control, variable substitution or line continuation leaves units
unresolved. This is deliberately not a full LAMMPS interpreter. Even a resolved
label has confidence `recipe_literal_not_engine_verified`, never an engine unit
validation claim. A property-specific protocol must verify its applicable unit
and dimensional model before assigning physical units to observables.

Thermo values remain raw and are not converted or silently assigned physical
units. In particular the smoke test's `units lj` density is dimensionless and
must never be treated as g/cm3. Energy normalization is not inferred from column
labels. Explicit LAMMPS `ERROR:`/`ERROR on proc ...:` log lines fail the run;
generic host diagnostics containing the word ERROR do not.

## Smoke recipe and test evidence

`examples/lammps_smoke.in` constructs 108 dimensionless Lennard-Jones particles
on an FCC lattice, runs 20 NVE steps, writes `smoke.data`, and prints the explicit
completion marker. Its numerical output is an engine integration check, not a
polymer density or equilibrated material benchmark. Run zero would also be legal
LAMMPS execution, but would provide no dynamical equilibration evidence.

`tests/test_lammps_adapter.py` uses temporary Python process fixtures to test
failure contracts. Those executables are not LAMMPS and their numeric rows are
clearly labeled test fixtures. Real installed-engine runs must be reported and
retained separately with executable provenance. If an installation needs a
single-rank MPI runtime setting, configure it externally and record it; the
adapter does not silently alter the scientific/runtime environment.

## Primary documentation

- Command-line `-h`, `-in`, `-log`, `-screen`: https://docs.lammps.org/Run_options.html
- Single-line thermo and custom observables: https://docs.lammps.org/thermo_style.html
- Run completion and legal `run 0`: https://docs.lammps.org/run.html
- Print output to screen/log: https://docs.lammps.org/print.html
- Clear resets settings to their defaults: https://docs.lammps.org/clear.html
- Unit systems and reduced LJ quantities: https://docs.lammps.org/units.html

Equations and engine behavior are attributed to the LAMMPS documentation. The
process/provenance/validation adapter is original project code.
