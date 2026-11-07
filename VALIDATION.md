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

## Coordinates and cell preparation v0.3 — executed milestone

Executed 7 November 2026 in Türkiye (6 November UTC), Linux, Python 3.12.14.

- Offline editable installation of `0.3.0.dev0` succeeded.
- **188 unit/contract tests passed, zero skipped**, including actual RDKit
  geometry generation, isotope/charge/lineage preservation, reconstruction of
  specified stereochemistry from coordinates, mirrored/inverted negative cases,
  periodic contacts, self-image bounds, image flags and XYZ views.
  [evidence/tests_v03.json](evidence/tests_v03.json) records the actual run.
- Full dependency-absent discovery under `python -S`: 188 discovered,
  **142 passed and 46 intentionally skipped**; neither RDKit nor NumPy was
  imported. These are the same tests in a different dependency environment,
  not additional independent validation cases. See
  [evidence/core_only_v03.json](evidence/core_only_v03.json).
- A new wheel was built offline, then imported from a separate temporary
  directory under `python -S`. New builder/geometry/IO modules and all 19
  catalogue families loaded from the wheel itself; no editable/source path or
  optional scientific package was used. See
  [evidence/wheel_v03.json](evidence/wheel_v03.json).
- CLI conformer generation, four-molecule periodic preparation and XYZ export
  executed successfully with the checked-in recipes. The demonstration contains
  128 explicit atoms in a user-selected 24 angstrom cube, with separate
  conformer seeds and a 1.5 angstrom geometric intermolecular cutoff.
- Independent review of an actual two-molecule PE DP3 preparation checked
  explicit neighbouring lattice images, wrapped/unwrapped identity, every
  intramolecular distance and signed-volume preservation. This was an algorithm
  check on 40 atoms, not a simulation or an additional material benchmark.
- The reusable `scripts/verify_geometry_v03.py` executed successfully. Its
  [actual report](evidence/geometry_v03/report.json) retains source/recipe/payload
  hashes and declared tolerances. The four-molecule cell passed 165888 explicit
  intermolecular image-distance evaluations, 106496 self-image distances and
  1984 intramolecular pair-distance checks. These are geometric checks, not
  independent statistical samples.
- Additional real finite-graph conformers succeeded: PE DP20 (C40H82, 122 atoms),
  PVDF DP8 (C16H18F16, 50 atoms) and PS DP4 (C32H34, 66 atoms). Their source
  recipes, seeds, coordinates and checks are in
  [evidence/geometry_v03](evidence/geometry_v03/README.md). Tacticity was not
  specified or inferred; no material-property comparison was performed.

The initial structures are finite capped oligomer graphs. Random rigid packing
does not establish an amorphous polymer phase. Construction density is
mass/selected volume arithmetic; no equilibrium density, parameter accuracy or
production MD/QM property was validated. The geometry tests use optional RDKit
ETKDGv3 as an established backend, not a newly invented conformer method.

TASK-006A is complete within its documented bounded geometry scope. TASK-006B
(production parameterization) and TASK-006C (actual polymer density validation)
remain open. Hosted Python 3.11/3.12 and real-engine regression status is reported
by GitHub checks separately from these executed local records.

## Sourced n-alkane assignment and parity v0.4 — executed milestone

Executed 7 November 2026 in Türkiye (6 November UTC), Linux, Python 3.12.14.

- Offline editable installation of `0.4.0.dev0` succeeded.
- **243 unit/contract tests passed, zero skipped**. New checks cover source-unit
  and RB/Fourier conversions, complete n-alkane graph typing, charges and pair
  scaling, analytical derivatives against independent formulas/finite differences,
  strict engine data/force parsing, malformed/tampered assignments, failed
  comparisons and the pinned backend's near-planar numerical screen. See
  [evidence/tests_v04.json](evidence/tests_v04.json).
- Dependency-absent `python -S` discovery: **197 passed, 46 intended skips**;
  RDKit and NumPy were not imported. The same suite is recorded in
  [evidence/core_only_v04.json](evidence/core_only_v04.json), not counted as a
  second set of independent physical-validation cases.
- A regular wheel was built offline and imported from a separate temporary
  directory under `python -S`. The wheel's packaged alkane JSON loaded and its
  own code assigned/evaluated 32 explicit atoms without source-module paths,
  RDKit or NumPy. See [evidence/wheel_v04.json](evidence/wheel_v04.json).
- The exact GROMACS v2025.3 CT/HC source snapshot was verified against actual
  upstream Git blobs. Source SHA-256 hashes, pinned rows, units and conversion
  rules accompany the parameter definitions. See
  [the numerical source audit](docs/PARAMETER_SOURCE_AUDIT.md).
- `scripts/verify_alkane_v04.py` executed actual LAMMPS run-zero single points
  for ethane, butane and decane using the checked-in recipe. All separate
  energy components and every force component passed the **predeclared**
  absolute tolerances of 1e-6 kJ/mol and 1e-5 kJ/(mol angstrom).
  [Actual inputs, raw scientific outputs and reports](evidence/alkane_v04/README.md)
  retain sources, identities, coordinates and artifact/code hashes.
- Maximum energy-component difference was 2.842170943040401e-14 kJ/mol.
  Force differences were about 1.49e-13/1.63e-13 for ethane/butane and
  **5.104751160445176e-6 kJ/(mol angstrom) for decane**. The latter is not
  roundoff agreement: its nearly trans CCCC torsion activates the pinned
  LAMMPS OPLS sine derivative safeguard. Independent selected-coordinate energy
  finite differences supported the analytical reference; source inspection and
  interaction localization identified the engine treatment. The fixture,
  reference equations and tolerances were retained. Each parity record now
  includes an informational, source/version-scoped geometric diagnostic.
- Independent review found and fixed a frozen-model integrity issue: modified
  serialized charges/parameters/topology/pair scales are rejected. Large finite
  engine-force mismatches also persist a finite failed record rather than
  overflowing the RMS diagnostic.

These results establish parameter identity, conversion/derivative behavior and
actual **isolated finite-cutoff** engine parity for supported neutral explicit
n-alkanes. They do not establish a periodic production polymer model, relaxed
geometry, equilibrium, bulk density, transport, QM charge fitting or experimental
predictivity. Zero MD integration steps were performed in the parity recipe.
TASK-006B is complete in its documented restricted scope; TASK-006C1 periodic
assembly/protocol preparation and the complete density milestone remain open.

Hosted Python 3.11/3.12 and actual-engine CI checks are separate from this local
execution record; workflow configuration alone is not a successful check.

## Periodic Ewald and preparation v0.5 — executed milestone

Executed 7 November 2026 in Türkiye (6 November UTC), Linux, Python 3.12.14.

- Offline editable install succeeded; **295 tests passed, zero skipped**.
  New checks cover global/source identity, image reconstruction, full periodic
  topology and pair rules, analytical Ewald derivatives, finite reciprocal
  enumeration/bounds, special-pair corrections, translated representations,
  numeric/resource gates, periodic data/force parsing and protocol contracts.
  [evidence/tests_v05.json](evidence/tests_v05.json) records the run.
- Under `python -S`, **249 passed and 46 intentionally skipped**, with neither
  RDKit nor NumPy imported. These are the same tests in a dependency-absent
  environment, not additional material-validation cases. See
  [evidence/core_only_v05.json](evidence/core_only_v05.json).
- A regular wheel was built offline. In a separate temporary directory under
  `python -S`, its own packaged code/data assigned a 128-atom cell, evaluated
  an explicitly coarse finite reciprocal sum and produced a no-execution
  protocol plan. The run directory remained absent. This verifies packaging,
  not Ewald convergence; see [evidence/wheel_v05.json](evidence/wheel_v05.json).
- The actual periodic auditor passed for the original 128-atom cell and a
  boundary-crossing translation of that same state. Maximum component-energy
  error was **6.96818874530436e-5 kJ/mol** and maximum force-component error
  **2.5901448821308293e-5 kJ/(mol angstrom)**, within predeclared 1e-4 absolute
  tolerances. Translation invariance passed the tighter 1e-7 criteria; translated
  cells are not independent preparation samples.
- Source-level review quantitatively explained residuals using pinned LAMMPS
  erfc approximations and its near-planar OPLS sine safeguard. Reciprocal engine
  energy agreed to about 4e-13 kJ/mol in the inspected case. Reference equations,
  seed, Ewald controls and declared tolerances were retained.
- The same auditor executed actual CG minimization plus **100 NVT + 100 NPT
  steps** at 0.25 fs, with time/iteration/box/bond-image guards. The process
  completed and final identity/unit/state checks passed. Minimization stopped
  at its 200-iteration budget (400 evaluations), rather than force convergence.
- Final instantaneous temperature/pressure were 190.8469707 K and
  -397.0057981 bar; these are not equilibrium averages. Canonical final density
  was 0.0684353798999353 g/cm^3 and raw native engine density
  0.06843537387703591 g/cm^3. Their different constants are explicitly audited.
  The result remains **insufficient_sampling**, **production_ready: false** and
  **density_prediction_validation: not_performed**.

[Actual public inputs, outputs and reports](evidence/periodic_v05/README.md)
retain source/model/artifact hashes and numerical contexts. The 50 fs pilot is
execution/preparation evidence, not an equilibrated polymer density or external
property benchmark. TASK-006C1 is complete in its restricted periodic domain;
TASK-006C2 sampling/state qualification and the full density milestone remain
open. Hosted checks are recorded separately by GitHub Actions.


## Restart and density review foundations v0.6 — executed milestone

Executed 7 November 2026 in Türkiye (6 November UTC), Linux, Python 3.12.14.

- Offline editable installation succeeded. **370 tests passed, zero skipped**;
  see [evidence/tests_v06.json](evidence/tests_v06.json). Meaningful new checks
  cover native checkpoint integrity, binary staging, runtime identity, evolved
  state/velocity/image restoration, typed topology multisets, density sampling,
  independent-cell accounting, capacity/resource gates and CLI input protection.
- Under `python -S`, **324 passed and 46 intentionally skipped**. RDKit, NumPy
  and LAMMPS were not imported. These are dependency-absence checks on the same
  suite, not extra material simulations; see
  [evidence/core_only_v06.json](evidence/core_only_v06.json).
- An offline wheel was exercised from a new temporary directory under
  `python -S`. Its own packaged code/data assigned the actual 128-atom source,
  planned a native continuation without reading checkpoints or creating a run
  directory, produced a campaign plan and rejected a short sampling window as
  insufficient. Packaged source bytes matched the checkout. See
  [evidence/wheel_v06.json](evidence/wheel_v06.json).
- The final actual restart auditor completed four bounded engine invocations in
  **8.886093162 s** on the measured single-rank CPU host. Preparation executed
  100 NVT + 100 NPT steps; the same checkpoint then continued for 40 steps versus
  20+20. Maximum endpoint differences were **1.0658141036401503e-14 angstrom**
  for positions, **1.93421667571414e-15 angstrom/fs** for velocities,
  **2.5501378786429996e-11 kJ/(mol angstrom)** for forces, zero for box
  components and **5.684341886080801e-13 kJ/mol** for energy, within the
  original predeclared tolerances. Atom/molecule/charge/image identities matched.
- Source/parameter/engine identities, initial restoration, NPT fix state,
  canonical/native units and output hashes were checked. All nine recorded
  source hashes match this final runtime source. Public text artifacts are
  byte-preserved; native binary checkpoint hashes are retained without committing
  the binaries. [Restart evidence](evidence/restart_v06/README.md) includes two
  initial failed read-back attempts and their source-scoped fixes. No scientific
  comparison tolerance was widened.
- Branches repeat one physical path and remain **one prepared cell**. There were
  280 integration steps across all branches; the longest path was 240 steps,
  **60 fs including preparation**. The five-frame continuous NPT window spans
  **0.01 ps** and correctly reports `insufficient_sampling`. No density
  uncertainty or external material comparison is invented.
- The executed no-engine campaign plan declares three independent cell recipes,
  paired numerical sensitivities and a resource budget. Its mass/cutoff capacity
  check rejects the candidate liquid-density range for the current small cell;
  planning does not execute the 45 requested runs. See
  [the plan](evidence/density_campaign_plan_v06.json) and
  [review of the earlier actual pilot](evidence/density_sampling_v06.json).

TASK-006C2A restart/diagnostic/planning foundations are implemented. The complete
state/sampling qualification and validated density milestone remain open.
TASK-006C2B must establish larger admissible periodic cells with verified numerical
scaling before a liquid-scale sampling campaign. Hosted checks are reported by
GitHub Actions independently of these local records.


## Compact larger-cell v0.7 — executed milestone

Executed 7 November 2026, Linux, Python 3.12.14. Required pinned RDKit/LAMMPS/
MPICH dependencies were installed in the fresh execution environment.

- **434 tests passed, zero skipped**; [local report](evidence/tests_v07.json).
  Under `python -S`, **388 passed and 46 intentionally skipped**, with no RDKit,
  NumPy or LAMMPS imported; [core-only report](evidence/core_only_v07.json).
  These are code/equation/contract checks, not additional material trajectories.
- Offline wheel creation and isolated `python -S` execution succeeded. Packaged
  source bytes matched the checkout. The wheel constructed/assigned the actual
  2304-atom geometry, planned its engine/preparation and compact campaign without
  creating run directories, and explicitly rejected the large independent
  reference calculation; [wheel report](evidence/wheel_v07.json).
- The original 128-atom legacy assignment canonical SHA256 and scientific
  LAMMPS setup/data hashes remain unchanged. Every legacy pair scale matches
  compact exceptions/defaults. Small-cell full finite Ewald energy components
  and forces match exactly; this is representation equivalence, not convergence.
- The independently checked ordered construction contains **72 C10H22 molecules,
  2304 atoms**, box **26 x 30.6 x 30.6 angstrom** and mass/volume density
  **0.6987592396767455 g/cm^3**. A separate motif-image reference establishes
  nearest periodic intermolecular contact **1.609940515 angstrom**, above the
  unchanged 1.5 angstrom geometric floor. The builder's generic linked-bin
  output remains a lower-bound check, not a global-minimum estimator.
- Compact assignment stores **12,384 special pairs** instead of enumerating
  **2,653,056 global unordered pairs**. Counts are exact; local molecular
  typing/rigidity work remains explicitly bounded. Supported capacity is at most
  8192 atoms, 512 per molecule and 2,000,000 summed local pairs. Actual engine
  evidence concerns 2304 atoms, not every system at those ceilings.
- Five actual engine invocations completed in **19.433876717 s** on the recorded
  single-rank CPU host: three same-geometry single points, bounded minimization
  plus 10 NVT + 10 NPT steps, and 10 native continuation steps. Full restoration,
  topology/parameter/unit/state audits passed. Total MD duration is **7.5 fs**;
  minimizer iterations do not add MD time.
- The first actual preparation attempt with CG200 reached its 60 s timeout
  before MD. Its failed report, input and partial log are retained. A separate
  resource-calibrated pilot was declared before its execution: CG10/eval100,
  with the same physical settings, review limits and MD steps. It stopped at
  max iterations (10 iterations, 18 evaluations), **not force convergence**.
- Fixed-geometry cutoff 10 to 12 angstrom changed total energy by
  **94.6021900313 kJ/mol**, or **0.0410599783 kJ/mol per atom**, exceeding the
  predeclared development review limit 0.02. Force-component change was
  0.00730841671 kJ/(mol angstrom). The recorded outcome is **requires_review**;
  no tolerance was widened. This changes the truncated Hamiltonian and is not
  an error estimate against an exact bulk model.
- Reciprocal bound 16 to 20 gave zero reported component-energy change and
  maximum force-component change **1.4210854715e-14 kJ/(mol angstrom)** for
  this fixed geometry. Finite-setting diagnostics passed; infinite Ewald or
  bulk-property convergence is not established.
- The continuation window has two thermo frames spanning **0.0025 ps** and
  correctly reports **insufficient_sampling**. Final instantaneous values
  (127.6384 K, 7097.8161 bar, canonical density 0.6987018285 g/cm^3) are execution
  output, not equilibrium averages or validated density predictions. The
  ordered repeated motif remains one constructed cell with unassigned phase.

[Public evidence](evidence/large_cell_v07/README.md) includes fixed construction
inputs, one shared engine data file, input/logs, three static single-point force
frames and a final velocity-bearing checkpoint frame. Other large outputs,
full assignment and native binaries are hash-only, explicitly reproducible with
the auditor. All ten source hashes and referenced public artifact hashes match.

TASK-006C2B's bounded representation/cell/engine foundations are implemented.
TASK-006C2C must qualify numerical cutoff and preparation history, then independent
state/phase sampling under measured resources. The full polymer density milestone
and external material validation remain open. Hosted checks are reported separately
by GitHub Actions.
