# Actual periodic Ewald and bounded preparation evidence

Executed with the checked-in recipe:

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_periodic_v05.py --output-dir artifacts/periodic-v05-audit
```

`report.json` identifies the actual source/parameter/input hashes, software,
random seed, resource limits and copied artifact hashes. The source is the
four-C10H22, 128-atom starting cell in a 24 angstrom cube. Its prescribed-volume
construction density is not a bulk density prediction.

## Numerical checks that passed

The original cell and a translated boundary-crossing representation both passed
actual LAMMPS periodic energy/force comparison. Absolute tolerances were declared
before execution: **1e-4 kJ/mol** for each energy component and
**1e-4 kJ/(mol angstrom)** for each force component.

- Maximum component-energy difference: 6.96818874530436e-5 kJ/mol.
- Maximum force-component difference: 2.5901448821308293e-5 kJ/(mol angstrom).
- Translating the same state changed energies by at most about 2.16e-12 kJ/mol
  and forces by at most about 1.67e-11 kJ/(mol angstrom), below the separately
  declared 1e-7 invariance tolerances.
- The translated representation includes 42 atoms with nonzero image flags and
  ten bonds crossing periodic faces. It is **not an independent replicate**.

The model is explicit 3D periodic, tin-foil Ewald: real cutoff 10 angstrom,
alpha 0.30/angstrom, reciprocal bound 16, geometric LJ mixing, no LJ energy shift
or tail correction, and the frozen OPLS pair exclusions/1–4 rules. The reference
uses a finite integer cube; the cited engine screens a sphere. Finite-list
bounds are recorded separately from any infinite-sum convergence claim.
Pinned backend erfc and near-planar OPLS numerical approximations are recorded;
neither the reference equations nor tolerances were changed to imitate them.

## Pilot execution and unresolved sampling

The actual preparation executed CG minimization followed by **100 NVT + 100 NPT
steps**, with a 0.25 fs timestep: 50 fs (0.05 ps) of MD, plus separate minimizer
iterations. The process and final identity/box/bond/unit audits passed.

| Quantity | Requested control or limit | Actual recorded outcome |
| --- | --- | --- |
| Minimization | At most 200 iterations / 1000 evaluations | 200 iterations / 400 evaluations; `budget_exhausted` |
| Temperature control | 300 K | Final instantaneous value 190.8469707 K |
| Pressure control | 1 bar | Final instantaneous value -397.0057981 bar |
| Canonical density | No validated density target | Final mass/volume value 0.0684353798999353 g/cm^3 |
| Native engine density | Pinned real-unit conversion | 0.06843537387703591 g/cm^3 |
| Equilibrium | Requires separate sampling evidence | `insufficient_sampling`; `production_ready: false` |

Instantaneous temperature/pressure fluctuate; these final values are not an
equilibrium estimator or a general convergence test. The short duration and
declared preparation scope do not establish the target thermodynamic state,
phase, chain relaxation or density accuracy. Minimization also did not satisfy
its force tolerance within the given iteration budget. The observed final
density is **not an equilibrated polymer prediction**.

Canonical density uses the exact SI Avogadro constant. The pinned LAMMPS real
units use a slightly different native conversion factor; both values and tight
mass/volume consistency checks are retained, rather than conflating conventions.

## Files and reproduction

- `original/` and `boundary_translated/`: source cell, compact assignment summary,
  actual input/data, scientific engine log, force dump and public parity report.
- `preparation_pilot/`: actual protocol inputs/log, prepared data, final snapshot
  and execution/stage report.
- `report.json`: aggregate checks, translation invariance, actual step/status
  evidence, source and artifact hashes.

Full all-pair assignments are regenerated from the source cell and frozen model;
their canonical hashes are recorded in the compact summaries. Private process
and redundant validation-input directories are excluded. Only host-specific
paths/incidental help diagnostics are replaced in public JSON; scientific engine
artifacts are copied byte-for-byte. Use a new output directory for reproduction.
