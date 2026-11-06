# Density protocol: implemented analysis, pending polymer simulation

TASK-006 is **not complete**. This revision supplies original mass/volume
conversion and density-series diagnostics. It does not prepare or run a polymer
cell, assign a production force field, or establish an externally validated
polymer density. An engine smoke test on Lennard-Jones particles is engine
integration evidence, not polymer-density validation.

## Implemented API

```python
from polychemmd.analysis import analyze_density, density_from_mass_volume

rho = density_from_mass_volume(masses_g_mol, volume_angstrom3)
diagnostics = analyze_density(
    densities_g_cm3,
    block_size=100,
    discard=1000,
    min_blocks=8,
    drift_relative_limit=0.02,
)
```

The geometric conversion is

`rho [g/cm^3] = sum(atomic molar masses [g/mol]) / (N_A * volume [angstrom^3] * 1e-24)`.

The Avogadro constant is exactly `6.02214076e23 mol^-1` in SI. Molecular masses
and cell volume must be supplied for the same finite cell. This identity does
not test a simulation protocol.

The series API accepts positive finite numeric densities at a uniform sampling
interval. It discards the explicitly requested prefix and reports an omitted
incomplete trailing block. The mean and variability use only complete retained
blocks. At least eight complete blocks are required by default; a user can
choose a minimum of five or more.

## What the statistical output means

- Report the sample mean, observed block-mean standard deviation and conditional
  standard error `s(block means)/sqrt(number of blocks)`.
- The optional nominal 95% interval uses a normal approximation **assuming
  independent block means**. The result explicitly says `calibrated: false`.
  Neither independence nor coverage follows from choosing a block size.
- Report split-half relative drift, adjacent-block Pearson correlation and the
  change in standard error on doubling the block size. They are diagnostics,
  not automatic convergence proofs or estimates of effective sample size.
- The default 2% drift limit is a visible, configurable review heuristic, not a
  validated physical convergence threshold. Fix suitable acceptance limits
  before a final reference comparison.
- Constant density, or identical block means despite within-block variation,
  produces an inconclusive result with no uncertainty interval. Fixed-volume
  density may be exactly constant without establishing material equilibration.
- Every output retains `equilibrium_status: not_assessed` and
  `polymer_density_validation: not_performed`. A nonflagged series is labeled
  `diagnostics_only`, never a validated material-property prediction.

## Required evidence before TASK-006 can close

1. Supported chemical connectivity, chain/cell coordinates and preparation
   lineage, including independent construction seeds.
2. A fully assigned, compatible and sourced production parameter set. The
   topology-coverage contract alone is insufficient.
3. Actual engine, version, inputs, outputs and current-run artifacts under a
   documented ensemble, temperature and pressure.
4. An NPT production window justified independently of the data selected for
   its final comparison. Check block length, autocorrelation, preparation,
   chain length and cell size; assess chain relaxation separately.
5. Independently prepared replicate cells with their separate results and a
   suitable experimental/reference comparison matched to chemistry, phase,
   conditions and preparation history.

Longer polymer campaigns must have an explicit resource budget. Short smoke
tests cannot meet the production sampling requirements by being repeated many
times, and correlated trajectory frames are not independent cells.

## Method references

- SI defining constants and Avogadro constant: https://www.bipm.org/en/measurement-units/si-defining-constants
- Molecular-simulation error estimates and correlated sampling:
  https://doi.org/10.1063/1.457480 (Flyvbjerg and Petersen, 1989).

These references support the conversion/statistical context. The implementation
here is independent and the simple diagnostics do not implement an automatic
Flyvbjerg-Petersen blocking plateau selection algorithm.
