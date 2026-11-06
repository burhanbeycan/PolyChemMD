# Scientific contracts and evidence levels

## Scope

Develop an independently implemented, general polymer platform linking chemical connectivity and microstructure to molecular interactions, dynamics and observables. No property is supported merely because its name appears in a registry or a roadmap.

The starter implements generic analytical and metadata functionality. Its illustrative Fourier profiles are analytical examples, not measured properties, force-field fits or simulated data for a real polymer.

## Units and chemical identity

- Use explicit units for dimensional inputs and outputs. Canonical energy is kJ/mol and canonical angular calculations use radians; engine conversions must be documented and tested.
- A material record must distinguish chemical identity, phase, thermodynamic conditions and preparation history.
- Formal molecular charge, partial-charge model and protonation/ionization assumptions are distinct fields.
- Track parent atom identities through transformations and retain the relationship between free monomers, repeat units, capped fragments and final chains.
- Mark unknown values as unknown. Do not convert an assumption into an observation.
- Store target and realized ensemble statistics separately. Mean molecular weight and dispersity do not uniquely specify the entire chain-length distribution.

## Canonical torsion convention

The analytical convention is

\[
V(\phi)=\sum_j K_j\,[1+\cos(n_j\phi-\delta_j)].
\]

Here \(K_j\) is in kJ/mol, \(n_j\) is a positive integer periodicity and \(\phi,\delta_j\) are radians. The convention does not contain an implicit factor of one half. Importers for other functional forms must explicitly convert prefactors, phases and units. A coefficient is not, in general, the full rotational barrier of a multi-term profile.

Minimum scientific checks:

1. Analytical single-term energies and derivatives at independently specified angles.
2. Periodicity and phase equivalence, including angles outside one principal interval.
3. Additivity of multiple Fourier terms and retention of all terms.
4. Angle-unit and energy-unit conversions with known reference values.
5. Convergence of grid-dependent extrema or population summaries when grid resolution changes.
6. Stable normalization of thermal weights where such weights are implemented.

For isolated torsional coordinates, Boltzmann weights describe the specified model profile and measure. They do not by themselves establish condensed-phase polymer conformer populations, coupled-torsion entropy or a macroscopic material property. Label sampled extrema as sampled; a finite grid does not prove an exact global extremum.

## Independent descriptor research

The proposed research extensions are full torsional-landscape features, chemically conditioned backbone/sidechain/junction representations, and sequence/contour correlations. Separate the following feature sets:

| Feature set | Availability | Valid use |
| --- | --- | --- |
| Structure and microstructure | Before simulation, subject to known inputs | Fast prediction and candidate screening |
| FF/QM features | After their actual parameter/reference calculations | Screening with the stated computational cost |
| Trajectory/dynamical features | After MD and sampling validation | Mechanistic analysis or explicitly costed multi-fidelity learning |

Do not include a target observable or an equivalent transformation of that target as an input to claim independent predictive accuracy. Train/test splits must separate related records, repeated simulations and literature sources as required by the scientific question.

## Property contracts

Each implemented property must define an estimator, units, phase, ensemble, preparation protocol, required time/length regime, independent replicates, uncertainty method and acceptance criteria. Numerical tolerances must be selected before the final comparison, not adjusted to produce agreement.

| Family | Required validity checks |
| --- | --- |
| Density and equilibrium thermodynamics | State/ensemble match, stationarity, autocorrelation-aware uncertainty and independent preparation |
| Tg and thermal transitions | Preparation, heating/cooling rates, fit definition and rate sensitivity; apparent transition versus equilibrium claim |
| Elasticity | Stress/strain convention, positive/negative small perturbations, finite-temperature averaging and isotropy before isotropic reductions |
| Reology and viscoelasticity | Linear-response amplitude, diffusive/terminal regime where needed, thermostat and rate/frequency sensitivity |
| Diffusion | Unwrapped trajectories, component identity, reference frame, diffusive window and finite-size assessment |
| Sorption and permeation | Chemical potential/loading, equilibrium/metastable state, solution–diffusion assumptions and pure/mixture selectivity |
| Ionic transport | Collective correlations and reference frame; label Nernst–Einstein estimates as approximations |
| Interfaces | Definition of free energy, collective variable, sampling overlap, surface state and finite-size effects |
| Electrical/optical response | Orientational versus electronic contributions, field linearity, polarization model and double-counting control |

Ordinary fixed-topology MD does not establish bond scission, degradation temperature, electronic conductivity or triboelectric electron transfer. Additional models are required for those claims.

## Evidence levels

1. **Specified:** design and equations exist.
2. **Implemented:** runnable code exists.
3. **Analytically tested:** independent limiting cases or exact examples pass.
4. **Engine tested:** actual executable and generated artifacts have been checked.
5. **Numerically validated:** sampling, finite-size/rate and protocol checks pass for a stated domain.
6. **Externally evaluated:** comparison with independent experimental or appropriate reference evidence.

Evidence is specific to a method, version, chemistry and state. Passing a unit test does not promote a feature directly to external validation. A failed or insufficiently sampled run must remain failed or inconclusive.

## Publication evidence plan

- Compare proposed descriptors with chemistry-only, existing force-field distribution and appropriate graph baselines.
- Include matched-composition systems with changed sequence, tacticity or linkage environment when data permit.
- Evaluate unseen chemical environments and architectures, and report computation cost with prediction quality.
- Separate analytical examples, generated simulation data and experimental references in every table.
- Preserve identifiers, units, source locations, protocol versions, seeds and file hashes.
- New code is an implementation contribution; global method novelty requires a specific literature comparison and demonstrated benefit.

## Primary references

- Kernel-mean force-field representations: https://doi.org/10.1103/PhysRevB.108.134107
- General organic polymer parameterization: https://doi.org/10.1021/acs.jcim.3c01691
- Psi4 SAPT interaction analysis: https://psicode.org/psi4manual/master/sapt.html
- CP2K QM/MM: https://manual.cp2k.org/trunk/CP2K_INPUT/FORCE_EVAL/QMMM.html
- Polymer persistence-length analysis: https://docs.mdanalysis.org/2.10.0/documentation_pages/analysis/polymer.html
