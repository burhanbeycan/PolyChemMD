# PolyChemMD architecture

Status: independent development core, 7 November 2026.

PolyChemMD is intended to become a general polymer chemistry and molecular-simulation library. Its domain objects, construction algorithms, orchestration, descriptors and analysis code are developed independently. RadonPy is a comparison target and literature source, not a runtime dependency. General numerical, cheminformatics and simulation software can be external dependencies with explicit versions and attribution.

## Implemented boundary

The core contains input-record validation, mapped repeat-unit contracts, supported linear molecular-graph assembly, full Fourier-torsion analysis/conversion, assigned-parameter coverage, a bounded LAMMPS adapter and density diagnostics. A real reduced-unit LJ engine smoke was executed. It does not yet build polymer coordinates/cells, assign a production polymer force field, run end-to-end polymer MD/QM, or provide experimentally validated polymer-property predictions. Graph and process verification are separate from material validation.

The names below describe architectural responsibilities; they do not assert that all listed Python classes or modules already exist.

## Layers and ownership

| Layer | Responsibility | Main contracts |
| --- | --- | --- |
| Chemistry | Repeat units, atom maps, attachment sites, termini and reaction transformations | Explicit graph, formal charge, stereochemistry and atom lineage |
| Microstructure | Sequence, tacticity, chain-length distribution, branches and networks | Target distribution, realized distribution, sampling assumptions and seed |
| Builders | Chains, packed cells, solutions, films, interfaces and crystals | Geometry, topology, periodic cell, collision and chemical checks |
| Quantum | Conformers, charge models, torsion scans, response and interaction references | Method, basis/model, charge, multiplicity, geometry and convergence |
| Force fields | Parameter lookup, assignment, fitting, transfer and coverage | Functional forms, units, exclusions, scaling and parameter provenance |
| Engines | Isolated interfaces to LAMMPS, QM engines and later other solvers | Structured inputs, executable/version discovery and fresh artifacts |
| Protocols | Property-specific preparation, equilibration, production and recovery | Explicit state transitions and finite retry budgets |
| Analysis | Observables, uncertainties, spectra, tensors and mechanism descriptors | Validated inputs, estimator definition and sampling diagnostics |
| Informatics | Feature tables, prediction, interpretation and candidate selection | Split identity, feature availability, applicability and cost |
| Validation | Analytical checks, engine parity, numerical convergence and external evidence | Test level, acceptance margin and independent reference |

## Data flow

```mermaid
flowchart TD
    A["Chemical specification"] --> B["Microstructure realization"]
    B --> C["Coordinates and topology"]
    C --> D["QM references and parameterization"]
    D --> E["Prepared simulation system"]
    E --> F["Property-specific protocol"]
    F --> G{"Quality criteria met?"}
    G -->|No| H["Bounded recovery or failure"]
    H --> F
    G -->|Yes| I["Observable and provenance"]
    D --> J["Pre-simulation descriptors"]
    I --> K["Interpretation and learning"]
    J --> K
```

Recovery must stop when its budget is exhausted. An unresolved failure is a recorded outcome, not a request to substitute a plausible value.

## Core object boundaries

1. **Chemical specification:** identity and intended connectivity. A free monomer, polymer repeat unit and capped oligomer are distinct objects.
2. **Material specification:** polymer ensemble, phase, composition and environmental state. A repeat-unit SMILES alone is not a complete material specification.
3. **Realization:** a particular finite graph and geometry drawn from that specification; its seed and realized statistics are retained.
4. **Parameter set:** a complete potential definition with explicit unit and functional-form conventions. Atom typing is not a validation result.
5. **Run specification:** immutable scientific inputs, engine requirements, protocol and resource request. Scheduler settings do not silently change scientific inputs.
6. **Run record:** actual engine identity, start/end information, file hashes, process status and convergence evidence.
7. **Observable result:** value or array, units, method, conditions, uncertainty, sampling limitations and evidence level.

## Dependency rules

- Keep domain validation independent of installed MD/QM executables.
- Provide optional engine and cheminformatics dependencies. Importing the lightweight core must not launch an engine or require a GPU.
- Keep engine adapters separate from scientific estimators. An estimator should be testable with a controlled trajectory/reference array.
- Preserve atom maps across chemical transformations, coordinate generation and engine export/import.
- Parse subprocess arguments as argument lists. Do not construct executable shell strings from chemical labels or file names.
- Never infer successful execution from a zero exit code alone: require expected, nonempty, current-run artifacts and explicit convergence checks.
- Persist large trajectories outside Git; retain small manifests, hashes, input recipes and retrieval instructions.
- Define property-specific suitability rules. Density equilibration alone cannot certify viscosity or ion transport.

## Independent development workflow

Use one branch or worktree per bounded change. Assign exclusive file ownership when tasks run in parallel; integrate through reviewed changes. Do not let multiple tasks redesign the same shared schema simultaneously. Freeze a schema version before independent builder, engine and analysis tasks consume it.

The first dependency chain is chemical contract → mapped-chain builder → parameter contract → engine adapter → one validated end-to-end observable. Broad property families remain on the roadmap until their own protocol and acceptance evidence exist.

## References

- RadonPy, comparison software: https://github.com/RadonPy/RadonPy
- Polymer parameterization precedent: https://doi.org/10.1021/acs.jcim.3c01691
- BigSMILES: https://doi.org/10.1021/acscentsci.9b00476
- LAMMPS finite-temperature elasticity: https://docs.lammps.org/stable/Howto_elastic.html
- LAMMPS viscosity methods: https://docs.lammps.org/Howto_viscosity.html

These references establish context and reusable scientific methods; they are not evidence that the planned PolyChemMD modules have been implemented or validated.
