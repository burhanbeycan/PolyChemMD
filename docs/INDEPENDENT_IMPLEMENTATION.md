# Independent implementation and attribution

## Project commitment

Develop the PolyChemMD architecture, domain model, algorithms, orchestration and analysis implementation as original project code. Do not import RadonPy, vendor its source, copy its functions, or rename an existing codebase as this project. Prior software and papers may be read and cited to understand scientific methods, limitations and comparison baselines.

Independent implementation is not a claim that established methods such as Fourier torsions, RESP, Green–Kubo relations or kernel embeddings were invented here. Scientific novelty must be stated at the level demonstrated by experiments and comparisons.

## External dependencies

General-purpose Python libraries and external scientific engines are allowed dependencies, with documented roles, versions and citations. Examples are NumPy/SciPy, optional RDKit, LAMMPS, Psi4/PySCF and CP2K. Their presence does not make the PolyChemMD application code a reimplementation of those engines.

Separate three kinds of contribution:

1. **Project implementation:** original code written for this repository.
2. **Scientific method or parameter:** published equation, algorithm, model or numerical value with source attribution.
3. **External software service:** an executable or library performing a computation with a recorded version.

Published force-field parameters require their own source/license handling. Original code must not imply that imported parameter sets are newly developed or universally valid.

## Provenance and review practices

- Record literature and software sources next to the method specification.
- For parameter fitting, retain original QM inputs, actual outputs, fitting settings and independent test references.
- Make copied third-party assets or permitted snippets explicit if ever introduced; do not label them as original project implementation.
- Record substantial AI-assisted implementation and the human scientific validation process according to the eventual journal's policy.
- Attribute benchmark methods and datasets even when no source code is copied.
- Keep a contribution and release history. Do not create fabricated dates, experiments, author contributions or execution records.

## Comparison boundary

RadonPy and other packages can be run independently as benchmark tools in a separately documented comparison environment. Preserve each tool's scientific conventions and report equivalent conditions where a comparison is claimed. A lower error obtained from changing the model or reference state is not automatically evidence of a better implementation.

## Claim review

Before a research release, identify which contributions concern software engineering, which concern polymer-specific representation or parameter transfer, and which concern new physical insight. Global originality statements need a focused current literature review. Broad code ownership, a larger feature list and a passing unit-test suite are not substitutes for that evidence.
