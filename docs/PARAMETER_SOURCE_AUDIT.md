# Frozen alkane numerical source and conversion audit

The first assigned model is the selected CT/HC numerical subset distributed
with **GROMACS v2025.3**, commit
`820db9e499e4edc7f0dc734e6baeb0e9c3f87859`. It is not an original fit and is not
defined by an unspecified generic OPLS-AA name. Full upstream files were read
for verification; they are not vendored or imported by PolyChemMD. The packaged
JSON is an independently structured transcription of the selected numerical
facts, source locations and explicit conversions.

The fetched source text was checked against these actual upstream Git blobs:

| File in `share/top/oplsaa.ff/` | Git blob SHA-1 | Selected rows |
| --- | --- | --- |
| `ffbonded.itp` | `944e9c82913333def13ab6b17a7bb5d8248016fc` | Bonds 192/196; angles 873/881/929; torsions 1596/1599/1824 |
| `ffnonbonded.itp` | `20e5b31d9562f4c368c3b9ece318ed5892712e5b` | Atom types 145/146/150 |
| `forcefield.itp` | `e69c4728d7117662a32fa392409e751ea6caf03d` | Defaults 18–20 |

Full SHA-256 hashes and pinned URLs are in
`forcefields/alkane_parameters.json` and accompany each assignment. Per-type
parameter definitions identify the particular file/row. The selected records
use terminal CH3 `opls_135`, internal CH2 `opls_136` and H `opls_140`.

## Convention checks

| Quantity | Frozen upstream convention | Project convention |
| --- | --- | --- |
| Bond length | nm | angstrom: multiply by 10 |
| Bond prefactor | `k/2 * (r_nm-r0_nm)^2` | `K * (r_A-r0_A)^2`, with `K=k/200` |
| Angle prefactor | `k/2 * (theta-theta0)^2` | No-half `K=k/2`; radians internally |
| Torsion | RB polynomial in `cos(psi)`, `psi=phi-pi` | Complete Fourier terms with trans `phi=pi` |
| OPLS Fourier amplitude | `V_n/2`, kcal/mol | `K_n=4.184*V_n/2`, kJ/mol; odd phase 0, even phase pi |
| Nonbonded mixing | Geometric sigma and epsilon | Preserved |
| Bonded pair distances | 1/2 excluded, 3 scaled | Both 1–4 LJ and Coulomb scaled by 0.5, once |

The exact CCC angular coefficient is 488.273 kJ/mol/rad² in the upstream
half-prefactor convention, hence 244.1365 in the project convention. Replacing
it with a value recomputed from a rounded kcal/mol table would silently change
the selected snapshot. Hydrogen-ended proper torsions must all be retained:
ethane has nine, not zero.

The chosen HC-CT-CT-HC RB row explicitly records its November 1999 update.
Consequently the selected numerical subset must not be presented as a verbatim
transcription of the 1996 paper's initial tables. Later OPLS variants, including
OPLS/2020, also change numerical values; they are not silently mixed into this
model. In particular, agreement with an engine cannot establish that historical
OPLS-AA accurately predicts long-chain/bulk polyethylene properties.

`tests/test_alkane_data.py` checks upstream-form versus canonical bond/angle
energies, full RB/Fourier equivalence and the explicit charge/LJ/1–4 rules.
The independent evaluator and LAMMPS exporter then compare actual energies and
forces for the same declared isolated finite-cutoff verification model.

## Primary sources

- [Pinned bonded data](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/ffbonded.itp)
- [Pinned nonbonded data](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/ffnonbonded.itp)
- [Pinned defaults](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/forcefield.itp)
- [GROMACS functional forms and OPLS/RB conversion](https://manual.gromacs.org/documentation/2025.3/reference-manual/functions/bonded-interactions.html)
- Jorgensen, Maxwell and Tirado-Rives (1996), [original OPLS-AA paper](https://doi.org/10.1021/ja9621760).
- [Source-author alkane refinement and OPLS/2020 context](https://pmc.ncbi.nlm.nih.gov/articles/PMC9939004/).

These sources establish numerical identity and method context. No new QM charge
calculation, parameter fit, production trajectory or experimental validation was
performed by transcribing the subset.
