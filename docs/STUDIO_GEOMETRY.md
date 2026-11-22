# Surface geometry and polymer identity in the MXene studio

The geometry job builds a starting model in ångströms. Its `structure.json` contains the cell lengths, periodically tiled slab, unwrapped atom-mapped polymer template, the polymer covalent bonds, actual termination counts, and the accessible *geometric* volume. The JSON box is periodic along x and y. Each slab ID and face is retained in atom records. `structure.xyz` is a viewable preview; `polymer.sdf` contains the corresponding chain with bond orders and atom maps. Its coordinates agree with the preview to SDF rounding. A build displays one chain template even if multiple copies are requested; packing places the requested copies.

## Ti₃C₂Tₓ starting coordinates

The basal lattice is a hexagonal Ti₃C₂ basis with a three-plane titanium and two-plane carbon stack. Terminations occupy CCP hollow sites above the middle-Ti registry according to Li et al., section 3.1, [ACS Omega 2022](https://doi.org/10.1021/acsomega.2c05913). An ASE integer cell transformation creates an orthorhombic box with lengths `nx × a` and `ny × √3a`, containing `2 × nx × ny` Ti₃C₂ units per sheet.

The imposed basal spacing and atomic layer heights are derived from the deposited 47%-OH mixed aqueous interface of Armstrong et al., [ACS Nano 2026](https://doi.org/10.1021/acsnano.5c20046). Source member `Stack3x3_H2Osys_OH_47.data` has SHA-256 `ec194edb8b6a5d1923269a1dee4bdf2ba047964ad930a5d52a887baa112fa873`. The in-plane length is `54.8927/18 Å`; the vertical separations below are the upper-sheet upper-face medians of the corresponding atom-type z coordinates minus the upper-sheet middle-Ti median. This gives one specified idealized template for sampling compositions.

| Layer / element | Normal displacement from middle Ti, Å |
| --- | ---: |
| C | 1.18934 |
| Outer Ti | 2.41927 |
| OH oxygen | 3.46238 |
| F | 3.48834 |
| OH hydrogen | 4.36813 |

The O–H bond length is 0.940 Å (source bond type 9). A seeded in-plane direction supplies the remaining projection. The bottom face is the inversion partner of the top basal stack. OH and F fractions are assigned independently on each face; discrete counts use nearest-integer rounding with half-up ties. The patch mode selects the sites nearest a seeded center under in-plane minimum-image distances. Random mode samples a fixed number of sites without replacement. These constructions impose layers, occupancy and local registry; vacancies, oxidation, edge chemistry, dynamic protonation and structural relaxation require separate models. No partial charges or force-field atom types are inferred for the generated substrate.

The single-surface water height begins above its highest termination nucleus. For the slit, `gap_A` is the distance between the inner faces' outermost termination nuclei, and `center_spacing_A` is the resulting center-to-center *initial* repeat. `accessible_volume_A3` is base area times gap (or water-column height); it has no solvent-radius or polymer displaced-volume correction. A chosen spacing is not a predicted XRD result. Actual dry/wet diffraction comparison requires an independently relaxed structure or measured trajectory and a matched control.

## Finite chain definitions

The common PolyChemMD repeat-unit contracts and neutral linear builder assemble a valence-checked graph. The studio applies the specified capping or pendant-group modifications, sanitizes valences, and checks connectedness and formal charge. Neutral chains are embedded with ETKDGv3 and relaxed in their chosen isolated organic force field (MMFF94 or UFF). The core geometry verifier checks bonded distances, finite positions and retention of specified stereochemistry. The written isomeric SMILES and atom maps identify the realized microstructure. The intrinsic isolated-template energy in kJ/mol is not an adsorption or solvation free energy.

| Identifier | Chemical construction |
| --- | --- |
| PEO | Diol `HO–(CH₂CH₂O)ₙ–H`; DP counts C₂H₄O repeats. |
| PVA | Vinyl-alcohol backbone with one H and one methyl cap. Acetates replace selected hydroxyl protons. Residue indices are 1-based in `acetate_residues_1based`. |
| PVP | N-vinyl-2-pyrrolidone backbone; intact five-membered lactams. H and methyl caps. |
| PEtOx | Poly(2-ethyl-2-oxazoline) in its ring-opened tertiary-amide form with ethyl groups on acyl carbon; methyl head/H tail caps. |
| PEI | Neutral linear polyethyleneimine, methyl head/H tail caps. The protonated polymer in a specific buffer requires a different chemical model. |
| PCE | Defined free-acid acrylic acid/methoxy-oligo(ethylene glycol) acrylate comb. Grafts occur at residues 1, 1+`graft_every`, ... . Each graft is `COO–(CH₂CH₂O)ₘ–CH₃`, where m=`side_dp`. Commercial PCE architecture, distribution, ionization and actual PEG molar mass must be characterized independently. |
| bare | No polymer atoms and zero chain copies. |

PVA hydrolysis specifies a target fraction. The actual degree for a finite chain is reported after rounding the number of acetate sites. In clustered mode these sites form a contiguous block; in random mode sites are sampled without replacement. Vinyl polymers use one seeded realization of still-unspecified stereogenic backbone centers. The reported graph provides the exact stereochemical sequence; it does not estimate batch tacticity.
