# Attributed OPLS-AA assignments for finite n-alkanes

TASK-006B first supports a deliberately restricted chemical graph: one neutral,
all-explicit, acyclic n-alkane containing at least two carbon atoms. The carbon
skeleton must be a single path. Every carbon has four single bonds and is CH3
at an endpoint or CH2 internally; every H has exactly one carbon neighbour.
Methane, branches, cycles, unsaturation, aromatic atoms/bonds, heteroatoms,
charged atoms, explicit isotope labels and nontrivial bond stereo are rejected.
Unsupported chemistry never falls back to a guessed type.

```python
from polychemmd.forcefields.alkanes import assign_alkane

assigned = assign_alkane(geometry, max_atoms=256)
record = assigned.to_dict()
```

`geometry` is an immutable `MolecularGeometry` with explicit atoms and bonds.
Assignment does not require RDKit or an MD engine. It identifies environments
from the complete graph, not from polymer labels, SMILES substrings, atom names,
coordinates, or repeat-unit provenance. Graph ordering does not change the
chemical assignment. The input geometry and full lineage remain unchanged.

## Frozen parameter and charge model

The exact model is `OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET`, transcribed as selected
numerical facts from the official GROMACS v2025.3 snapshot at commit
`820db9e499e4edc7f0dc734e6baeb0e9c3f87859`. The packaged parameter data record
source files, source locations, conventions and hashes. They are not an original
fit, a generic claim about every OPLS version, or a copy of an upstream workflow.

| Environment | Assigned atom type | Fixed partial charge / e |
| --- | --- | --- |
| Terminal saturated CH3 | `opls_135` | -0.18 |
| Internal saturated CH2 | `opls_136` | -0.12 |
| Explicit H bonded to alkane C | `opls_140` | +0.06 |

Carbon and hydrogen parameter masses are 12.011 and 1.008 g/mol. Supplied masses
must agree within the explicit absolute tolerance of 1e-6 g/mol; assignment does
not silently replace a geometry's masses. No isotope-specific mass policy is
inferred. Charge assignment is the attributed empirical OPLS atom-type model,
not RESP, ESP, Gasteiger or a QM calculation performed by PolyChemMD. Total
partial charge must sum to zero within 1e-12 e.

## Complete graph topology

The constructor retains every bond, all unordered neighbour pairs about each
central atom, and **every** proper four-atom path about each C-C bond. Hydrogen
endpoints are essential: ethane has nine HC-C-C-HC torsions, and n-butane has
27 proper paths (one C-C-C-C, ten C-C-C-H and sixteen H-C-C-H), not just its one
backbone torsion. There are no impropers for this saturated alkane subset.
Reversed paths are the same interaction and are not duplicated.

Assigned bonded IDs are `bond_CC`, `bond_CH`, `angle_CCC`, `angle_CCH`,
`angle_HCH`, `torsion_CCCC`, `torsion_CCCH` and `torsion_HCCH`. Their parameter
definitions carry source references. Atom assignment records retain their source
atom ID/map, selected type and partial charge. Per-interaction records link atom
IDs to a parameter type; the parameter set and source snapshot accompany them.

Every unordered atom pair also receives an independently determined shortest
graph distance. Distances one and two are excluded from LJ and Coulomb pair
energy; distance three uses 0.5 scaling for **both** LJ and Coulomb; larger
distances use unit scaling. Each pair appears once. The source model uses
geometric sigma and epsilon mixing. These are OPLS-specific scaled 1-4 pairs:
the generic warning about excluding 1-4 interactions for some other
Ryckaert-Bellemans models must not be applied to this OPLS model.

The existing assigned-type coverage contract must report complete coverage.
Its generic `chemical_typing: not_performed` field intentionally refers only to
that coverage checker. Actual graph typing is separately reported in the outer
`validation` object. Coverage is not parameter accuracy or sampling validation.

## Resource and evidence limits

The default 256-atom limit bounds explicit all-pair enumeration. An invalid or
exceeded budget fails before assignments are returned; callers can explicitly
choose a different positive integer limit. Each output is immutable JSON text,
and `to_dict()` returns fresh nested data. No timestamps or execution state are
invented by assignment.

The public serialized constructor `AssignedAlkane(payload_json)` independently
regenerates the frozen-model assignment from the embedded geometry before
accepting a round trip. Changed charges, coefficients, topology, pair scaling,
source identity, validation claims or unexpected fields are rejected rather
than inheriting the model's label. Invalid, ambiguous or nonfinite JSON is
rejected. The constructor uses the same 256-atom default budget; a larger
positive `max_atoms` can be requested explicitly. Coordinate-only changes can
remain valid assignment inputs, but neither their geometry nor supplied source
provenance gains an additional physical-validation claim. The normal assignment
factory avoids redundant regeneration through a private validated-data path.

This layer produces a sourced assignment for its supported finite graph. It
does not evaluate energy/forces, export an engine file, relax coordinates,
choose cutoffs or long-range electrostatics, prepare a thermodynamic state,
equilibrate chains, or establish accurate polymer density/transport. Therefore
`production_ready` remains false. Static energy/force parity is a separate
numerical check; material-property validation requires a later protocol and
appropriate independent evidence. The historical alkane model's transferability
to long-chain/bulk polyethylene is not established by successful assignment.

## Primary source locations

- [Bonded numerical data at the pinned snapshot](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/ffbonded.itp):
  bond rows 192/196; angle rows 873/881/929; torsion rows 1596/1599/1824;
  improper context 2547-2549.
- [Nonbonded numerical data](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/ffnonbonded.itp):
  atom types 135/136/140 at rows 145/146/150.
- [Defaults and 1-4 conventions](https://github.com/gromacs/gromacs/blob/820db9e499e4edc7f0dc734e6baeb0e9c3f87859/share/top/oplsaa.ff/forcefield.itp):
  rows 18-20.
- [Official functional forms and OPLS/RB conversion context](https://manual.gromacs.org/documentation/2025.3/reference-manual/functions/bonded-interactions.html).
- Jorgensen, Maxwell and Tirado-Rives, JACS 118 (1996), 11225-11236,
  https://doi.org/10.1021/ja9621760 . The pinned implementation snapshot, rather
  than this citation alone, defines the exact selected numerical variant.
