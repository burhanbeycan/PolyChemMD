# Periodic assembly of sourced n-alkane assignments

TASK-006C1 joins the existing per-molecule sourced OPLS-AA assignments into an
orthorhombic periodic starting cell. This independently written layer uses no
cheminformatics backend or simulation engine. It selects neither an electrostatic
method nor a relaxation/ensemble protocol; those are separate explicit layers.

```python
from polychemmd.forcefields.periodic import assign_periodic_alkanes

assigned = assign_periodic_alkanes(cell_record, max_atoms=512)
record = assigned.to_dict()
```

`cell_record` is the payload returned by `pack_cell(...).to_dict()` or the
`pack-cell` CLI result's `payload`. A result envelope is not silently treated as
a packed cell. Known additive CLI preparation metadata is preserved, together
with other JSON source metadata, but does not become verified evidence.

## Input identity and geometry checks

The source must explicitly identify its schema, angstrom coordinates, zero-origin
orthorhombic box, positive finite lengths/volume and three periodic axes. Atom,
molecule and bond counts must match their arrays. Global atom IDs/maps are unique;
the separate source IDs/maps identify the original finite molecule. Molecule
indices are explicit and source/global atom correspondence is checked exactly.

Each molecule's `MolecularGeometry` is reconstructed and passed through the
existing restricted `assign_alkane` graph/type/charge contract. Consequently all
the earlier restrictions remain: neutral explicit n-alkanes with at least two
carbons, CH3/CH2/HC environments, source-compatible masses, no branches, rings,
unsaturation, ions, heteroatoms or isotope labels. The source graph digest is
checked for agreement between its records, not treated as independent proof of
an unseen parent structure.

Every original element, isotope, formal charge, mass and complete lineage record
must agree with its packed atom. Every source bond must be present once with its
correct global atom IDs/maps and chemical attributes. Missing, extra or
intermolecular covalent bonds are rejected. The source molecule's recorded rigid
placement must contain an orthonormal rotation with determinant +1, the actual
source geometric centroid and a finite translation. Its stored transform must
reproduce the serialized unwrapped coordinates. All intramolecular pair distances
are also checked directly, so very large finite translations cannot silently
collapse a molecule through floating-point precision loss.

Each atom must carry a wrapped coordinate inside `[0,L)` and three integer image
flags, satisfying `unwrapped = wrapped + image * L`. Whole-molecule translation,
including crossing a periodic face or changing its lattice image, is accepted
when placement and image records remain coherent. Wrapped coordinates must not
be used to infer unwrapped bond vectors. Coordinate/rigidity comparisons use an
explicit absolute tolerance of `1e-8` angstrom; rotation orthogonality and
determinant comparisons use `1e-10`. Unrepresentable data fail explicitly.

The source packing's arbitrary validation flags and material/phase claims are
not trusted. Assignment does not repeat the random packing procedure, certify
its historical distance-budget claims or re-establish nonbonded contact quality.
Its own output states that source contact validity remains unassessed here.
The downstream potential evaluator must still reject singular/coincident or
otherwise unsupported configurations. This starting-cell rigidity contract does
not assert that a later thermal trajectory remains rigid.

## Assembled topology and pair conventions

The output contains one frozen parameter set and source snapshot. Each cell atom
retains the complete original packed record, with its sourced `type_id` and
`charge_e` added. Per-molecule complete sourced assignments accompany the result.
Bond, angle and proper-dihedral endpoints refer to **global cell atom IDs**;
interaction IDs are prefixed by molecule index. No bonded term crosses molecules.
All hydrogen-ended proper torsions remain present. No impropers are introduced
for the supported saturated alkane environments.

Every unordered pair of cell atoms appears exactly once. Within one molecule,
shortest graph distances one/two have zero LJ/Coulomb scaling, distance three
has one-half scaling for both terms, and longer distances have unit scaling.
Pairs belonging to distinct molecules have `same_molecule: false`,
`topological_distance: null` and unit LJ/Coulomb scaling. Topological exclusions
therefore never remove interactions merely because source atom names repeat in
different molecules. Geometric sigma/epsilon mixing remains in the frozen
parameter set; it is not changed during assembly.

The per-molecule neutral charges must yield a neutral cell. Total supplied mass,
formal charge and construction density are recomputed and checked against the
source cell's values. Construction density means mass divided by the prescribed
box volume, not an equilibrium density prediction.

## Immutable serialized boundary and budgets

The default 512-atom limit bounds explicit all-pair enumeration and direct
geometry checks. A caller may explicitly choose another positive integer limit.
No partially assigned cell is returned on an error. The output has immutable
serialized storage and `to_dict()` returns a fresh copy.

`AssignedPeriodicAlkaneCell(payload_json)` is a strict round-trip boundary: it
reconstructs and validates the embedded source cell, regenerates the complete
frozen assignment, and rejects changed charges, parameters, pair scales,
interaction topology, source snapshot or validation claims. Ambiguous/nonfinite
JSON and extra assigned-record fields are rejected. The normal factory uses a
private already-validated path to avoid repeating the same assembly work.

The output model ID is the same exact attributed numerical subset used by
TASK-006B. No new parameter fit or new physical force field is claimed. Successful
assembly establishes identities and topology; `production_ready` remains false,
with no simulation, equilibrium or material-property accuracy claim. Periodic
cutoff/long-range electrostatics, bonded images and the eventual MD protocol must
be defined and verified separately.

Source parameter identities, functional forms and numerical conversions are
documented in [the alkane assignment contract](ALKANE_ASSIGNMENT.md) and
[the source audit](PARAMETER_SOURCE_AUDIT.md).
