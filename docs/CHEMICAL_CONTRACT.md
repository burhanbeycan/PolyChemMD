# Repeat-unit and finite-chain chemical contract, version 1

TASK-001 defines an already transformed repeat-unit graph. TASK-002 assembles
identical copies into a finite head-to-tail linear graph. This is independently
written PolyChemMD code using RDKit for SMILES parsing, molecular graph objects,
sanitization and serialization. It does not import another polymer builder.

## Chemical identity and public API

```python
from polychemmd.chemistry import RepeatUnitSpec
from polychemmd.builders import build_linear

unit = RepeatUnitSpec.from_mapping({
    "label": "polyethylene repeat graph",
    "smiles": "[*:1]CC[*:2]",
    "head_map": 1,
    "tail_map": 2,
})
chain = build_linear(unit, 5, head_cap="H", tail_cap="methyl")
payload = chain.to_dict()
mol = chain.mol  # a defensive RDKit copy; all hydrogens are explicit
```

`RepeatUnitSpec` accepts exactly the four illustrated fields; the constructor and
`from_mapping` both validate them. `to_dict()` losslessly preserves the original
input, and `to_mol()` returns a freshly parsed graph with its attachment dummies.
The two attachment map numbers must be distinct positive Python integers.
Booleans and numeric strings are rejected. Labels and SMILES must be nonempty.
An invalid graph raises `ValueError`; supported syntax outside the implemented
chemistry boundary raises `UnsupportedChemistryError`, a `ValueError` subclass.
SMILES names and CXSMILES extensions are rejected rather than silently discarded;
enhanced stereochemistry groups are outside this version's preservation contract.
An unavailable optional dependency raises `ChemistryDependencyError`, a
`RuntimeError` subclass. Importing the package, chemistry contracts or builders
does not itself import RDKit or require an MD/QM executable.

A free monomer such as ethylene `C=C`, an already transformed repeat unit
`[*:1]CC[*:2]`, and the finite H-capped dimer `CCCC` are different objects. This
builder does not infer reaction chemistry, remove leaving groups, change bond
orders, predict conversion, or model a reaction trajectory. Condensation repeat
units must already encode the intended polymer connectivity. Material metadata
still has its separate state/microstructure contract in `records.py`.

## Attachment and supported-chemistry rules

- The input is one connected graph with exactly two dummy atoms, each bonded to
  one non-dummy, non-hydrogen atom by a single bond. The dummy atoms are neutral,
  nonisotopic, hydrogen-free, and carry the declared head/tail maps.
- Every supplied nonzero atom map is unique within a repeat unit, including
  optional maps on real atoms. Unmapped real atoms retain source-index lineage.
- Copies are assembled by connecting tail anchor `i` to head anchor `i+1` with a
  single bond. Every dummy is deleted. Each remaining terminus is explicitly
  capped with either `H` (one H) or `methyl` (one C and three H). Other strings,
  including the alias `CH3`, are rejected rather than interpreted.
- Sidegroups and sidegroup rings are permitted; they do not create additional
  attachment sites. Branch/network construction, head-to-head placement,
  alternating/random/block sequences, cyclic-chain closure, aromatic/double/
  triple attachment bonds, radicals and other caps are unsupported.
- Formal charges, explicit source H atoms and isotopes are retained. Neutral
  caps add no formal charge, so product formal charge is DP times repeat charge.
  This is charge bookkeeping, not a partial-charge/protonation model or proof of
  realistic ionic thermodynamics.
- RDKit valence sanitization, one-component connectivity, removal of all dummies,
  exact atom/isotope balance, formal-charge balance and unique identities are
  checked after assembly. These are graph checks, not proof of chemical
  synthesizability, steric viability, stable geometry or correct material state.

## Stereochemistry boundary

Specified tetrahedral centers away from the two attachment anchors and internal
E/Z bonds are retained. Attachment-anchor stereochemistry is rejected because
replacement of its neighbor requires an explicit transformation-aware stereo
contract. E/Z definitions whose stereo reference atoms include an attachment
dummy are also rejected, as is non-tetrahedral atom stereochemistry. The builder
rejects a specified tetrahedral center if the chosen assembly/caps make that
center undefined or alter its retained tag. This includes a short capped chain
whose formerly distinct open arms become chemically identical. The code does
not assign missing stereochemistry or realize tacticity. Preserving repeated
input stereochemistry is not a claim of isotacticity for arbitrary repeat units.

## Atom lineage and graph serialization

`chain.lineage` is an immutable tuple of `AtomLineage` records. Each atom has:

| Field | Meaning |
| --- | --- |
| `atom_id` | Stable origin identifier, such as `repeat:0:atom:1`, `cap:head:H` or a parent ID followed by `:H:0` |
| `atom_index` | Zero-based index in `chain.mol` |
| `atom_map` | Positive unique map in the assembled molecule; currently `atom_index + 1` |
| `origin` | `repeat`, `cap` or `added_hydrogen` |
| `repeat_index` | Zero-based copy index; added H inherits its parent's copy, cap atoms have null |
| `source_atom_index` | Zero-based index in the parsed input SMILES; null for newly added atoms |
| `source_atom_map` | Original positive atom map if provided; null for unmapped/new atoms |
| `parent_atom_id` | Parent of an automatically added H atom; otherwise null |
| `terminus` | `head`/`tail` for cap atoms and their added H atoms; otherwise null |
| `structural_role` | `backbone`, `sidechain`, `terminus` or `unassigned` |

Source explicit H atoms remain source atoms and preserve their indices, maps and
isotopes. Both source H atoms and automatically completed H atoms inherit their
bonded parent's structural role, independently of explicit/implicit H notation.
Deleted attachment dummy atoms do not appear in the product lineage. Output maps
are a current-graph index; they intentionally differ from the source attachment
maps. Source indices and IDs are stable for the same parsed input/copy, not under
arbitrary rewriting of equivalent input SMILES. Atom-map lineage survives JSON
round trips via `AtomLineage.to_dict()` and `AtomLineage.from_mapping()`.

Structural roles use the unique acyclic path between the two attachment anchors.
Nonpath real atoms are sidegroups. If a path edge is in a ring, multiple backbone
routes exist: the assignment method becomes `ambiguous_cyclic_path` and all real
atoms are `unassigned`. The algorithm does not guess a chemically privileged
route through aromatic rings. `backbone_assignment` records the method and source
path, enabling later descriptors to reject or handle ambiguity explicitly.

`canonical_smiles` clears atom maps and ordinary explicit H atoms to represent
molecular identity. Isotopic H remains explicit. `mapped_smiles` preserves the
complete explicit-H output graph and unique maps. Canonical SMILES atom order
need not match the `atoms` table; joins must use `atom_map`, not serialization
position. The `atoms` table contains lineage, element, isotope, formal charge and
chiral tag; bonds refer to atom maps. `.mol` returns a defensive copy so editing
it cannot silently invalidate the stored result.

The result identifies `operation=repeat_unit_assembly` and
`evidence_level=graph_validated`. It contains no coordinates, force-field
parameters, energies, density, reaction history or physical-property predictions.

## Limits, verification and compatibility

DP, `max_degree_of_polymerization` and `max_atoms` are strict positive Python
integers. Defaults limit DP to 1000 and the complete explicit-H graph to 100000
atoms. The predicted atom count is checked before copying/assembling; callers may
explicitly choose different limits. There is no stochastic step or random seed.
Determinism is within the same RDKit version; cross-version canonicalization is
not guaranteed, and the actual RDKit version is recorded in results.

The version-1 contract is additive and does not reinterpret or migrate the
existing metadata-only `records.RepeatUnit` class. The stricter graph object is
named `RepeatUnitSpec` and explicitly requires attachment-role fields. Existing
material records remain metadata and are not retroactively graph-validated.

Executed verification belongs in the task change report. The accompanying tests
cover manually derived polyethylene formulas `C_(2n)H_(4n+2)`, methyl/mixed cap
balances, exact dimer identities, polyether linkages, sidegroups, charges,
isotopes, malformed graphs/maps/valence, deterministic origin mapping, optional
dependency absence and explicit stereochemical support boundaries. No MD or QM
engine is run by this task.

RDKit primary documentation for the library operations used here:

- [Getting started in Python](https://www.rdkit.org/docs/GettingStartedInPython.html)
- [SMILES parser and serialization API](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolfiles.html)
- [Graph operations, hydrogen handling and sanitization](https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html)
- [RDKit stereochemistry discussion](https://www.rdkit.org/docs/RDKit_Book.html#stereochemistry)
