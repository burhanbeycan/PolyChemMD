# Bounded periodic rigid-molecule packing

This TASK-006A component independently constructs an initial orthorhombic
periodic cell from supplied `MolecularGeometry` objects. It runs no MD, assigns
no force field and does not establish an amorphous or equilibrium polymer
state. Conformer generation and parameterization are separate layers.

```python
from polychemmd.builders.packing import pack_cell

cell = pack_cell(
    geometries,  # one immutable MolecularGeometry per whole molecule
    box_lengths_angstrom=(24.0, 24.0, 24.0),
    seed=4101,
    min_distance_angstrom=1.5,
    max_attempts_per_molecule=500,
    max_atoms=10000,
    max_distance_checks=5000000,
)
record = cell.to_dict()
```

## Construction and geometric checks

The constructor centres each molecule on its geometric centroid, draws a proper
rotation, and translates its centroid uniformly within the prescribed box.
Four independent zero-mean unit-variance Gaussian draws are normalised into a
unit quaternion. The resulting rotation samples SO(3) uniformly before the
packing rejection step. Accepted placements are **not** an equilibrium or
uniform sample of the complete packed-system configuration space; sequential
rejection introduces order and acceptance biases.

Rotations preserve distances and handedness; no reflections or internal
deformations are used. Each rigid placement is checked for bond-length
preservation. Full source coordinates and all source atom/lineage records are
retained. Input internal nonbonded contacts are not newly certified by packing.

All distinct-molecule atom pairs use the orthorhombic minimum-image distance.
The cutoff applies equally to every element and includes explicit hydrogens.
It is a user-selected geometric exclusion in angstrom, **not** a van der Waals
radius rule, force-field energy threshold or experimentally derived contact
criterion. An independent full-pair audit rechecks the accepted cell and reports
the actual minimum intermolecular distance. This quantity is null for a single
molecule, for which there is no pair of distinct molecules.

To exclude contacts with a molecule's own periodic images, every rotated
molecular axis span must satisfy `span_axis <= L_axis - cutoff`. For any
nonzero lattice translation, at least one component of every inter-image atom
pair is therefore at least `L_axis - span_axis >= cutoff`. This sufficient
condition covers all nonzero images, including translations beyond neighbouring
cells. It is conservative: long chains, winding periodic molecules and some
otherwise valid orientations are unsupported. An arbitrary orientation cannot
be used to fit a chain whose diameter exceeds the available box diagonal.

The method is intended for small initial structures and dilute/low-overlap
preparation examples. Dense melt packing or long-chain construction requires a
more capable, separately validated algorithm; raising the attempt budget does
not establish physical equilibration.

## Coordinates and identity

Each atom receives a unique `cell_atom_id` prefixed by its molecule index and a
new sequential global `atom_map`. The distinct fields `source_atom_id` and
`source_atom_map` identify the original geometry. Its original lineage object
is retained unchanged and describes that source identity; it is not presented
as the new global identity. Global bonds reference cell IDs/maps, and complete
per-molecule source geometries accompany the placements.

`coordinates_angstrom` retain the whole molecule unwrapped. Each atom also has
`wrapped_coordinates_angstrom` within `[0, L_axis)` and integer `image_flags`,
with `unwrapped = wrapped + image * L` componentwise. The box origin is zero and
all three axes are periodic. Molecular topology must never be reconstructed by
assuming the wrapped coordinate separation is the unwrapped bond vector.

## Budgets, reproducibility and density

The RNG is a private `random.Random(seed)` instance; global Python random state
is not modified. The seed, Python RNG name, PolyChemMD/Python versions and
implementation, actual placements and attempt counts
are retained. Bitwise reproducibility is expected for the same input, Python
implementation/version and platform; it is not promised across library/runtime
changes. No stochastic physical ensemble is implied by multiple packing seeds.

Finite atom, per-molecule attempt and total pair-distance budgets are enforced.
The total distance budget includes the final independent audit. Exhaustion
raises `PackingError`; no partly packed result is returned. Unsupported input
or nonrepresentable numeric ranges raise explicit errors. The naive all-pair
implementation intentionally favours auditable correctness for small cells over
large-system speed; it is not a production-scale packer.

`construction_density_g_cm3` is total supplied atomic molar mass divided by the
**user-prescribed** cell volume, using the exact SI Avogadro constant. It is an
input-dependent construction quantity, not a density prediction or MD result.
The record keeps `simulation_performed: false`, `equilibrium_status:
not_assessed`, `phase: not_assigned` and `force_field_assigned: false`. Mixtures
of supplied molecular geometries can be placed, but mixing, solubility and phase
behaviour are not determined by this geometric procedure.

## Method context

- Shoemake, K., “Uniform Random Rotations”, *Graphics Gems III* (1992),
  https://www.sciencedirect.com/science/article/pii/B9780080507552500361 .
  Quaternion sampling is an established mathematical method; this implementation
  does not copy the chapter's code.
- SI defining constants, including the Avogadro constant:
  https://www.bipm.org/en/measurement-units/si-defining-constants .

The minimum-image expression and the self-image bound are documented above so
the checked geometry can be reproduced independently of the code.
