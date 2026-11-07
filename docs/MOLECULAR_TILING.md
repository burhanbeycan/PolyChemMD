# Ordered molecular motif tiling

`builders/tiling.py` independently constructs a deterministic orthorhombic grid
of one supplied rigid molecular geometry. It complements random rejection
packing when an explicit ordered starting configuration is useful for numerical
checks. It does not identify a crystalline phase, sample a liquid or produce an
amorphous polymer ensemble. Every copy retains the same conformer and supplied
orientation, creating a preparation bias that later state qualification must
address.

```python
from polychemmd.builders.tiling import tile_molecular_geometry

cell = tile_molecular_geometry(
    molecular_geometry,
    repeats=(2, 6, 6),
    spacing_angstrom=(13.0, 5.1, 5.1),
    rotation_matrix=None,  # identity; otherwise an explicit proper 3x3 matrix
    min_distance_angstrom=1.5,
    max_atoms=8192,
    max_distance_checks=5_000_000,
)
record = cell.to_dict()
```

The geometry is structurally revalidated independently of its provenance flags.
The supplied rotation must be orthonormal with determinant +1, preserving
handedness. Each source centroid is placed at `(i+0.5)*spacing` on each axis.
The global box has lengths `repeats*spacing`, zero origin and three-axis periodic
boundaries. Wrapped coordinates and integer image flags accompany unwrapped
molecular coordinates. Global atom IDs/maps are distinct from the preserved
source IDs/maps and lineage; copied bonds retain their chemical attributes.

The resulting object uses the existing
`orthorhombic_rigid_packed_cell` coordinate schema, so a strict parameter
assignment can independently revalidate every source-to-cell transform. The
record's method explicitly distinguishes motif tiling from stochastic packing.
`seed` and `random_generator` are null; no random generator is invoked.

## Geometry and finite work

The builder rejects nonfinite inputs, boolean numerical values, malformed
repeats/rotations, nonpositive or unrepresentable derived box volumes, coincident
source positions, and atom counts beyond the requested budget. This bounded
implementation has a hard ceiling of 8192 atoms. Every actual intramolecular
pair distance is compared before and after placement with a 1e-8 angstrom
absolute tolerance, detecting floating-point resolution loss under large
translations. The declared distance budget includes those checks and all actual
intermolecular candidate-distance evaluations. A mandatory rigidity workload
already larger than the budget fails before constructing copies.

Each rotated molecular axis extent must be no larger than the corresponding
global box length minus the geometric contact floor. This conservative
sufficient bound excludes all nonzero self images. It can reject a geometry
that a more elaborate self-image algorithm could permit.

Intermolecular contacts use a sparse periodic linked-bin dictionary. Bin widths
are at least the geometric contact floor; therefore a forbidden pair must be
in the same or neighboring bins. The same bin plus its 26 neighbors are searched
using minimum-image distances. Wrapped neighbor keys are deduplicated when an
axis has one or two bins, and same-molecule pairs are handled by the independent
rigidity/self-image checks. The number of bins per axis is capped by the atom
count; the implementation avoids overflowing `box_length/contact_floor`
ratios and corrects rounded bin widths at integer boundaries.

Spatial bins are an established neighbor-search method; see the
[primary LAMMPS neighbor-list documentation](https://docs.lammps.org/Developer_par_neigh.html)
for general context. This project implements its own Python geometric contact
audit. It does not import or copy the solver's neighbor-list implementation.

The contact workload is bounded by the supplied number of actual evaluations,
and scales with the candidates encountered, rather than materializing every
global atom pair. No universal linear-time bound is claimed for pathological
point distributions. Intramolecular rigidity checks remain quadratic in the
number of atoms per molecule.

The linked-bin pass establishes **absence of contacts below the declared
floor**. It does not measure the global nearest intermolecular distance.
Accordingly:

- `contact_distance_lower_bound_angstrom` reports the tested floor.
- `minimum_checked_candidate_distance_angstrom` reports only the smallest
  distance among candidates actually examined.
- `minimum_inter_molecular_distance_angstrom` is null, and
  `global_minimum_distance_status` explicitly says it was not measured.

The geometric contact floor is not an element-specific van der Waals radius or
a force-field overlap energy. Passing it does not prove low initial energy.
Failure returns no partial cell, writes no files and launches no engine.

## The larger decane starting geometry

The first C10H22 geometry stored in
`evidence/geometry_v03/cell.json` can be tiled with the declared example above.
Its 72 molecules give **2304 atoms** in a **26 x 30.6 x 30.6 angstrom** box,
with prescribed construction density **0.6987592397 g/cm^3**. All box axes are
longer than twice a 12 angstrom numerical cutoff. This establishes a geometric
and cell-capacity condition; the density is supplied by molar mass divided by
the chosen box volume, rather than an MD prediction.

The initially considered spacing `(14, 4.9, 4.9)` fails the unchanged 1.5
angstrom geometric floor (a nearest motif-image pair is about 1.46760
angstrom). The declared `(13, 5.1, 5.1)` construction passes. These spacings
are explicit geometry inputs selected before engine execution; no force-field
or scientific comparison tolerance is adjusted.

There is **one** constructed material cell. The 72 identical motif copies,
their molecular frames, translations and any later trajectory fragments are
not 72 independent material samples. The output retains
`phase="not_assigned"`, `equilibrium_status="not_assessed"`,
`simulation_performed=false` and `material_validation="not_performed"`.
Phase assessment, preparation relaxation, independent cells and sufficient
sampling remain necessary for any bulk density benchmark.

Tests compare small cells against an independent 27-image all-pair calculation,
exercise one/two-bin periodic cases and fixed synthetic random point sets,
check chirality/rigidity/images/identity and reject malformed or exhausted
inputs. The dense decane test uses stored geometry without requiring RDKit or
an MD executable; it verifies construction, not a material property.
