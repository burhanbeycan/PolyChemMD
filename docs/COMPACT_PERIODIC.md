# Compact periodic n-alkane assignment

The schema-2 compact assignment removes global all-pair storage for the same
frozen OPLS-AA n-alkane model. It extends the supported starting-cell capacity
without changing charges, potential coefficients, graph rules or end groups.
It is an assignment and representation implementation, not a different force
field or a certificate that a cell is a liquid or equilibrated polymer.

```python
import json
from polychemmd.forcefields.periodic import (
    CompactAssignedPeriodicAlkaneCell, assign_compact_periodic_alkanes,
)
assigned = assign_compact_periodic_alkanes(
    packed_cell, max_atoms=8192, max_atoms_per_molecule=512,
    max_intramolecular_pairs=2_000_000,
)
record = assigned.to_dict()
round_trip = CompactAssignedPeriodicAlkaneCell(json.dumps(record))
```

`packed_cell` is an actual supported rigid starting-cell record, not a trajectory
frame. The inherited molecular-source checks verify geometry, global/source
atom identity, maps, masses, lineage, proper rotations, rigid distances,
wrapped/unwrapped positions, image flags, graph bonds and parameter applicability.
The numerical periodic cutoff/long-range model is a separate decision made by
an engine setup or protocol.

## Two explicitly different representations

| Contract | Schema 1 | Schema 2 |
| --- | --- | --- |
| Class | `AssignedPeriodicAlkaneCell` | `CompactAssignedPeriodicAlkaneCell` |
| Factory | `assign_periodic_alkanes` | `assign_compact_periodic_alkanes` |
| Default global atom budget | 512 | 8192 |
| Global pair representation | Every unordered pair in `nonbonded_pairs` | Unit-scale default plus `special_pairs` |
| Short topology | Complete rows mixed into all pairs | Only graph distances 1, 2 and 3 retained globally |
| Molecular assignment record | Full local assignment | Recomputed local assignment SHA256 and atom count |
| Parameters | Frozen sourced OPLS-AA subset | Exactly the same subset and source snapshot |

Schema 1 remains unchanged. In particular, the canonical SHA256 of the original
128-atom four-decane assignment remains
`07cf4aff9366ad00b4da1972ded417c7cd283a5a8037e9cfd129218a0bfaa307`.
Existing serialized objects and native checkpoint identities are not migrated or
relabeled in place. Choosing schema 2 creates a separately identified assignment.

Schema 2 uses kind `compact_assigned_periodic_opls_aa_n_alkane_cell` and does not
contain a global `nonbonded_pairs` array. Its `nonbonded_rules` states:

- Every omitted unordered pair has LJ scale 1 and Coulomb scale 1.
- `special_pairs` contains each same-molecule pair at graph distance 1 or 2 with
  both scales zero, and distance 3 with both scales 0.5.
- Every intermolecular pair is governed by the unit-scale default.
- Same-molecule graph distances 4 and above are also governed by that default.

The exception rows retain source-global atom endpoints, `same_molecule=true`,
actual graph distance and both scales. Endpoints follow the original global
atom order, rows are deterministically ordered, and each exception appears once.
This preserves the scaled-pair-once convention; exclusions and 1–4 scaling must
still be applied by the actual periodic numerical model, including the proper
long-range electrostatic correction.

## Capacity and complexity

For `N` total atoms and molecular sizes `n_m`, the exact number of unordered
global pairs is `N*(N-1)/2`. The exact intermolecular pair count is

\[
\frac{N(N-1)}{2} - \sum_m\frac{n_m(n_m-1)}{2}.
\]

These counts are calculated algebraically rather than enumerating global
intermolecular atom pairs. Counts of excluded, scaled and unscaled local pairs
come from the independently recomputed local molecular assignments. They sum
with the intermolecular count to the global pair count. Source geometry and
topology are retained, so a compact record does not silently drop molecules or
replace them by a united-atom model.

Short graph-distance neighborhoods have bounded size in the supported saturated
C/H graph. The number of stored special pairs therefore scales with atom count.
Molecular typing and actual rigid-distance verification still use local quadratic
work; this release does not claim linear complexity in arbitrarily long chain
length. Limits are checked before copying the source or assigning molecular
graphs:

| Budget | Hard supported ceiling |
| --- | ---: |
| Global atoms | 8192 |
| Atoms in one explicit molecule | 512 |
| Sum of local unordered pair counts | 2,000,000 |

Callers may reduce any budget and cannot increase it beyond this release's
supported ceiling. Sixteen 512-atom molecules would require 2,093,056 local pairs
and fail the aggregate gate despite meeting the global atom limit. These are
resource/applicability bounds, not statements about polymer molecular weight or
the minimum size needed for a bulk observable.

Full local assignments are transiently recomputed from each retained source
geometry. Their canonical SHA256 values are stored in the compact molecule
records; duplicated local parameter sets and local all-pair arrays are not
persisted again. The global parameter set is retained once. Source geometry
records can still repeat across molecules, so the JSON is not an interned or
binary trajectory format.

## Frozen serialized boundary and evidence

The public schema-2 constructor regenerates the entire expected compact record
from `source_cell` and compares canonical JSON exactly. Modified default rules,
exceptions, pair counts, local hashes, charges, parameters, topology or claimed
validation status are rejected. Metadata hashes bind reproducibility records;
they are not signatures authenticating an untrusted party's provenance.

Tests compare every small-cell pair against an independent breadth-first graph
oracle, including the unscaled longer-range intramolecular pairs. An instrumented
640-atom case rejects any accidental global quadratic pair enumeration; only
local eight-atom geometry pairs may be enumerated. Serialized-size and stored
exception counts are also checked across increasing molecule counts. Other
tests exercise malformed inputs, work budgets, image preservation and tampered
serialized records.

The independent analytical Ewald reference remains a separately bounded small-
cell verification kernel. Compact assignment does not make a pure-Python global
pair sum or reciprocal sum scalable. Larger actual engine calculations, numerical
sensitivity, cell-size assessment, phase preparation, chain relaxation and
independent-cell sampling require their own evidence. Construction density
remains input mass divided by chosen volume, never a predicted material result.
