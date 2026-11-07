# Compact periodic model export and bounded engine evaluation

The compact periodic assignment changes how the same frozen OPLS-AA n-alkane
model is represented. It does not change its force-field values, molecular
identities, topology, exclusions, 1–4 scaling, mixing rule or electrostatics.
`engines/periodic_alkane.py` accepts both validated assignment classes. On the
same small source cell they generate identical LAMMPS scripts and data files.

The scientific setup validates the special short-topology image/cutoff regime
by visiting only the compact record's `special_pairs`. Default unscaled pairs
are represented implicitly and are handled by LAMMPS neighbor lists; the exporter
does not reconstruct a quadratic global pair table. Existing preparation and
native continuation protocols share this exporter and its validation boundary.
The bounded exporter accepts at most 8192 atoms. This ceiling is a software
resource limit, not a demonstrated material length or time convergence regime.

## Independent reference remains a small verification kernel

`evaluate_periodic_alkanes` accepts either representation but still supports
at most 512 atoms. Its additional `max_pair_evaluations` budget defaults to
130816, the number of unordered pairs in a 512-atom cell, and cannot exceed
that hard ceiling. Pair and reciprocal-work limits are checked before numerical
evaluation. A larger compact assignment produces an explicit unsupported-size
error; the kernel does not silently attempt a large Python Ewald calculation.

For a compact input within those limits, the kernel streams unordered numeric
pairs and obtains their scales from the sparse exceptions or the explicit
default. It never materializes a new global pair array. Bonded, LJ, real,
reciprocal, self and special-pair correction equations are unchanged. Tests on
the same small source compare complete energy components, forces and pair
counts with the original exhaustive representation. Those tests verify
representation equivalence, not bulk force-field accuracy.

## Actual single-point execution without a reference claim

```python
from polychemmd.engines.periodic_alkane import run_periodic_alkane_singlepoint

result = run_periodic_alkane_singlepoint(
    assigned, run_directory="runs/large-singlepoint",
    timeout_seconds=60, cutoff_angstrom=10,
    ewald_alpha_per_angstrom=0.30, reciprocal_nmax=16,
)
```

The run directory must be new. This command executes one actual run-zero
evaluation, bounded to at most 120 seconds. It requires fresh engine artifacts,
exactly the prescribed finite thermo row, atom/molecule/image/coordinate
identity in `forces.dump`, a component-energy sum consistent with the reported
potential energy, and an actual `singlepoint.data` header confirming real units.
Canonical energies and forces are explicitly converted from LAMMPS kcal/mol
and kcal/(mol angstrom) to kJ/mol and kJ/(mol angstrom). Nonfinite conversions
fail and retain a failed output-audit record.

Its report contains energy components, ordered forces, source atom IDs,
assignment and generated-input hashes, actual engine metadata and unit/state
audit outcomes. It deliberately records
`independent_reference_status="not_performed"`,
`native_parameter_readback_status="not_performed_by_singlepoint_adapter"` and
`numerical_convergence_status="not_established"`. The single-point adapter
checks the supplied frozen assignment and generated input; its native-data
header check alone is not a full coefficient/charge/topology read-back audit.
The preparation/native-continuation layer has a separate full read-back audit,
which can also be applied explicitly by a reproducibility runner.

Success means `completed_singlepoint`, zero MD integration steps and no
material-property validation. Engine failures and output-audit failures retain
their distinct phases. An unavailable binary or stale existing directory
cannot be converted into a numeric result.

## Numerical sensitivity boundary

Larger-cell tests may evaluate the **same** coordinates with a larger cutoff or
reciprocal bound and compare component energies and atom forces. These are
paired numerical evaluations of one state; they are not independent prepared
cells, density estimates or a thermodynamic finite-size study. A comparison
must keep its original tolerances and retain failed or unresolved outcomes.
Two finite settings alone do not establish infinite-sum or bulk convergence.

The selected model remains three-dimensional periodic Ewald with explicit
alpha and reciprocal index limits, geometric LJ mixing, unshifted finite LJ
cutoff and no LJ tail correction. It does not silently switch to PPPM, another
force field or a smaller cutoff to make the cell fit. Existing
[Ewald reference documentation](PERIODIC_EWALD_REFERENCE.md) and
[engine verification documentation](PERIODIC_LAMMPS.md) describe the equations,
finite cube versus engine sphere difference and pinned backend approximations.
The primary engine contracts are
[kspace_style](https://docs.lammps.org/kspace_style.html),
[kspace_modify](https://docs.lammps.org/kspace_modify.html), and
[write_data](https://docs.lammps.org/write_data.html).
