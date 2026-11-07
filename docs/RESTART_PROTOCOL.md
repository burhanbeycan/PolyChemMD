# Native NPT checkpoint and continuation

Status: implemented bounded engine integration for the restricted periodic
OPLS-AA n-alkane model. This does not qualify an equilibrium polymer density.

The starting `AssignedPeriodicAlkaneCell` remains the immutable chemistry,
parameter and rigid-packing source. An evolved MD configuration has different
coordinates and often a different box; it is a native checkpoint with its own
lineage, rather than a replacement rigid-packed cell. Restarting preserves that
distinction and never regenerates a conformer or replaces source atom identities.

## Public API

```python
from polychemmd.protocols import (
    PeriodicAlkaneContinuationSpec,
    plan_periodic_alkane_continuation,
    run_periodic_alkane_continuation,
)

# `prepared` is a successful run_periodic_alkane_protocol result, and `assigned`
# is the same validated original assignment used in that preparation.
spec = PeriodicAlkaneContinuationSpec(steps=40, timeout_seconds=60,
                                      thermo_stride_steps=10)
plan = plan_periodic_alkane_continuation(
    assigned, prepared["checkpoint"], spec,
    checkpoint_directory="runs/preparation", run_directory="runs/continuation",
)
result = run_periodic_alkane_continuation(
    assigned, prepared["checkpoint"], spec,
    checkpoint_directory="runs/preparation", run_directory="runs/continuation",
)
```

Planning validates the manifest, scientific controls and assignment identity
without reading checkpoint files, creating directories or discovering an engine.
It explicitly defers artifact checks until execution. A continuation is limited
to 1–10,000 integration steps and at most 120 seconds. It inherits the original
temperature, pressure, timestep, cutoff and Ewald controls; changing those is a
new protocol, outside this continuation API.

## Checkpoint contract and execution gates

Preparation writes `checkpoint.restart` while `polychemmd_npt` is still active,
then writes a reviewable `checkpoint.dump` containing positions, molecule/atom
IDs, charges, image flags, forces and velocities. The binary retains native
thermostat/barostat variables and atom velocities. The manifest records both
artifacts' SHA256 values and sizes, the original assignment's canonical SHA256,
model, units, step, physical time, box, controls and actual producer identity.
Each subsequent checkpoint also names its parent's canonical SHA256 and chain
index. The chain is reproducibility metadata, not a cryptographic signature or
authentication system for untrusted checkpoints.

Before any engine integration, execution checks the exact manifest schema,
source/model identity and artifact checksums. It rejects symlink artifacts and
files larger than 64 MiB before reading them into memory. The checkpoint bytes
are staged under `source.restart` into a new exclusive run directory. The solver
never reads mutable inputs directly from an older run directory.

LAMMPS restores the native state and recreates `fix polychemmd_npt all npt` with
the original ID, style, group, constant targets and 100/1000 fs damping times.
The continuation explicitly restores Ewald and neighbor settings because these
are not all retained in a native restart. It performs a run-zero restoration
audit followed by exactly the requested steps. There is no minimization,
velocity creation, timestep reset or new stochastic seed on continuation.

The output audit requires the engine's explicit confirmation that it restored
the saved global NPT fix information. It checks source versus restored atom
IDs, molecule IDs, charges, images, positions, velocities and box. Position,
velocity and box components use an absolute tolerance of 1e-11 in their stated
units. Restored forces are independently compared with the saved full-precision
force snapshot at an absolute component tolerance of 1e-7 kJ/(mol Å). A mismatch
is a failed continuation with retained process/output records.

Two native `write_data` read-backs verify real units, masses, types, charges,
potential coefficients and bonded topology. Pinned LAMMPS coefficient writers
use C `%g`, so the text coefficient comparison uses the corresponding six
significant digit serialization. The underlying dynamics still use native
binary doubles, and the unchanged binary has its own SHA256 gate. This text
audit does not claim to resolve coefficient differences below the writer's
precision. The full-precision restored-force comparison supplies an additional
behavior check. Physical parity tolerances are not relaxed to accommodate text
serialization.

LAMMPS may renumber bonded-interaction rows in local-atom traversal order when
writing data. Those row numbers are not molecular identities. Bond, angle and
proper-dihedral audits compare full multisets of type and atom endpoints,
including multiplicity, under full endpoint reversal. Full reversal is
equivalent for this model's harmonic bonds/angles and OPLS proper terms.
Molecular atom IDs, parameter-type IDs and actual connectivity remain strict.

Runtime guards retain `Lx,Ly,Lz > 2*cutoff` and `bondmax < cutoff/3`. The latter
is a sufficient bound for the short excluded/scaled topology endpoints to
remain inside the selected minimum-image cutoff regime. Thermo and final
snapshot checks repeat the box/bond checks. Canonical exact-SI density and raw
LAMMPS real-unit density remain distinct, using the documented native constant.
The returned `thermo_rows` contains only the integration block, so the run-zero
audit row is not duplicated in sampling input.

## Native compatibility and scope

Compatibility is deliberately restricted to the same engine build, platform,
configuration and direct single-process layout. The process record stores the
entrypoint hash, engine help/configuration hash, platform and selected runtime
environment. For the recognized PyPI console entrypoint it also hashes the
installed native `lmp` and `liblammps` files, since hashing the small Python
wrapper alone would miss native engine updates. A successful help response can
include incidental timestamp/PID-bearing UCX host diagnostics; raw help and its
hash remain in provenance, while the stable configuration fingerprint begins
at the LAMMPS banner and excludes explicit timestamped UCX ERROR/WARN lines.
Changing scientific configuration still changes the fingerprint and is rejected.

The fingerprint is conservative and is not an exhaustive hash of every linked
system library. External launchers are identified by their actual entrypoint
and help configuration; their transitive environment is not fully inspected.
Native restarts are not a portable archival interchange format. The text data
snapshot supports inspection, but cannot reconstruct active NPT fix state.

Successful execution returns `completed_npt_continuation`,
`equilibrium_status: insufficient_sampling`, `production_ready: false` and
`density_prediction_validation: not_performed`. A process restart is not an
independent preparation replicate. Multiple continuation segments or correlated
frames must not be counted as independent cells. Phase qualification, adequate
physical duration, observable-specific sampling gates, cutoff/long-range and
finite-size sensitivity, and state-matched external evaluation remain separate.

## Primary engine references

- LAMMPS native restart contents, compatibility and fix ID restoration:
  https://docs.lammps.org/read_restart.html
- Native checkpoint writing and active fix information:
  https://docs.lammps.org/write_restart.html
- Nose–Hoover NVT/NPT dynamics and restart behavior:
  https://docs.lammps.org/fix_nh.html
- Pinned harmonic-angle binary versus `%g` text coefficient serialization:
  https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/MOLECULE/angle_harmonic.cpp

These are engine contracts and numerical implementation sources, not physical
validation of the force field or any polymer material property.
