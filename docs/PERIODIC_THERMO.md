# Static periodic pressure and energy evaluation

`engines/periodic_thermo.py` adds bounded run-zero evaluation of energy, forces
and the full static virial-pressure tensor. It accepts either validated periodic
assignment representation and reuses the existing scientific setup. The older
periodic parity/single-point adapters and their inputs remain unchanged.

```python
from polychemmd.engines.periodic_thermo import (
    plan_periodic_alkane_thermo_singlepoint,
    run_periodic_alkane_thermo_singlepoint,
)

options = dict(run_directory="runs/static-pressure", timeout_seconds=60,
               cutoff_angstrom=16, ewald_alpha_per_angstrom=0.30,
               reciprocal_nmax=20)
plan = plan_periodic_alkane_thermo_singlepoint(assigned, **options)
result = run_periodic_alkane_thermo_singlepoint(assigned, **options)
```

The plan exposes the exact input script, staged data, SHA256 values and process
request. It discovers no engine and creates no directory. Execution requires a
new exclusive directory and is limited to at most 120 seconds per process.
The returned process evidence is distinct from numerical or material validation.

## Bounded neighbor-list storage

A measured 4704-atom, 17 angstrom cutoff case exceeded LAMMPS' default
per-atom neighbor-list storage. Larger supported cells therefore receive an
explicit storage request from the shared periodic setup and native restart:

```text
one = max(2000, 2*N)
page = max(100000, 10*one)
```

For `N <= 512`, no storage command is added and earlier input bytes remain
unchanged. For `513 <= N <= 8192`, `neigh_modify one ... page ...` sets at most
16384 neighbors per atom and 163840 neighbor entries per page. These are bounded
storage allocations; a page setting is not a bound on total process memory.
They do not alter the 2 angstrom skin, rebuild policy, force-field parameters,
cutoff, Ewald settings or scientific comparison tolerances. Native continuation
redeclares the same values after reading the checkpoint rather than relying on
saved neighbor configuration. The policy supports the measured finite cases;
it is not a guarantee against every pathological geometry or periodic-image
population. An engine overflow remains an explicit retained failure. See the
[primary neighbor-list storage definitions](https://docs.lammps.org/neigh_modify.html).

## Pressure convention and static state

The new script explicitly sets every atomic velocity to zero, then uses

```text
compute polychemmd_virial all pressure NULL virial
thermo_modify press polychemmd_virial norm no format float %.17g
```

LAMMPS' `virial` selection includes the supported pair, bond, angle, dihedral
and reciprocal electrostatic contributions while excluding kinetic pressure.
The pressure is therefore a static response of the supplied configuration.
It is not an isothermal pressure estimate, a time average or an assertion that
the configuration represents a relaxed zero-temperature phase.

The scalar pressure equals one third of the Cartesian tensor trace. The native
six-component order is `xx, yy, zz, xy, xz, yz`; the API returns the corresponding
symmetric **3 by 3** matrix in `pressure_tensor_bar`. Its sign convention is
positive-compression pressure. This is not the opposite-signed per-atom stress
convention. In a periodic system the complete virial must include periodic and
reciprocal contributions; a naive sum of supplied unwrapped atom positions
times net forces is not used as an independent pressure reference.

Primary engine definitions are the [pressure compute documentation](https://docs.lammps.org/compute_pressure.html),
[thermodynamic fields](https://docs.lammps.org/thermo_style.html), and
[pressure compute selection](https://docs.lammps.org/thermo_modify.html).

## Output and identity checks

One prescribed step-zero thermo row must report zero temperature and zero
kinetic energy. A full-precision `forces.dump` contains atom ID, molecule ID,
partial charge, wrapped coordinates, image flags, force and velocity. Every
source identity, fixed charge, image and coordinate is checked. Each output
velocity must be exactly zero. These checks prevent a thermal or displaced
configuration from being silently compared as the same static state.

The adapter checks finite energies/forces, the sum of all six energy components
against potential energy, and scalar pressure against tensor trace divided by
three. The energy components retain the names `bond`, `angle`, `dihedral`, `lj`,
`coulomb_realspace`, `coulomb_kspace` and `total`. It checks thermo volume against
the actual assigned box and native density against the pinned engine's
mass/volume convention. The actual `singlepoint.data` header must confirm
`units = real`. This header check is not a full native parameter/topology
read-back; that separate audit remains explicitly unperformed by this adapter.

Canonical units are kJ/mol, kJ/(mol angstrom), bar, angstrom cubed and g/cm cubed.
Native LAMMPS real-unit energies/forces are converted with 4.184 and pressure
with **1 atm = 1.01325 bar**. Native engine density is retained separately from
canonical mass/volume density with the exact SI Avogadro constant. The pinned
native density factor is `1/0.602214129`; see
[LAMMPS real units](https://docs.lammps.org/units.html) and the
[pinned constants](https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L175).

The result carries the complete assignment's canonical SHA256 and an additional
`geometry_canonical_sha256`. The latter hashes the complete box and ordered
atom fields `cell_atom_id`, `molecule_index`, `charge_e`, `coordinates_angstrom`,
`wrapped_coordinates_angstrom` and `image_flags`, using finite, sorted-key,
compact UTF-8 JSON. This lets a numerical comparison require identical source
state across settings. Equal geometry hashes alone do not authenticate an
untrusted record or establish physical suitability.

## Evidence limits and failure behavior

Successful execution reports `completed_static_thermo_singlepoint` with
`pressure_interpretation="static_virial_not_equilibrium_average"`, zero MD
integration steps, `equilibrium_status="not_assessed"`,
`production_ready=false` and no independent reference or bulk convergence
claim. Varying cutoffs can change the truncated Hamiltonian itself; differences
in these static responses are sensitivity observations, not uncertainty in an
equilibrated material property.

Engine and output-audit failures retain separate failure phases and no invented
pressure result. Unavailable engines and stale directories do not create or
overwrite a successful result. Synthetic unit fixtures verify tensor ordering,
trace/sign/unit conventions, native/canonical density separation and independent
identity/velocity/resource failure cases. Actual engine execution is recorded
separately by the bounded audit runners.
