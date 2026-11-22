# Running a topology matched to a packed system

A geometry contains atom identities, connectivity and coordinates. A user-supplied parameterization adds numerical atom types, charges and coefficients with documented conventions. Registering those data locally enables a bounded LAMMPS calculation on that **specific** packing job. The Studio has separate source-pinned IFF demonstration systems. Their topology files should only be used with their exact original geometry and atom ordering.

## Quick numerical check

After installing the Studio with LAMMPS, make sure its `lmp` executable is on PATH. The four-atom example is a Lennard–Jones numerical fixture, with one fixed atom and three mobile atoms. It has no MXene, polymer or water.

```bash
python -m polychemmd.studio.register_topology add \
  --packed examples/studio_registered_lj/packed.json \
  --data examples/studio_registered_lj/typed.data \
  --manifest examples/studio_registered_lj/profile.json
python -m polychemmd.studio.register_topology list
```

Use the ID printed by `add`:

```bash
python -m polychemmd.studio.register_topology pilot \
  --source-dir examples/studio_registered_lj --id YOUR_PRINTED_ID \
  --out artifacts/registered-lj-static --steps 0
```

For the square with four sides of length 4 Å and two diagonals of length √32 Å, the canonical target is `4·U(4) + 2·U(√32) = −0.9838608497228445 kcal/mol`, where `U(r)=4ε[(σ/r)^12−(σ/r)^6]`, ε = 0.238 kcal/mol, σ = 3.405 Å, unshifted, cutoff 6 Å. A second test runs an eight-step, neutral flexible triatomic plus frozen probe to exercise PPPM slab corrections, real partial charges and harmonic bonds/angles. These systems test the engine/format translation only. To run both engine regressions, set `POLYCHEMMD_TEST_REGISTERED_LAMMPS=1` and run the Studio registration unit tests.

## Register a real study

1. Finish a GUI packing job. Export its `packed.json` unchanged. Generate a fully typed **ASCII LAMMPS data file** for those exact coordinates and IDs using a legally available, scientifically reviewed parameter set. Publishable values need independent parameter, ensemble and sampling checks.
2. Prepare a profile based on `examples/studio_registered_lj/profile.json`. Include the actual reference, scope and compatible water model. `packed_sha256` hashes the file bytes, not just the coordinates. Calculate it with `sha256sum packed.json` or, cross-platform, `python -c "import hashlib; print(hashlib.sha256(open('packed.json','rb').read()).hexdigest())"`.
3. Run the same `add` command above with those three files. The operator can set `POLYCHEMMD_TOPOLOGY_REGISTRY` to an absolute private or shared local directory; its default is `~/.polychemmd/studio-topologies`. Configure the same directory for the GUI and worker. The derived topology ID includes both the parameter file and manifest content hashes. The numerical file is copied into the local registry. The GUI accepts the ID of a registered entry and a completed packing job. It accepts no shell command, executable path or LAMMPS input script from a browser.
4. Select the **registered** MD profile, source packing job and matching ID; start with zero steps. Inspect the saved profile, generated script, input hash, initial potential, numerical forces and LAMMPS record before a short trajectory. If the data or selected packing changes, re-register the new combination. The CLI `pilot` is a local alternative for the same checks.
5. For a registered native Studio MXene geometry, the downstream trajectory selection resolves the middle Ti planes and accessible facing terminations. A single sheet selects its top O/F termination atoms. A slit selects the lower sheet's top and upper sheet's bottom O/F atoms. All slab atoms are immobilized by the present protocol. A measured wet/dry XRD change requires a separate model with mobile sheets and a converged hydration state.

## Supported parameter conventions

| Field | Required meaning |
| --- | --- |
| `schema` | `studio-full-slab-v1` |
| `label`, `parameter_source`, `applicability` | Original parameter set, citation or version, chemistry and conditions for which the operator considered it appropriate |
| `type_elements`, `bond_orders` | String LAMMPS type ID mapped to element or numeric chemical bond order, including bonds to explicit H |
| `water_model`, `water_constraints` | Human-readable model, plus `flexible` for an explicit three-site water population or `no_water` |
| `mixing` | `arithmetic` (LAMMPS geometric ε/arithmetic σ), `geometric` (geometric ε/σ), or `explicit` (complete PairIJ matrix) |
| `lj14`, `coul14` | 1–4 LJ/Coulomb factors in [0,1]. 1–2 and 1–3 are both excluded |
| `cutoff_A`, `dt_fs` | Nonbonded cutoff 2–20 Å and below half both periodic edge lengths; 0.01–1 fs pilot step |
| `packed_sha256` | Exact SHA-256 of the GUI packing artifact |

The adapter restricts the executable model to LAMMPS real units and atom_style full; orthorhombic, zero-origin, periodic in x and y, nonperiodic in z. Atom IDs are one-based in the data file and follow the order in packed.json. Å coordinates must agree to 0.00201 Å, exact graph bonds/orders and mobile molecular membership must agree, and mobile molecular and global charges must balance. Masses are checked against element identity. Three-site flexible water requires explicit bonded H–O–H angles. For polymer graphs, the basic covalent-angle inventory and proper torsion coverage of internal heavy-atom bonds are checked. This logic cannot verify chemical transferability or the numerical quality of a supplied Hamiltonian.

Allowed data sections are Masses, Pair or PairIJ Coeffs, and the coefficient and topology sections for harmonic bonds, harmonic angles, OPLS Fourier dihedrals and harmonic impropers. Excluded formats include SHAKE/rigid water, virtual sites, class-II terms, CHARMM-specific nonbonded 1–4 coefficients, arbitrary `fix` scripts, label-based coefficient files, charge equilibration and reactive models. Treat these as new adapter implementations with explicit convention tests. Check which parameter family matches the allowed functional forms; a citation or a label alone never supplies missing interactions. An added inorganic bond graph that is absent in the source packing is rejected.

The static run reports initial potential energy and per-atom force components in `engine/trajectory.dump`. Nonzero runs initialize mobile velocities at the declared temperature and integrate up to 200 steps in NVE with a fixed XY cell and fixed substrate. A charged cell uses `lj/cut/coul/long`, PPPM (10^{-5}) and slab extension 3.0; a cell with every partial charge zero uses the equivalent `lj/cut` potential. No artificial confining wall or thermostat is imposed. A free-z particle leaving the original box invalidates this bounded protocol. Logs and underlying engine records are kept, and initial/dynamics energies are reported in kcal/mol. No equilibration, adsorption free energy, dispersion lifetime or hydration-dependent basal spacing is certified by this pilot.

[Data format](https://docs.lammps.org/read_data.html) · [Nonbonded mixing](https://docs.lammps.org/pair_modify.html) · [Slab long-range electrostatics](https://docs.lammps.org/kspace_modify.html)
