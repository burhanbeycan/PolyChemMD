# Reproduce a designed PVA–MXene interface calculation

This example checks the physical transfer from a controlled PVA repeat sequence and patterned Ti₃C₂Tₓ faces to packed atoms, a generated partial-charge topology and native LAMMPS forces. The experiment uses one sparse water box and a frozen surface. Its numerical tests address geometry, units, force consistency and integration.

From a full Studio environment with the pinned IFF source installed, run:

```bash
OMP_NUM_THREADS=1 python scripts/verify_interface_design.py --output-dir artifacts/designed-interface-acceptance
```

Use a new directory. The script creates its own job store, topology registry and two displaced single-point calculations under that directory. Full inputs, atom-parameter mapping, native LAMMPS logs, trajectories and analysis files remain available there. [Recorded output](../evidence/studio/interface_design_v2.json)

## What is checked

| Stage | Independent check |
| --- | --- |
| PVA microstructure | Read all C=O bonds back from the packed graph, map their repeat IDs to positions along the chain, count neighbouring acetates and compare with the builder's reported sequence. |
| Surface microstructure | Reconstruct the six neighbours of each OH or F site from the Cartesian coordinates with periodic x/y. Compute the correlation of adjacent OH indicator pairs and compare with the builder's face statistics. |
| Topology transfer | Match the type map to every packed atom by ID, element and molecule. Check individual molecule and sheet neutrality. Verify the digests along the packing → topology → native MD route. Compare every initial engine coordinate with the packed configuration using minimum images. |
| Energy and forces | Run a static energy evaluation. Recompute the maximum mobile force from the native per-atom force dump. Choose the polymer oxygen with the largest absolute Cartesian component, perturb that coordinate by ±0.001 Å, then compare the native force with the negative central energy difference. |
| Dynamics and analysis | Propagate mobile atoms for four 0.25 fs NVE steps. Verify fixed substrate positions and mobile displacements. Integrate the calculated water oxygen number-density profile over the simulation cell and recover the number of water molecules submitted. |

Each face has 70 termination sites in this particular cell. Measured correlations are reconstructed from the 210 distinct neighbouring pairs per face. Correlation targets are continuous inputs on a finite lattice; the report contains the achieved discrete values. The C=O check applies to this specific PVA/VAc chemistry.

The four-step trajectory and force derivative probe provide numerical verification of the assembled candidate Hamiltonian. Adsorption free energies, hydration equilibria and colloidal stability require additional sampling and independent experimental benchmarks.
