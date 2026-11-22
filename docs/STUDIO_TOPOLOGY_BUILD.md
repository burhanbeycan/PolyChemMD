# Molecular dynamics from a newly packed interface

The topology job compiles the atoms, coordinates and mobile bond graph of a completed Studio packing. Its output is registered against the SHA-256 of that exact packing. A second job can run LAMMPS on it. Changing a sequence, surface pattern, salt content or molecular packing creates a new system.

## Running the workflow

Use the full Studio installation with the Armstrong IFF source prepared. Moltemplate 2.23.0 and its OPLS/2020 table are required for polymer-containing cells. Select a zero-salt PVA/VAc or PEO/PEG-diol system or a bare hydrated control. Prepare both periodic lateral edges longer than twice your chosen cutoff. At the default 8 Å cutoff the 6 × 4 Studio cell satisfies this condition. The 12 Å setting requires correspondingly larger lateral dimensions. Add vacuum headspace above an open solvent surface using `vacuum_padding_A`; this preserves assigned water and nominal solvent volume. The pilot applies an open z boundary.

In the GUI, pack the system and select the completed packing in the topology builder. Use the resulting topology ID with the same packing in the registered MD panel. Select zero MD steps for static forces, or a positive number (up to 200) for an NVE integration. Before integration, the optional fixed-substrate conjugate-gradient minimization accepts 0–100 iterations. The chosen initial velocities apply after minimization.

To replay the packaged interface recipe from an activated full environment:

```bash
polychemmd studio --data-dir ./artifacts/interface interface-run examples/studio/designed_interface.json
```

The working artifacts include `typed.data`, `profile.json`, `atom-types.csv`, `topology-build.json` and a generated `oligomer/polymer.lt` when a chain is present. The CSV maps LAMMPS atom IDs to source atom types, charges, molecule identifiers and project atom keys. The JSON records the packed and geometric input fingerprints, the native and OPLS sources, selected site occupancies and polymer microstructure. Geometry previously assembled by the operator has its own [topology registration workflow](STUDIO_TOPOLOGY_REGISTRY.md).

## Surface, water and cross interactions

The native reference file is `Stack3x3_H2Osys_OH_47.data`. The compiler reads and verifies its mass, charge and nonbonded coefficients. It uses the two independent Ti types for inner and outer titanium. In the immobilized model the four carbon environments can share one nonbonded type because their Lennard–Jones parameters and charges in this source are equal. Each surface site must have one F or one O–H group. In the source each of those terminations has total charge −0.500 e. A Ti₃C₂T₂ unit is neutral at every permitted OH/F occupancy. Partial charges, neutral molecule totals, explicit atom identities, composition and the original bonding are audited before registration. [Native model and source installation](STUDIO_IFF.md) · [Armstrong et al.](https://doi.org/10.1021/acsnano.5c20046)

This is a frozen surface field, with two-dimensional periodicity and long-range electrostatics. The surface exerts its Lennard–Jones and Coulomb forces on mobile atoms. An explicit pair matrix uses geometric mixing of epsilon and sigma for mobile–mobile and mobile–surface interactions. Substrate–substrate LJ terms are set to zero and the rigid substrate has no bonded terms. The electrostatic self contribution of the fixed charged stack remains in the total energy and is constant for a given frozen configuration. Compare forces on mobile atoms, within-system energy changes or quantities with an explicitly matched reference. An absolute potential-energy subtraction across different surfaces includes different substrate offsets.

Deposited flexible three-site water uses the native O–H bond and H–O–H angle. Optional NaCl, oxide terminations, vacant or reactive sites, edges and mobile crystal mechanics require additional model scope. The current compiler requires zero requested salt, zero realized salt and zero ions, including boxes where a small requested concentration rounded to zero pairs. A bulk pH or salt condition cannot be inferred from these neutral, fully occupied nonreactive systems.

## Polymer chemistry

The compiler classifies atoms on the explicit-H molecular graph. A local Moltemplate run with strict force-field checks assigns the valence coefficients. Each chemically identical copy has its own coordinates and molecule ID. The compiler checks the source template against every copy before transferring parameters, including distinct OAc positions, termini, atom map, bond list and bond order. OPLS Fourier dihedrals, CVFF impropers and 0.5 scaling on 1–4 electrostatic and LJ interactions are retained. [OPLS/2020 source parameters](https://doi.org/10.1021/acs.jpcb.3c06602) · [Moltemplate](https://doi.org/10.1016/j.jmb.2021.166841)

- PVA/VAc includes alcohols, acetate esters and their connecting carbon environments. The two established DP12 sequence fixtures reproduce their archived atom-type assignments. New sequences are typed from their actual bond graphs.
- PEO uses the Studio’s HO–(CH₂CH₂O)ₙ–H geometry, with diol end groups and separate dialkyl-ether environments. Source transfer assumes neutral solution chemistry.
- A single neutral capped PEtOx repeat can be compiled. A longer PEtOx chain requires a sourced inter-repeat NM–CT–CT–NM torsion. The strict assignment stops on this missing term. Other families can use an independently prepared, complete operator topology.
- Multiple copies in one packing currently share a chemical template. Parameter and site completeness are checked per input. Equilibrium adsorption, bulk properties and the predictive accuracy of these mixing choices need independent evidence.

## Force and coordinate acceptance

LAMMPS reports the largest Euclidean **mobile-atom** force, in kcal mol⁻¹ Å⁻¹, before and after optional minimization and at the final step. Thermodynamic totals use kcal/mol for the simulation cell. Minimization has a finite step budget, and the achieved force is reported even if further relaxation would be needed. The initial temperature is specified for the free NVE segment. A geometry with a smaller nearest-atom distance can have a different force distribution depending on chemical identity and orientation.

A dedicated end-step snapshot gives the final coordinates joined by stable atom ID even when the visualization stride misses the last MD step. Original packed and final snapshot hashes are recorded independently. LAMMPS loss of atoms or a final atom outside the nonperiodic boundary fails the pilot. The measured plane separation remains imposed by the substrate constraint.

With the source cache populated and full dependencies active, run:

```bash
POLYCHEMMD_TEST_REGISTERED_LAMMPS=1 python -m unittest discover -s tests -p 'test_studio_topology_build.py' -v
POLYCHEMMD_TEST_REGISTERED_LAMMPS=1 python -m unittest discover -s tests -p 'test_studio_registration.py' -v
python -m unittest discover -s tests -p 'test_studio_polymer_opls.py' -v
```
