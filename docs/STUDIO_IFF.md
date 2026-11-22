# Explicit-water Ti₃C₂Tₓ MD in the studio

The MD screen runs an independently selected, fully typed canonical system. A model packed in the geometry screen contains starting coordinates and molecular identity. To propagate a newly packed chemistry, all charges, bonded and cross parameters must be established for that chemistry. Its geometry controls are never silently applied to a canonical MD run.

## Prepare the data once

In the activated project environment on Linux or WSL:

```bash
python -m pip install -e '.[studio,studio-analysis,studio-models,md-smoke]'
python -m polychemmd.studio.iff_sources --fetch --opls
python -m polychemmd.studio.iff_sources --status
```

`--fetch` downloads the original 114 MB Supporting Information for Armstrong et al., *ACS Nano* (2026), DOI [10.1021/acsnano.5c20046](https://doi.org/10.1021/acsnano.5c20046), [original archive](https://pmc-oa-opendata.s3.amazonaws.com/PMC13450432.1/nn5c20046_si_002.zip). `--si local-download.zip` serves offline installations. The downloader verifies the complete archive and one specific LAMMPS data member against pinned SHA-256 values. The model is the two-layer `Stack3x3_H2Osys_OH_47` from the liquid-interface example.

`--opls` rebuilds the two PVA/VAc sequence fixtures using Moltemplate 2.23.0's `oplsaa2024.lt` and the published OPLS/2020 parameters ([Jorgensen et al.](https://doi.org/10.1021/acs.jpcb.3c06602); [Moltemplate](https://doi.org/10.1016/j.jmb.2021.166841)). The supplied project inputs describe mapped, individually typed oligomer geometries, matched coordinates near the outer face and 70 removed whole water molecules. The assembler calls the independently installed Moltemplate compiler, merges the resulting topology with the native solid/solvent data, and verifies the final uncompressed coordinates and parameters byte-for-byte against the earlier Atlas fixtures. Both shell helpers supplied with Moltemplate must be on PATH. This preparation was exercised with Python 3.12 under Linux. Once prepared, jobs need only the local cache, a working LAMMPS executable and MDAnalysis.

For owners of the prior PVA or full MXene–Polymer Atlas archive, the alternate command below imports the two OPLS and two MATCH-CGenFF molecular assignments. No third-party MATCH or CGenFF library is redistributed. These exact source files are necessary for the MATCH profiles.

```bash
python -m polychemmd.studio.iff_sources --atlas MXene_Polimer_Atlasi_Arastirma_Paketi.zip
```

`PVA_MXene_Arastirma_Paketi.zip` is also recognized. Prepare the bare source with `--fetch` or `--si` separately. Files live at `~/.cache/polychemmd/iff`, or set `POLYCHEMMD_IFF_CACHE` for a different location. `--status` lists installed inputs. Actual MD always rechecks each digest. A missing, corrupt or unsupported input stops the run. Native Windows users can import the prebuilt Atlas in Python; the POSIX Moltemplate rebuild path and this LAMMPS configuration were tested in Linux.

## Molecular models and boundaries

| Profile | Active atoms and composition | Candidate force-field convention |
| --- | --- | --- |
| iff-water | 5,148 fixed solid atoms, 3,456 flexible water molecules; 15,516 atoms | Deposited IFF/solvent terms, geometric mixing, 1–4 Coulomb/LJ scaling 1.0 |
| iff-opls | 5,148 fixed solid atoms, 3,386 water molecules and one 104-atom, 75%-hydrolyzed DP12 poly(vinyl alcohol-co-vinyl acetate); 15,410 atoms | Same native pair terms, geometric mixing for OPLS, 1–4 Coulomb/LJ scaling 0.5 |
| iff-charmm | As above, same initial contacts as OPLS for the same sequence | MATCH 0.01 CGenFF v2a oligomer; explicit retained native unlike pairs; arithmetic σ for polymer-containing pairs, geometric ε, 1–4 scaling 1.0 |

The dispersed sequence has acetate units 3,6,9; the clustered sequence has 5,6,7. Input atom IDs, finite molecular mass and termini describe the oligomer. Salt concentration is zero in these presets. The source surface has 342 F and 306 OH sites on the two exterior faces collectively. The chain begins on an exterior basal surface. Across profiles the frozen stack has the same deposited configuration and flexible water is propagated at constant volume. In the original native water literature recipe the solid was flexible; here we freeze it to provide a common restraint convention. The fixed central-plane separation is geometric. Predicting hydrated film spacing requires a mobile stack, appropriate hydration and experimental comparators.

## Engine and analysis

LAMMPS uses real units, triclinic periodic boundaries, 12 Å cutoff, PPPM tolerance 10⁻⁵, through 25 conjugate-gradient minimization iterations, a 298.15 K Nosé–Hoover thermostat on the mobile group, 100 fs damping and 0.5 fs step. Each job requests 1–2,000 steps. A bare pilot and both polymer profiles were executed with LAMMPS 22 July 2025 Update 4. A connected worker can extend an existing MD result by specifying `restart_job` and additional `md_steps`. Source parameters, binary checkpoint SHA-256, engine identity and the parent input digest are checked. Restart restores positions, velocities, the fixed cell and integrator/thermostat state. The same force-field coefficients are reapplied explicitly after restart. The PVA profiles are model candidates without independent aqueous adsorption or free-energy validation. Total energies can differ between Hamiltonians, including fixed stack self-terms.

The output directory retains the exact engine input, full atom table, LAMMPS log and warnings, sorted trajectory, final binary restart, per-frame hydration CSV, and the last sampled structure for the GUI. `trajectory-groups.json` maps central Ti, external termination sites, water O and oligomer heavy atoms to stable IDs. The MDAnalysis 2.11.0 check uses minimum-image triclinic distances and counts water O neighbors within 3.5 Å separately around F and OH oxygen. Both exterior faces enter the statistic. This is a geometric count over output frames, including step zero. Use independently prepared systems and longer sampled trajectories with uncertainty analysis to study adsorption. A ≤1 ps engine pilot cannot establish equilibrium, a polymer ranking, or 24 h dispersion stability. Absolute thermodynamic data from LAMMPS real units are labeled kcal/mol in these engine artifacts, with area/density observables labeled separately.

## Verification

After all five inputs are available, run in a new output directory:

```bash
python scripts/verify_studio_iff.py --output-dir ./iff-verification
```

The script runs three matched 10-step + 10-step continuations and direct 20-step controls with seed 15991. It pairs all final atoms by ID under triclinic periodicity and checks thermo and coordinate discrepancies. [Recorded measurements](../evidence/studio/iff_engine_checks.json) include hashes, frame counts, per-site water contacts and observed split errors. On the recorded Linux CPU run, maximum final displacement was below 2.71 × 10⁻⁷ Å at the printed precision; simulations were only 0.010 ps long. The code accepts 2 × 10⁻⁶ Å for the maximum final coordinate discrepancy and 10⁻⁵ K and kcal/mol for final temperature and full-cell potential energy, respectively. Hardware, MPI layout and executable configuration matter. In a restricted single-process sandbox requiring UCX shared-network isolation, launch with `UCX_TLS=self OMP_NUM_THREADS=1`. Standard single-process LAMMPS installations can use local defaults.
