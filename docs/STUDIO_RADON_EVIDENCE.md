# External oligomer calculations and experimental records

## Optional RadonPy process

The Studio launches one registered external calculation in its own Python process. The core package and its numerical kernels have no RadonPy dependency. The reference workflow calls RadonPy's linear graph builder for eight- or user-selected two- to sixteen-repeat neutral chains, with one methyl substituent at each model-chain attachment. The PEO model uses these caps rather than the PEG-diol termination. PVA uses fully hydrolyzed repeat units in this workflow. Stereochemistry at initially undefined stereocenters is one embedded realization. The separate surface builder controls partially hydrolyzed PVA.

Use a separate Python 3.11 or 3.12 virtual environment. With a POSIX shell, from the project root:

```bash
python -m venv .radonenv
.radonenv/bin/python -m pip install "radonpy-pypi==0.2.11" "rdkit==2026.9.1"
export POLYCHEMMD_RADON_PYTHON="$PWD/.radonenv/bin/python"
```

PowerShell equivalent after creating the virtual environment:

```powershell
.radonenv\Scripts\python.exe -m pip install "radonpy-pypi==0.2.11" "rdkit==2026.9.1"
$env:POLYCHEMMD_RADON_PYTHON = (Resolve-Path .radonenv\Scripts\python.exe).Path
```

Start the server and worker with that variable in their operating environments. The interpreter path is operator configuration, not a web request parameter. Without it, the application checks the current interpreter. The diagnostic verifies imports and reports the detected package versions. The optional Psi4 flag indicates discoverability only. Neither Psi4 nor LAMMPS is needed for this registered reference calculation.

A run uses RadonPy graph assembly, a seeded RDKit ETKDGv3 conformer, MMFF94 minimization and a RadonPy GAFF2 atom-type/parameter-assignment audit. It records an SDF, the bond graph and coordinates in Å, individual GAFF2 atom types, the explicit calculation model, software versions, UTC time and checksums. Force-field typing uses `charge=None`. The resulting SDF is a geometry/identity model, and has no validated polymer–surface cross interactions, compatible solvent/ion model or assigned numerical partial charges. The same job never silently launches a bulk property or hydrated interface simulation. Reported experimental outcomes and molecular structure calculations occupy separate views.

A live isolated RadonPy 0.2.11 / RDKit 2026.09.1 smoke test (seed 17, DP8) constructed and typed 64 atoms for each of PEO and PVA, 144 for PVP, and 136 for PEtOx. The PEO and PVA graphs have the same chemical formula under these caps, with different connectivity. A reproducibility regression tests repeated DP3 PEO SDF bytes with fixed software, verifies atom maps and excludes implicit partial-charge claims. Molecular size checks count heavy atoms, validate formal neutrality and check convergence and finite coordinates. Different numerical library builds can affect conformers. These calculations do not determine physical properties of a polymer film.

RadonPy method: Hayashi et al., *npj Computational Materials* (2022), [10.1038/s41524-022-00906-4](https://doi.org/10.1038/s41524-022-00906-4). [Upstream code](https://github.com/RadonPy/RadonPy), [published distribution](https://pypi.org/project/radonpy-pypi/).

## Literature observations

The three compact CSV files inside `polychemmd/studio/data/` are numeric secondary transcriptions and source locators from the prior MXene–Polymer Atlas. The source figures and publisher artwork are not bundled. The JSON API `/api/atlas` exposes source DOI, panel, specimen/condition, method, unit and numerical observations:

- 30 A750 observations span five redispersion media for five polymers and an untreated control from Park et al., [10.1126/sciadv.abl5299](https://doi.org/10.1126/sciadv.abl5299). The two reported digitization bounds describe reading error, not biological or batch replication. SW is the source's specific alkaline mixture. A750 is absolute absorbance. The unspecified panel time prevents direct conversion to 24 h optical retention.
- Five source-reported pairs of bare and coated dry/as-fabricated basal repeats and one wet pair appear separately. Each pair retains its original chemistry, composition and source. Xing et al., [10.1002/eom2.12300](https://doi.org/10.1002/eom2.12300), Table 2, is a dry/water-immersed source-matched quartet.
- Three washed-composite retained polymer mass fractions from Di et al., [10.1039/D4SC04935G](https://doi.org/10.1039/D4SC04935G), are TGA observations. PDDA and PAH share a mass-feed protocol. PIL has a different feed.

## Measured endpoints

For XRD, enter four source-matched basal repeats in nm per independent material batch. The program gives the coated-minus-bare difference separately for the dry and wet states and `(coated_wet − bare_wet) − (coated_dry − bare_dry)`. The included Xing quartet gives +0.060, −0.112 and −0.172 nm, respectively. It is labeled **published secondary analysis**; it supplies one study-level paired contrast. For a fresh measurement select the user-rows mode and provide raw d values. Aggregate outputs include the number of batch pairs, mean and sample SD when at least two numeric pairs are available. An unresolved diffraction peak must remain unquantified; enter only resolved, complete matched quartets.

For optics, enter time-zero and 24 h absorbance A and matched specimen-free background B per batch, medium and polymer. Inputs must follow the same measurement route, wavelength, vessel, sampling height and calibrated linear intensity range; metadata should be saved with the source experiment. `R24=(A24−B24)/(A0−B0)` requires a strictly positive background-corrected baseline. Values outside [0,1] are retained for scientific review. Contrasts pair PCE minus PEO within independent batch and medium. Conditions with an incomplete pair remain visible. Batch-level variation is distinct from digitization and technical replication. The program cannot distinguish oxidation from sedimentation from an optical ratio alone.
