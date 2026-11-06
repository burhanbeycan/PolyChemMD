# Executed restricted alkane energy/force parity

These are actual outputs of:

```bash
OMP_NUM_THREADS=1 UCX_TLS=self python scripts/verify_alkane_v04.py --output-dir artifacts/alkane-v04-final
```

The three cases use the checked-in recipe and its **predeclared** absolute
tolerances: 1e-6 kJ/mol for every energy component and 1e-5 kJ/(mol angstrom)
for every Cartesian force component. Each invokes actual LAMMPS 22 Jul 2025
Update 4 with `run 0`, not an MD trajectory. Sources, versions, numerical
conventions, inputs, hashes and actual error values are in `report.json`.

| Case | Atoms | Proper dihedrals | Maximum force-component difference, kJ/(mol angstrom) |
| --- | --- | --- | --- |
| Ethane | 8 | 9 | 1.4921397450962104e-13 |
| Butane | 14 | 27 | 1.6342482922482304e-13 |
| Decane | 32 | 81 | 5.104751160445176e-6 |

Maximum energy-component difference across these cases was
2.842170943040401e-14 kJ/mol. **Not every force agrees to roundoff.** The decane
geometry has one nearly trans C-C-C-C torsion. The pinned LAMMPS OPLS backend
floors the sine used in a derivative ratio near planar torsions; this explains
the localized force difference. Its diagnostic, affected atom IDs and source
are retained in `decane/parity.json`. The geometry, reference equations, seed
and acceptance tolerances were not changed to suppress the discrepancy.

Each case contains:

- `geometry.json` and `assignment.json`: actual coordinates, atom lineage,
  frozen sourced parameters, full topology and fixed charges.
- `input.lammps` and `system.data`: the actual generated single-point inputs.
- `log.lammps` and `forces.dump`: byte-preserved engine scientific artifacts.
- `parity.json`: actual component/force comparisons, engine identity/hash and
  numerical diagnostics; only host-specific paths and incidental help
  diagnostics are replaced, as described in its public-record note.

Private process directories are excluded. Public artifact hashes were checked
after copying. Report source-module hashes identify the code actually used;
timestamps and run timing mean byte-identical envelopes/logs are not expected
on a fresh execution. Use a new output directory when reproducing the recipe.

This verifies assignment identity, conventions, derivatives and engine export
for neutral explicit **n-alkanes in isolation**, using unshifted finite-cutoff
LJ/Coulomb interactions. It does not validate long-range periodic bulk
treatment, a relaxed geometry, equilibrium, long-chain polymer density,
transport, a QM-fitted charge model or experimental accuracy.
