# Explicit torsion convention conversion

Status: independently implemented and analytically tested conversion routines.
No engine execution, force-field assignment, material simulation, or experimental
validation is implied. The checked-in example coefficients are illustrative
analytical inputs, not parameters fitted for a named polymer.

## Canonical convention

PolyChemMD uses the entire term collection

\[
V(\phi)=\sum_j K_j[1+\cos(n_j\phi-\delta_j)]
\]

with molecular energy in kJ/mol and angles in radians. There is no implicit
one-half prefactor. Each coefficient, phase, periodicity and term order is
retained. This representation is a numerical functional-form contract, not an
atom-typing rule or a claim that a parameter is physically applicable.

## LAMMPS Fourier, `units real`

The LAMMPS expression has the same cosine convention. Its phase coefficients
are degrees and its energy coefficients use kcal/mol in `real` units. Conversion
therefore uses **1 kcal/mol = 4.184 kJ/mol** and radians = degrees × pi/180.

The input record is deliberately explicit:

```json
{
  "style": "fourier",
  "units": "real",
  "m": 2,
  "terms": [
    {"k_kcal_mol": 1.0, "periodicity": 1, "phase_deg": 0.0},
    {"k_kcal_mol": -0.2, "periodicity": 3, "phase_deg": 47.0}
  ]
}
```

`m` must equal the exact array length. Signed and zero coefficients are allowed.
Unknown/missing fields, boolean numeric values, non-finite values and numeric
conversions outside the floating-point range are errors. Unwrapped phase values
are retained, rather than silently normalized. Export fails if the requested
degree value cannot be represented as a finite number.

LAMMPS permits periodicity zero, while the current PolyChemMD core requires
positive periodicities. These constant-offset terms are **explicitly rejected**,
even if their coefficient is zero. They are not dropped or converted to an
arbitrary energy baseline. Extending the core to support constant offsets is a
separate contract change.

```python
from polychemmd.descriptors.conventions import (
    import_lammps_fourier, export_lammps_fourier,
)

potential = import_lammps_fourier(record)
canonical_grid = potential.describe(n_grid=720, temperature_k=300)
roundtrip = export_lammps_fourier(potential)
```

These are pure conversion functions. The caller must retain the original record,
source citation, parameter identity and input-file hash in the run provenance.
The conversion does not perform conformational-degeneracy or central-bond
coefficient redistribution: source parameters must already correspond to the
individual torsion being represented.

## LAMMPS OPLS, `units real`

The separately named importer uses the LAMMPS equation

\[
V(\phi)=\tfrac12K_1[1+\cos\phi]+\tfrac12K_2[1-\cos2\phi]
+\tfrac12K_3[1+\cos3\phi]+\tfrac12K_4[1-\cos4\phi].
\]

Its schema is `{"style":"opls", "units":"real", "k_kcal_mol":[K1,K2,K3,K4]}`.
The canonical amplitudes are 4.184 × K/2; phases are 0, pi, 0, pi. All four terms,
including zero coefficients, are retained. An OPLS exporter is not provided;
the generic canonical form may contain terms/phases that OPLS cannot represent.

```python
from polychemmd.descriptors.conventions import import_lammps_opls

potential = import_lammps_opls({
    "style": "opls", "units": "real", "k_kcal_mol": [1., 2., 3., 4.]
})
```

Converting this torsion equation does not import or validate an entire OPLS
force field. Units, nonbonded exclusions, 1–4 scaling, mixing rules and parameter
coverage remain separate scientific contracts.

## Verification and scope

`tests/test_torsion_conventions.py` checks independent energy equations at
nontrivial phases/angles, hand-calculated unit/prefactor cases, signed and zero
terms, term-order retention, Fourier round trips and unsupported input rejection.
Equality is floating-point numerical agreement, not bitwise recovery of decimal
source files. No torsion profile is claimed to predict a condensed-phase polymer
population or material property merely because its conversion is correct.

Primary equations checked on 7 November 2026:

- [LAMMPS Fourier dihedral](https://docs.lammps.org/dihedral_fourier.html)
- [LAMMPS OPLS dihedral](https://docs.lammps.org/dihedral_opls.html)
- [LAMMPS units](https://docs.lammps.org/units.html)

The implementation is written for PolyChemMD from the cited functional-form
definitions. No upstream polymer workflow code or parameter database was copied.
