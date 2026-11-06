# Assigned-parameter contract and coverage

Status: independently implemented TASK-004 contract. It validates supplied
parameter definitions and checks assigned interaction type IDs. It does not
perform chemical atom typing, fit a force field, generate topology interactions,
assign atomic partial charges, export an MD system or establish physical accuracy.

## API and input

```python
from polychemmd.forcefields import ParameterSet, TopologyRequirements, check_coverage

parameters = ParameterSet.from_mapping(data['parameter_set'])
topology = TopologyRequirements.from_mapping(data['topology'])
report = check_coverage(parameters, topology)
report.require_complete()  # raises MissingParametersError when supplied types lack definitions
output = report.to_dict()
```

`examples/parameter_contract_demo.json` is an intentionally synthetic schema
example. Its arbitrary unit-valued parameters are not a physical atom model or
parameters for any named material. The demonstration topology is a list of type
requirements, not a molecular connectivity graph.

The parameter-set object has exactly `schema_version` (`0.1`), `name`, `units`,
`source`, `nonbonded`, and `parameters`. Required units are energy `kJ/mol`, length
`angstrom`, angle `radian`, mass `g/mol`, and charge `e`. Unit conversion must occur
before import; a mixture of radians/degrees or kcal/kJ fields is rejected.

Each source contains `kind` and `reference`. Supported source kinds are
`analytic_demo`, `literature`, `qm_fit`, and `user_supplied`. Both the parameter set
and every individual definition require provenance. Provenance is recorded, not
independently authenticated by the validator.

## Supported functional forms

All five registry arrays (`atoms`, `bonds`, `angles`, `dihedrals`, `impropers`) are
required, including explicit empty arrays. Every definition contains `type_id`,
`form`, `values`, and `source`. A type ID is a nonempty string, unique within its
interaction kind. The accepted field names make prefactors and units explicit:

| Kind / form | Definition | `values` fields |
| --- | --- | --- |
| atoms / `lj12_6` | `4 epsilon [(sigma/r)^12 - (sigma/r)^6]` | `mass_g_mol`, `epsilon_kj_mol`, `sigma_angstrom` |
| bonds / `harmonic_k_delta_r_squared` | `K (r-r0)^2` | `k_kj_mol_angstrom2`, `r0_angstrom` |
| angles / `harmonic_k_delta_theta_squared` | `K (theta-theta0)^2` | `k_kj_mol_rad2`, `theta0_rad` |
| dihedrals / `fourier_sum_k_1_plus_cos` | `sum K_n [1+cos(n phi-phase)]` | `terms`, with every term's `k_kj_mol`, positive integer `periodicity`, and optional `phase_rad` |
| impropers / `harmonic_k_delta_chi_squared` | `K (chi-chi0)^2` | `k_kj_mol_rad2`, `chi0_rad` |

No implicit factor of one half is introduced. Improper ordering and angular
geometry must be preserved by a future graph/engine adapter; merely calling an
interaction an improper does not determine its correct convention. An engine's
different form, ordering or parameter unit must be converted and tested explicitly.

Mass, sigma, bond/angle/improper force constants and bond equilibrium length must
be positive; epsilon may be zero. Angle equilibrium lies in `(0, pi]`; improper
equilibrium lies in `[-pi, pi]`. Fourier coefficients may be signed. NaN, infinity,
booleans used as numbers and unrepresentable numbers are rejected.

## Nonbonded and 1–4 conventions

The `nonbonded` object requires:

- `mixing_rule`: `lorentz_berthelot` (arithmetic sigma, geometric epsilon) or
  `geometric` (geometric sigma and epsilon), retained without conversion.
- `excluded_bonded_distances`: exactly `[1, 2]`, representing 1–2 and 1–3 exclusions.
- `one_four`: explicit `lj_scale`, `coulomb_scale`, and `handling`.

The handling must be `scaled_pair_once`, with both scales finite in `[0,1]`.
No default force-field family is inferred. Excluding graph distance 3 while also
requesting scaled 1–4 pairs, or specifying overlapping explicit/scaled handling,
is rejected. Other models, explicit pair overrides and many-body force fields
are outside this version's contract. A later adapter must generate special pairs
from molecular connectivity and prove that each pair is counted once.

## Coverage and limitations

The topology object requires all five arrays, each containing exactly `{id,
type_id}` records. Positive integer IDs are normalized to strings; IDs must be
unique within each kind. Atom records must be nonempty. Unassigned or misspelled
types appear in `missing_types` and `missing_interactions`; defined but unused
types appear separately in `unused_types`. Coverage checks angles, dihedrals and
impropers as well as atoms and bonds.

`require_complete()` gates **assigned-type coverage only**. A caller must not
omit interactions from its topology list and then infer a complete physical
model from the report. Connectivity completeness, chemical typing, atomic partial
charges, valence, units at engine boundaries, sampling and external evaluation
need their own evidence. Consequently the report always keeps `production_ready`
false and marks chemical typing and physical accuracy as unassessed.

## Equation and convention sources

Equations are established scientific models, implemented independently here:

- Harmonic bonds: https://docs.lammps.org/bond_harmonic.html
- Harmonic angles: https://docs.lammps.org/angle_harmonic.html
- Fourier torsions: https://docs.lammps.org/dihedral_fourier.html
- Harmonic improper geometry/prefactor: https://docs.lammps.org/improper_harmonic.html
- Special pair weighting and 1–4 caveats: https://docs.lammps.org/special_bonds.html

The recorded definitions are not evidence of a new force field, a force-field fit,
or production readiness for any polymer family.
