# Executed coordinate and starting-cell evidence

These files are actual outputs of:

```bash
python scripts/verify_geometry_v03.py --output-dir artifacts/geometry-v03-review
```

They were copied byte-for-byte after the auditor passed, preserving its artifact
hashes. The report records UTC time, software, input hashes, the auditor source
hash, numerical tolerances, measured geometry checks and explicit limitations.
Use a new output directory when reproducing this run.

| File | Actual content |
| --- | --- |
| `conformer.json` | One ETKDGv3 C10H22 realization with atom/bond identity and lineage |
| `cell.json` | Four separately seeded C10H22 starting molecules, 128 explicit atoms, 24 angstrom cube; complete source geometries and image flags |
| `cell.xyz` | Unwrapped Extended XYZ view of `cell.json`; atom order matches the JSON |
| `pilot-pe-dp20.json` | H-capped PE repeat graph, C40H82, 122 atoms |
| `pilot-pvdf-dp8.json` | H-capped PVDF repeat graph, C16H18F16, 50 atoms |
| `pilot-ps-dp4.json` | H-capped PS repeat graph, C32H34, 66 atoms; tacticity unspecified |
| `report.json` | Executed CLI steps, independent image/rigidity/identity audit and pilot summaries |

The nominal cell construction density is about 0.06837 g/cm^3 because that mass
was placed in that prescribed volume. It is **not** a density prediction. The
1.5 angstrom exclusion is geometric, not a force-field contact model. No phase,
equilibrium, production parameterization, energy-ranked conformational ensemble,
MD/QM calculation or material-property accuracy is validated by these records.

Byte hashes of envelopes change when UTC timestamps change. Canonical payload
hashes are also recorded so equivalent geometry payloads can be compared within
the same compatible software/runtime. Exact coordinates across RDKit/Python or
platform changes are not guaranteed. XYZ does not carry bonds, isotope masses,
atom identity or the full provenance; retain its `cell.json` companion.
