# Seeded molecular coordinates

TASK-006A adds one conformer for the supported finite linear-chain graph. The
PolyChemMD graph/lineage/checking layer is independently written; distance
geometry and experimental torsion preferences are supplied by optional RDKit.
No other polymer workflow package is imported or copied.

```python
from polychemmd.chemistry import RepeatUnitSpec
from polychemmd.builders.linear import build_linear
from polychemmd.builders.conformers import generate_conformer

unit = RepeatUnitSpec("PE repeat graph", "[*:1]CC[*:2]", 1, 2)
chain = build_linear(unit, 5)
geometry = generate_conformer(chain, seed=173, max_iterations=1000, max_atoms=2000)
geometry.to_dict()
```

The immutable `MolecularGeometry` contract records explicit atoms, complete
source lineage, map numbers, bonds, isotopes, formal charges, isotope-aware
atomic molar masses, Angstrom coordinates and a SHA-256 of the source graph's
canonical JSON. Coordinates retain the same order as the chain atoms. Loading
a geometry contract validates its structure; it does not by itself independently
reproduce its chemical or geometric checks.

## Construction and checks

- RDKit `ETKDGv3()` with the recorded integer seed, one thread,
  `enforceChirality=True`, explicit positive `maxIterations`, and a 30-second
  backend timeout. The timeout is a backend control, not a hard subprocess
  preemption guarantee. Failure produces an exception and no replacement
  coordinates. Seed values are integers from 0 through 2147483647.
- One realization only: no separate post-embedding force-field minimization,
  energy ranking, conformational search, production force-field assignment or
  molecular dynamics. ETKDG's internal coordinate-objective optimization remains
  part of the embedding method and is not a production MD relaxation.
- Exact mapped-graph identity is checked before/after embedding. All atoms and
  complete atomic lineage are retained; the source graph remains unchanged.
- Coordinates must be finite and marked 3D. No two atoms may be closer than
  0.35 Angstrom. Each bond must lie between 0.55 and 1.65 times the sum of its
  RDKit covalent radii. These intentionally broad rejection bounds detect
  collapsed or detached geometry; they are not a force field or equilibrium
  bond-length validation.
- For specified tetrahedral/E/Z stereo, a separate copy has its existing stereo
  tags and directions removed and reconstructed from the coordinates. Atom maps
  are cleared on both comparison copies, so maps do not distinguish otherwise
  equivalent explicit hydrogens. A mismatch is rejected even when the original
  copied graph tags still agree. Unsupported non-tetrahedral/atropisomeric
  assignments fail explicitly.
- Unspecified chemical stereochemistry stays unspecified in the source graph.
  A spatial realization does not establish tacticity, a stereosequence or a
  stereochemical sampling distribution.

The selected algorithm is explicitly **ETKDGv3**. It uses `ETversion=2` for its
standard torsion definitions, as required by RDKit's documented convention; the
names do not mean that ETversion should be set to 3. Metadata records the actual
RDKit version and embedding switches. Exact seeded reproduction is expected
only within compatible backend/runtime builds, not promised across versions.

## Scientific limits

Generation of a finite conformer does not establish a bulk polymer conformation,
an amorphous phase, relaxed packing, a thermodynamic ensemble, equilibrium,
density or a physical-property prediction. Long or difficult chains can fail
within the configured budget. Packed starting cells and production parameters
remain distinct preparation/validation stages.

Tests include hand-checked ethane dimensions and mass, graph/isotope/charge and
lineage preservation, same-seed reproduction, distinct-seed realization,
coordinate-derived stereo and mirrored/inverted negative cases, malformed
geometry, backend failures and resource/input bounds. Mock failures test only
error handling; passing coordinate tests actually invoke RDKit.

## Primary method and API sources

- RDKit distance geometry and embedding parameters:
  https://www.rdkit.org/docs/source/rdkit.Chem.rdDistGeom.html
- Coordinate-derived stereo assignment:
  https://www.rdkit.org/docs/source/rdkit.Chem.rdmolops.html
- RDKit conformer-generation conventions:
  https://www.rdkit.org/docs/RDKit_Book.html#conformer-generation
- ETKDG method: Riniker and Landrum, *J. Chem. Inf. Model.* **2015**, 55,
  2562–2574, https://doi.org/10.1021/acs.jcim.5b00654 .
- Small-ring/macrocycle improvements: Wang, Witek, Landrum and Riniker,
  *J. Chem. Inf. Model.* **2020**, 60, 2044–2058,
  https://doi.org/10.1021/acs.jcim.0c00025 .

These sources describe the backend method and API, not a novel conformer
generation algorithm or evidence of PolyChemMD material-property accuracy.
