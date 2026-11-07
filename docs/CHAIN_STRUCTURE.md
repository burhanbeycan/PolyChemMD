# Static chain dimensions and apolar order

`analysis/structure.py` independently calculates configuration observables for
the validated neutral explicit n-alkane model. It accepts either periodic
assignment representation and requires no RDKit, NumPy or simulation engine.

```python
from polychemmd.analysis.structure import analyze_chain_structure

initial = analyze_chain_structure(assigned)
evolved = analyze_chain_structure(assigned, preparation_record["final_state"])
```

The original assigned source is immutable. An evolved state retains its own
step, coordinates and periodic box. The analyzer checks every source atom ID,
engine sequence ID, molecule ID and fixed charge; supplied masses/types must
match the frozen model. Coordinates must be finite and satisfy their wrapped
coordinate/image-flag identity. Each unwrapped bond vector must use the unique
molecular image within half each box axis. A whole-molecule lattice translation
is allowed; an inconsistent individual atom image is rejected. Optional forces
and velocities are checked for finite three-component values even though these
structural formulas do not consume them. A public summary that omits atoms
cannot substitute for a full state.

## Chain dimensions and shape

The carbon-carbon bond graph must form the supported single, unbranched path.
Terminal carbons and backbone order are derived from connectivity, rather than
atom names. Terminal order is selected reproducibly by the original atom map.
For each molecule, with supplied atomic molar masses `m_i`, the mass-weighted
gyration tensor and radius are

\[
G_{\alpha\beta}=\frac{1}{M}\sum_i m_i
(r_{i,\alpha}-r_{\rm cm,\alpha})(r_{i,\beta}-r_{\rm cm,\beta}),
\qquad R_g=\sqrt{\operatorname{tr}G}.
\]

Calculations center coordinates relative to a molecular anchor before summing,
reducing cancellation under translation. Tensor entries are in angstrom squared;
off-diagonal entries may be negative and are not square-rooted. The same formula
is calculated separately over all atoms and over carbon-backbone atoms. Carbon
weights use their frozen molar masses, so all ordinary carbons have equal weight.
The endpoint distance `Re` is the Euclidean separation between terminal carbon
positions in molecularly unwrapped coordinates, not their minimum-image shortcut.

For nonnegative gyration eigenvalues `lambda_i`, relative shape anisotropy is

\[
\kappa^2=1-3\frac{\lambda_1\lambda_2+\lambda_1\lambda_3+
\lambda_2\lambda_3}{(\lambda_1+\lambda_2+\lambda_3)^2}.
\]

It is zero for an isotropic tensor and one for a linear point distribution.
The output retains all-atom eigenvalues and relative anisotropy, plus all-atom
and backbone tensors, radii and terminal-carbon vectors. Molecular averages
give equal weight to each chain; they are not a mass-weighted whole-cell radius.

The mass-weighted and unwrapped-coordinate conventions are documented by the
[primary LAMMPS gyration reference](https://docs.lammps.org/compute_gyration.html).
The equivalent gyration-eigenvalue anisotropy convention is given in its
[shape reference](https://docs.lammps.org/compute_gyration_shape.html), which
attributes polymer shape analysis to Theodorou and Suter, *Macromolecules* 18,
1206 (1985). These equations are implemented independently; no LAMMPS shape
routine or another polymer workflow is imported.

## Nematic tensor and director

Each chain contributes its carbon-endpoint unit vector `u_m` with equal weight.
This implementation uses the normalized traceless convention

\[
Q=\frac1{N_{\rm chains}}\sum_m\frac{3u_m u_m^{\mathsf T}-I}{2},
\qquad S=\lambda_{\max}(Q).
\]

Thus perfectly parallel or antiparallel chains give `S=1`, while three equally
weighted orthogonal axes give `Q=0` and `S=0`. Reversing a chain endpoint order
changes its vector sign but not its tensor contribution. The principal axis is
apolar: `n` and `-n` are equivalent. A deterministic presentation sign does not
assign chemical polarity. The
[primary freud nematic-method tutorial](https://freud.readthedocs.io/en/latest/gettingstarted/examples/module_intros/order.Nematic.html)
describes the orientation outer-product/eigenvalue method. Its intermediate
tensor `M-I/3` differs by a factor of `3/2` from the normalized tensor above;
the convention is stated explicitly to avoid silently mixing these matrices.
No freud package or implementation is a dependency.

An independently written symmetric 3x3 Jacobi solver uses normalized matrices,
at most 64 rotations, and residual/orthogonality checks. Degenerate largest
eigenvalues within `1e-12` leave the director null and report
`inconclusive_degenerate_director`; the scalar order remains available.
Consequently an isotropic or exactly planar-axis-degenerate configuration does
not receive an arbitrary meaningful director. If any chain's terminal-carbon
distance is zero, the entire orientation summary is inconclusive: that chain
is not silently excluded and no unit vector is fabricated. Dimensions for all
chains remain available.

## Evidence boundary

Both original assignment and supplied/synthesized state have canonical SHA256
identities in the result. Snapshot statistics cannot establish liquid,
crystalline or amorphous phase, equilibrium chain relaxation, or material
predictivity. Finite-cell `S` also fluctuates and is not a universal phase
classifier. A readiness protocol may use a predeclared orientational-bias
threshold as one diagnostic, with separate state/sampling requirements.

Repeated motifs are **one constructed cell**, not independent replicates.
All output keeps `production_ready=false`, `equilibrium_status=not_assessed`,
and `material_property_validation=not_performed`. No uncertainty is inferred
from one configuration or by counting its molecules as independent samples.

Tests include hand-computed mass-weighted/discrete-rod tensors, signed cross
terms, isotropic shape, orthogonal and antiparallel axes, endpoint degeneracy,
known eigenvalues, legacy/compact agreement, rigid motions, boundary crossings,
identity/image rejection and an actual stored ordered decane motif (`S=1`).
They verify formulas and representation handling, not an equilibrated polymer.
