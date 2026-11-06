"""Independent isolated-molecule OPLS energy and analytical force reference.

Project-written local forward-mode differentiation evaluates derivatives of the
explicit potential, independently of external engines. It is intentionally a
small verification kernel, not an MD integrator or a bulk electrostatics model.
"""

from __future__ import annotations

from dataclasses import dataclass
import math

from ..geometry import MolecularGeometry
from .contracts import ParameterSet


# The exact numerical conversion used by LAMMPS ``units real`` for engine
# convention parity, not a claim of an exact modern SI-derived constant.
COULOMB_KCAL_MOL_ANGSTROM_E2 = 332.06371
COULOMB_KJ_MOL_ANGSTROM_E2 = COULOMB_KCAL_MOL_ANGSTROM_E2 * 4.184
_SINGULAR_SINE_LIMIT = 1e-12


@dataclass(frozen=True)
class _Dual:
    value: float
    gradient: tuple[float, ...]

    def _constant(self, value):
        if isinstance(value, _Dual):
            if len(value.gradient) != len(self.gradient):
                raise ValueError("incompatible differentiation dimensions")
            return value
        return _Dual(float(value), (0.0,) * len(self.gradient))

    def __add__(self, other):
        other = self._constant(other)
        return _Dual(self.value + other.value,
                     tuple(a + b for a, b in zip(self.gradient, other.gradient)))

    __radd__ = __add__

    def __neg__(self):
        return _Dual(-self.value, tuple(-g for g in self.gradient))

    def __sub__(self, other):
        return self + (-self._constant(other))

    def __rsub__(self, other):
        return self._constant(other) - self

    def __mul__(self, other):
        other = self._constant(other)
        return _Dual(self.value * other.value,
                     tuple(a * other.value + self.value * b
                           for a, b in zip(self.gradient, other.gradient)))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self._constant(other)
        if other.value == 0:
            raise ValueError("singular zero denominator in force evaluation")
        value = self.value / other.value
        return _Dual(value, tuple((a - value * b) / other.value
                                  for a, b in zip(self.gradient, other.gradient)))

    def __rtruediv__(self, other):
        return self._constant(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int):
            raise ValueError("reference differentiation supports integer powers")
        if exponent == 0:
            return self._constant(1.0)
        value = self.value ** exponent
        slope = exponent * self.value ** (exponent - 1)
        return _Dual(value, tuple(slope * g for g in self.gradient))


def _sqrt(value):
    if not math.isfinite(value.value) or value.value <= 0:
        raise ValueError("zero, nonfinite or degenerate geometric norm")
    root = math.sqrt(value.value)
    return _Dual(root, tuple(g / (2 * root) for g in value.gradient))


def _cos(value):
    slope = -math.sin(value.value)
    return _Dual(math.cos(value.value), tuple(slope * g for g in value.gradient))


def _atan2(y, x):
    denominator = x.value*x.value + y.value*y.value
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("singular or unrepresentable angular geometry")
    return _Dual(math.atan2(y.value, x.value),
                 tuple((x.value * dy - y.value * dx) / denominator
                       for dy, dx in zip(y.gradient, x.gradient)))


def _subtract(first, second):
    return tuple(a - b for a, b in zip(first, second))


def _dot(first, second):
    return sum(a*b for a, b in zip(first, second))


def _cross(a, b):
    return (a[1]*b[2] - a[2]*b[1], a[2]*b[0] - a[0]*b[2],
            a[0]*b[1] - a[1]*b[0])


def _norm(vector):
    return _sqrt(_dot(vector, vector))


def _angle(points):
    left, right = _subtract(points[0], points[1]), _subtract(points[2], points[1])
    lengths = _norm(left) * _norm(right)
    cross_length = _norm(_cross(left, right))
    if cross_length.value / lengths.value <= _SINGULAR_SINE_LIMIT:
        raise ValueError("collinear or nearly collinear angle has undefined reference force")
    return _atan2(cross_length, _dot(left, right))


def _dihedral(points):
    first = _subtract(points[1], points[0])
    middle = _subtract(points[2], points[1])
    last = _subtract(points[3], points[2])
    n1, n2 = _cross(first, middle), _cross(middle, last)
    n1_length, n2_length = _norm(n1), _norm(n2)
    first_length, middle_length, last_length = _norm(first), _norm(middle), _norm(last)
    if (n1_length.value / (first_length.value * middle_length.value) <= _SINGULAR_SINE_LIMIT
            or n2_length.value / (middle_length.value * last_length.value) <= _SINGULAR_SINE_LIMIT):
        raise ValueError("collinear or nearly collinear torsion has undefined reference force")
    # Adjacent oriented plane normals give cis=0, trans=+-pi. The reference
    # retains the signed angle, including general explicit Fourier phases.
    y = _dot(_cross(n1, n2), tuple(v / middle_length for v in middle))
    return _atan2(y, _dot(n1, n2))


def _bond_energy(points, values):
    r = _norm(_subtract(points[0], points[1]))
    return values["k_kj_mol_angstrom2"] * (r - values["r0_angstrom"])**2


def _angle_energy(points, values):
    return values["k_kj_mol_rad2"] * (_angle(points) - values["theta0_rad"])**2


def _dihedral_energy(points, values):
    phi = _dihedral(points)
    return sum(term["k_kj_mol"] * (1 + _cos(term["periodicity"] * phi - term["phase_rad"]))
               for term in values["terms"])


def _lj_energy(points, epsilon, sigma, scale=1.0):
    distance = _norm(_subtract(points[0], points[1]))
    sixth = (sigma / distance)**6
    return 4 * scale * epsilon * (sixth**2 - sixth)


def _coulomb_energy(points, charge_product, scale=1.0):
    distance = _norm(_subtract(points[0], points[1]))
    return scale * COULOMB_KJ_MOL_ANGSTROM_E2 * charge_product / distance


def _differentiate(points, function):
    """Return value and exact chain-rule derivatives of one local interaction."""
    dimensions = 3 * len(points)
    variables = tuple(tuple(_Dual(value, tuple(float(k == 3*i + j) for k in range(dimensions)))
                            for j, value in enumerate(point))
                      for i, point in enumerate(points))
    try:
        evaluated = function(variables)
        if not isinstance(evaluated, _Dual):
            raise ValueError("empty or unsupported differentiated interaction")
        if not math.isfinite(evaluated.value) or any(not math.isfinite(g) for g in evaluated.gradient):
            raise ValueError("energy or force exceeds the finite numeric range")
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("energy or force is singular or numerically unrepresentable") from exc
    return evaluated.value, tuple(tuple(evaluated.gradient[3*i + j] for j in range(3))
                                   for i in range(len(points)))


def _coordinates(value, size):
    if not isinstance(value, (tuple, list)) or len(value) != size:
        raise ValueError("coordinates must match the assigned atom sequence")
    coordinates = []
    for point in value:
        if not isinstance(point, (tuple, list)) or len(point) != 3:
            raise ValueError("each coordinate must contain three components")
        row = []
        for item in point:
            if isinstance(item, bool) or not isinstance(item, (int, float)):
                raise ValueError("coordinates must be finite numbers, not booleans")
            try:
                number = float(item)
            except (ValueError, OverflowError) as exc:
                raise ValueError("coordinate is outside the finite numeric range") from exc
            if not math.isfinite(number):
                raise ValueError("coordinates must be finite")
            row.append(number)
        coordinates.append(tuple(row))
    return tuple(coordinates)


def evaluate_alkane(assigned, *, coordinates_angstrom=None, cutoff_angstrom=12.0):
    """Evaluate a restricted assignment in isolation, with no periodic images.

    Nonbonded terms are unshifted and included only at ``r < cutoff``. A point
    on the hard cutoff is rejected because the potential is discontinuous there.
    Reported forces are local analytic chain-rule derivatives, never finite
    differences. Optional positions retain the assigned atom ordering.
    """
    from .alkanes import AssignedAlkane
    if not isinstance(assigned, AssignedAlkane):
        raise ValueError("assigned must be an AssignedAlkane produced by chemical assignment")
    if isinstance(cutoff_angstrom, bool) or not isinstance(cutoff_angstrom, (int, float)):
        raise ValueError("cutoff_angstrom must be positive and finite")
    try:
        cutoff = float(cutoff_angstrom)
    except (ValueError, OverflowError) as exc:
        raise ValueError("cutoff_angstrom must be positive and finite") from exc
    if not math.isfinite(cutoff) or cutoff <= 0:
        raise ValueError("cutoff_angstrom must be positive and finite")
    data = assigned.to_dict()
    geometry = MolecularGeometry.from_mapping(data["geometry"])
    parameters = ParameterSet.from_mapping(data["parameter_set"])
    if parameters.nonbonded.mixing_rule != "geometric":
        raise ValueError("the restricted reference requires geometric mixing")
    if data["interactions"]["impropers"]:
        raise ValueError("impropers are outside the restricted n-alkane reference")
    atom_ids = [a.atom_id for a in geometry.atoms]
    index = {atom_id: i for i, atom_id in enumerate(atom_ids)}
    points = _coordinates(geometry.coordinates_angstrom if coordinates_angstrom is None
                          else coordinates_angstrom, len(atom_ids))
    definitions = {(p.kind, p.type_id): p.values for p in parameters.definitions}
    assignments = {a["atom_id"]: a for a in data["atom_assignments"]}
    forces = [[0.0, 0.0, 0.0] for _ in atom_ids]
    energies = {name: [] for name in ("bond", "angle", "dihedral", "lj", "coulomb")}

    def accumulate(component, ids, function):
        indices = [index[atom_id] for atom_id in ids]
        energy, gradient = _differentiate(tuple(points[i] for i in indices), function)
        energies[component].append(energy)
        for i, local in zip(indices, gradient):
            for axis, derivative in enumerate(local):
                forces[i][axis] -= derivative

    for kind, component, kernel in (("bonds", "bond", _bond_energy),
                                    ("angles", "angle", _angle_energy),
                                    ("dihedrals", "dihedral", _dihedral_energy)):
        for interaction in data["interactions"][kind]:
            values = definitions[kind, interaction["type_id"]]
            accumulate(component, interaction["atom_ids"], lambda p, v=values, f=kernel: f(p, v))

    included_pairs = excluded_pairs = beyond_cutoff_pairs = 0
    for pair in data["nonbonded_pairs"]:
        lj_scale, coulomb_scale = pair["lj_scale"], pair["coulomb_scale"]
        if lj_scale == 0 and coulomb_scale == 0:
            excluded_pairs += 1
            continue
        ids = pair["atom_ids"]
        distance = math.dist(*(points[index[identity]] for identity in ids))
        if not math.isfinite(distance) or distance <= 0:
            raise ValueError("nonbonded pair has coincident or unrepresentable positions")
        if abs(distance - cutoff) <= max(1e-10, 1e-12 * cutoff):
            raise ValueError("nonbonded distance lies on the discontinuous hard cutoff")
        if distance >= cutoff:
            beyond_cutoff_pairs += 1
            continue
        first, second = (assignments[identity] for identity in ids)
        left, right = (definitions["atoms", a["type_id"]] for a in (first, second))
        epsilon = math.sqrt(left["epsilon_kj_mol"] * right["epsilon_kj_mol"])
        sigma = math.sqrt(left["sigma_angstrom"] * right["sigma_angstrom"])
        if lj_scale:
            accumulate("lj", ids, lambda p, e=epsilon, s=sigma, a=lj_scale: _lj_energy(p, e, s, a))
        if coulomb_scale:
            product = first["charge_e"] * second["charge_e"]
            accumulate("coulomb", ids, lambda p, q=product, a=coulomb_scale: _coulomb_energy(p, q, a))
        included_pairs += 1
    try:
        components = {name: math.fsum(values) for name, values in energies.items()}
        components["total"] = math.fsum(components.values())
    except OverflowError as exc:
        raise ValueError("total energy exceeds the finite numeric range") from exc
    if any(not math.isfinite(v) for v in components.values()) or any(not math.isfinite(v) for row in forces for v in row):
        raise ValueError("total energy or force exceeds the finite numeric range")
    return {
        "schema_version": "1", "kind": "isolated_alkane_energy_force_reference",
        "model_id": data["model_id"], "atom_ids": atom_ids,
        "evaluated_coordinates_angstrom": [list(point) for point in points],
        "coordinates_override_used": coordinates_angstrom is not None,
        "energy_components_kj_mol": components,
        "forces_kj_mol_angstrom": forces,
        "energy_units": "kJ/mol", "force_units": "kJ/(mol angstrom)",
        "boundary_conditions": "isolated_nonperiodic",
        "nonbonded_model": "unshifted_lj12_6_plus_truncated_coulomb",
        "cutoff_angstrom": cutoff, "dielectric_relative": 1.0,
        "coulomb_constant_kj_mol_angstrom_e2": COULOMB_KJ_MOL_ANGSTROM_E2,
        "coulomb_constant_basis": "332.06371 kcal angstrom/(mol e^2), matching LAMMPS real-unit numerical convention",
        "lj_mixing_rule": "geometric_sigma_and_epsilon", "tail_correction": False,
        "long_range_electrostatics": "not_used; finite isolated verification model",
        "pair_counts": {"included": included_pairs, "excluded": excluded_pairs,
                        "beyond_cutoff": beyond_cutoff_pairs},
        "force_derivation": "original local forward-mode differentiation of explicit energy equations",
        "finite_difference_forces": False, "trajectory_integration_performed": False,
        "equilibrium_status": "not_assessed", "material_validation": "not_performed",
    }
