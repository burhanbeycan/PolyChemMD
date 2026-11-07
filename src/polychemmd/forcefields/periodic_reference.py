"""Independent finite Ewald energy/force reference for periodic alkane cells.

This is a small, explicit verification kernel, not an MD integrator or evidence
that a selected finite Ewald parameterization has converged for a bulk material.
"""

from __future__ import annotations

from itertools import combinations, product
import math

from .contracts import ParameterSet
from .reference import (
    COULOMB_KJ_MOL_ANGSTROM_E2 as COULOMB_CONSTANT,
    _angle_energy, _bond_energy, _coordinates, _differentiate, _dihedral_energy,
)


_MAX_ATOMS = 512
_MAX_PAIR_EVALUATIONS = _MAX_ATOMS*(_MAX_ATOMS-1)//2
_MAX_NMAX = 24
_MAX_RECIPROCAL_ATOM_EVALUATIONS = 20_000_000
_SPHERE_SOURCE = (
    "https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/"
    "src/KSPACE/ewald.cpp#L291-L310"
)


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be positive and finite")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be positive and finite") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _minimum_image(displacement, box):
    """Return the nearest orthorhombic image, rejecting unrepresentable math."""
    result = []
    for value, length in zip(displacement, box):
        ratio = value / length
        if not math.isfinite(value) or not math.isfinite(ratio):
            raise ValueError("periodic displacement is numerically unrepresentable")
        result.append(value - length * math.floor(ratio + .5))
    if not all(math.isfinite(value) for value in result):
        raise ValueError("periodic displacement is numerically unrepresentable")
    return tuple(result)


def _pair_budget(atom_count, maximum):
    if type(maximum) is not int or not 1 <= maximum <= _MAX_PAIR_EVALUATIONS:
        raise ValueError(f"max_pair_evaluations must be an integer in [1, {_MAX_PAIR_EVALUATIONS}]")
    count = atom_count*(atom_count-1)//2
    if count > maximum:
        raise ValueError("independent periodic reference pair-evaluation budget exceeded")
    return count


def _numeric_pairs(data):
    """Yield small-kernel numeric pair data, without changing the assignment.

    Unscaled graph distance is immaterial to these equations. The compact
    assignment's complete default applies to all pairs absent from its sparse
    exception table; an intra/inter-molecule flag is retained for diagnostics.
    No global pair array is materialized, even for this bounded reference.
    """
    if data["schema_version"] == "1":
        yield from data["nonbonded_pairs"]
        return
    exceptions = {frozenset(pair["atom_ids"]): pair for pair in data["special_pairs"]}
    for first, second in combinations(data["atoms"], 2):
        atom_ids = [first["cell_atom_id"], second["cell_atom_id"]]
        exception = exceptions.get(frozenset(atom_ids))
        if exception is not None:
            yield exception
        else:
            yield {"atom_ids": atom_ids,
                "same_molecule": first["molecule_index"] == second["molecule_index"],
                "lj_scale": 1.0, "coulomb_scale": 1.0}


def _pair_terms(displacement, charge_product, epsilon, sigma, lj_scale,
                coulomb_scale, alpha):
    """Return LJ, screened real, special correction and combined pair force.

    The supplied displacement points from the second atom to the first. Ewald
    real space is not multiplied by the special scale: its missing full-Coulomb
    part is removed separately, as required with reciprocal electrostatics.
    """
    r = math.hypot(*displacement)
    if not math.isfinite(r) or r <= 0:
        raise ValueError("coincident or unrepresentable periodic pair geometry")
    ar = alpha*r
    if not math.isfinite(ar):
        raise ValueError("Ewald pair argument is numerically unrepresentable")
    try:
        if lj_scale and epsilon:
            sixth = (sigma/r)**6
            lj = 4*lj_scale*epsilon*(sixth*sixth-sixth)
            lj_prefactor = 24*lj_scale*epsilon*(2*sixth*sixth-sixth)/(r*r)
        else:
            lj = lj_prefactor = 0.0
        erfc = math.erfc(ar)
        gaussian = math.exp(-ar*ar)
        cq = COULOMB_CONSTANT*charge_product
        real = cq*erfc/r
        correction = -(1-coulomb_scale)*cq/r
        electrostatic_prefactor = cq*(erfc + 2*ar/math.sqrt(math.pi)*gaussian
                                      - (1-coulomb_scale))/(r*r*r)
        force = tuple((lj_prefactor+electrostatic_prefactor)*v for v in displacement)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("periodic pair energy/force exceeds the numeric range") from exc
    if not all(math.isfinite(v) for v in (lj, real, correction, *force)):
        raise ValueError("periodic pair energy/force exceeds the numeric range")
    return lj, real, correction, force


def _reciprocal_budget(nmax, atom_count, max_atom_evaluations):
    if type(nmax) is not int or not 1 <= nmax <= _MAX_NMAX:
        raise ValueError(f"reciprocal_nmax must be an integer in [1, {_MAX_NMAX}]")
    if type(max_atom_evaluations) is not int or max_atom_evaluations <= 0:
        raise ValueError("reciprocal work budget must be a positive integer")
    mode_count = ((2*nmax+1)**3-1)//2
    if mode_count*atom_count > max_atom_evaluations:
        raise ValueError("reciprocal atom-evaluation budget exceeded")
    return mode_count


def _reciprocal_ewald(points, charges, box, alpha, nmax,
                      max_atom_evaluations=_MAX_RECIPROCAL_ATOM_EVALUATIONS):
    """Explicit cube of Ewald vectors; exploit the exact k/-k symmetry only."""
    mode_count = _reciprocal_budget(nmax, len(points), max_atom_evaluations)
    volume = math.prod(box)
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("periodic cell volume is outside the numeric range")
    unit_k = tuple(2*math.pi/length for length in box)
    if not all(math.isfinite(v) and v > 0 for v in unit_k):
        raise ValueError("reciprocal wavevectors are outside the numeric range")
    # The cited engine's sphere differs from the declared reference cube.
    # Triangle bounds quantify *that finite-list difference*, not the omitted
    # infinite Ewald tail or real-space truncation error.
    try:
        sphere_limit = 1.00001*max((nmax*v)**2 for v in unit_k)
        energy_factor = 4*math.pi*COULOMB_CONSTANT/volume
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("reciprocal scale is outside the numeric range") from exc
    if not math.isfinite(sphere_limit) or not math.isfinite(energy_factor) or energy_factor <= 0:
        raise ValueError("reciprocal scale is outside the numeric range")
    force_factor = 2*energy_factor
    if not math.isfinite(force_factor):
        raise ValueError("reciprocal force scale is outside the numeric range")
    energies, omitted_weights = [], []
    omitted_axis = [[], [], []]
    forces = [[0., 0., 0.] for _ in points]
    wrapped = [tuple(value % length for value, length in zip(p, box)) for p in points]
    for integer in product(range(-nmax, nmax+1), repeat=3):
        # First nonzero component positive selects exactly one of each +/- pair.
        first_nonzero = next((v for v in integer if v), 0)
        if first_nonzero <= 0:
            continue
        k = tuple(n*u for n, u in zip(integer, unit_k))
        try:
            k2 = math.fsum(v*v for v in k)
        except OverflowError as exc:
            raise ValueError("squared reciprocal wavevector exceeds the numeric range") from exc
        if not math.isfinite(k2) or k2 <= 0:
            raise ValueError("squared reciprocal wavevector is outside the positive numeric range")
        scaled = .5*math.sqrt(k2)/alpha
        weight = math.exp(-scaled*scaled)/k2
        if not math.isfinite(weight):
            raise ValueError("reciprocal weight exceeds the finite numeric range")
        if k2 > sphere_limit:
            omitted_weights.append(weight)
            for axis in range(3):
                omitted_axis[axis].append(weight*abs(k[axis]))
        phases = [math.fsum(a*b for a, b in zip(k, point)) for point in wrapped]
        if not all(math.isfinite(phase) for phase in phases):
            raise ValueError("reciprocal phase exceeds the finite numeric range")
        cosines, sines = [math.cos(p) for p in phases], [math.sin(p) for p in phases]
        real = math.fsum(q*c for q, c in zip(charges, cosines))
        imag = math.fsum(q*s for q, s in zip(charges, sines))
        energies.append(energy_factor*weight*(real*real+imag*imag))
        for i, (q, cosine, sine) in enumerate(zip(charges, cosines, sines)):
            coefficient = force_factor*weight*q*(real*sine-imag*cosine)
            for axis in range(3):
                forces[i][axis] += coefficient*k[axis]
    q_abs = math.fsum(abs(q) for q in charges)
    maximum_q = max(abs(q) for q in charges)
    energy = math.fsum(energies)
    bounds = {
        "cited_backend_screen_source": _SPHERE_SOURCE,
        "cited_sphere_squared_wavevector_limit_per_angstrom2": sphere_limit,
        "energy_triangle_bound_kj_mol": energy_factor*q_abs*q_abs*math.fsum(omitted_weights),
        "max_force_component_triangle_bound_kj_mol_angstrom": (
            force_factor*maximum_q*q_abs*max(math.fsum(values) for values in omitted_axis)),
        "scope": "finite reference cube versus cited engine sphere only; not an infinite-sum, real-cutoff or physical-accuracy error bound",
        "floating_point_note": "triangle expressions evaluated in ordinary double precision",
    }
    if not all(math.isfinite(v) for v in (energy, *[v for row in forces for v in row],
                                         bounds["energy_triangle_bound_kj_mol"],
                                         bounds["max_force_component_triangle_bound_kj_mol_angstrom"])):
        raise ValueError("reciprocal energy/force exceeds the numeric range")
    return energy, forces, {
        "shape": "full_integer_cube_with_exact_plus_minus_symmetry",
        "integer_nmax": nmax, "full_vector_count": 2*mode_count,
        "evaluated_symmetry_reduced_modes": mode_count,
        "atom_mode_evaluations": mode_count*len(points),
        "max_atom_mode_evaluations": max_atom_evaluations,
        "cited_backend_finite_list_difference": bounds,
    }


def evaluate_periodic_alkanes(
    assigned, *, cutoff_angstrom=10.0, ewald_alpha_per_angstrom=.30,
    reciprocal_nmax=16, coordinates_angstrom=None,
    max_reciprocal_atom_evaluations=_MAX_RECIPROCAL_ATOM_EVALUATIONS,
    max_pair_evaluations=_MAX_PAIR_EVALUATIONS,
):
    """Evaluate a neutral 3D periodic orthorhombic cell using finite Ewald sums.

    Bonded geometry is molecularly unwrapped. Nonbonded real-space terms use
    minimum images, while reciprocal electrostatics includes the declared cube
    of wavevectors. Tinfoil boundary conditions and zero net charge are required.
    """
    from .periodic import AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell
    if not isinstance(assigned, (AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell)):
        raise ValueError("assigned must be a validated periodic alkane assignment")
    data = assigned.to_dict()
    atoms = data["atoms"]
    if not 0 < len(atoms) <= _MAX_ATOMS:
        raise ValueError(f"periodic reference supports 1 to {_MAX_ATOMS} atoms")
    pair_evaluations = _pair_budget(len(atoms), max_pair_evaluations)
    _reciprocal_budget(reciprocal_nmax, len(atoms), max_reciprocal_atom_evaluations)
    box_data = data["box"]
    if (box_data["periodic"] != [True, True, True]
            or box_data["origin_angstrom"] != [0., 0., 0.]):
        raise ValueError("periodic reference requires three periodic zero-origin axes")
    box = tuple(_positive(v, "box length") for v in box_data["lengths_angstrom"])
    if len(box) != 3:
        raise ValueError("periodic reference requires three box lengths")
    cutoff = _positive(cutoff_angstrom, "cutoff_angstrom")
    alpha = _positive(ewald_alpha_per_angstrom, "ewald_alpha_per_angstrom")
    if cutoff >= min(box)/2:
        raise ValueError("cutoff must be strictly smaller than half every box length")
    volume = math.prod(box)
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("periodic cell volume is outside the numeric range")
    ids = [a["cell_atom_id"] for a in atoms]
    index = {identity: i for i, identity in enumerate(ids)}
    points = _coordinates([a["coordinates_angstrom"] for a in atoms] if coordinates_angstrom is None
                          else coordinates_angstrom, len(atoms))
    charges = tuple(float(a["charge_e"]) for a in atoms)
    if any(not math.isfinite(q) for q in charges) or abs(math.fsum(charges)) > 1e-12:
        raise ValueError("periodic Ewald reference requires finite neutral charges")
    parameters = ParameterSet.from_mapping(data["parameter_set"])
    if (parameters.nonbonded.mixing_rule != "geometric"
            or parameters.nonbonded.lj_14_scale != .5
            or parameters.nonbonded.coulomb_14_scale != .5):
        raise ValueError("periodic reference requires the frozen OPLS mixing/pair conventions")
    if data["interactions"]["impropers"]:
        raise ValueError("impropers are unsupported in the restricted periodic alkane model")
    definitions = {(p.kind, p.type_id): p.values for p in parameters.definitions}
    forces = [[0., 0., 0.] for _ in atoms]
    values = {name: [] for name in ("bond", "angle", "dihedral", "lj", "coulomb_screened", "coulomb_exclusion")}
    # All bond vectors must use their intended molecular image, rather than
    # silently turning a different periodic atom image into a bonded neighbour.
    for bond in data["interactions"]["bonds"]:
        first, second = (points[index[identity]] for identity in bond["atom_ids"])
        displacement = tuple(a-b for a, b in zip(first, second))
        if any(abs(d) >= length/2 for d, length in zip(displacement, box)):
            raise ValueError("bonded unwrapped displacement exceeds its unique periodic image")
    for kind, component, kernel in (("bonds", "bond", _bond_energy),
                                    ("angles", "angle", _angle_energy),
                                    ("dihedrals", "dihedral", _dihedral_energy)):
        for interaction in data["interactions"][kind]:
            indices = [index[identity] for identity in interaction["atom_ids"]]
            parameter = definitions[kind, interaction["type_id"]]
            energy, gradient = _differentiate(tuple(points[i] for i in indices),
                                               lambda p, f=kernel, v=parameter: f(p, v))
            values[component].append(energy)
            for i, derivative in zip(indices, gradient):
                for axis in range(3):
                    forces[i][axis] -= derivative[axis]
    pairs_inside = pairs_outside = special_pairs = intermolecular_inside = 0
    for pair in _numeric_pairs(data):
        first, second = (index[identity] for identity in pair["atom_ids"])
        raw = tuple(a-b for a, b in zip(points[first], points[second]))
        displacement = _minimum_image(raw, box)
        r = math.hypot(*displacement)
        if not math.isfinite(r) or r <= 0:
            raise ValueError("coincident or unrepresentable periodic pair geometry")
        special = pair["coulomb_scale"] != 1.0
        if special:
            if (any(abs(v) >= length/2 for v, length in zip(raw, box)) or r >= cutoff):
                raise ValueError("every excluded/scaled short pair must use its molecular minimum image inside the real cutoff")
            special_pairs += 1
        if abs(r-cutoff) <= max(1e-10, 1e-12*cutoff):
            raise ValueError("periodic pair lies on the hard real-space/LJ cutoff")
        if r >= cutoff:
            pairs_outside += 1
            continue
        left = definitions["atoms", atoms[first]["type_id"]]
        right = definitions["atoms", atoms[second]["type_id"]]
        epsilon = math.sqrt(left["epsilon_kj_mol"]*right["epsilon_kj_mol"])
        sigma = math.sqrt(left["sigma_angstrom"]*right["sigma_angstrom"])
        lj, real, correction, pair_force = _pair_terms(
            displacement, charges[first]*charges[second], epsilon, sigma,
            pair["lj_scale"], pair["coulomb_scale"], alpha)
        values["lj"].append(lj)
        values["coulomb_screened"].append(real)
        values["coulomb_exclusion"].append(correction)
        for axis, force in enumerate(pair_force):
            forces[first][axis] += force
            forces[second][axis] -= force
        pairs_inside += 1
        intermolecular_inside += not pair["same_molecule"]
    reciprocal, reciprocal_forces, truncation = _reciprocal_ewald(
        points, charges, box, alpha, reciprocal_nmax, max_reciprocal_atom_evaluations)
    self_energy = -COULOMB_CONSTANT*alpha/math.sqrt(math.pi)*math.fsum(q*q for q in charges)
    for target, source in zip(forces, reciprocal_forces):
        for axis in range(3):
            target[axis] += source[axis]
    try:
        sums = {name: math.fsum(entries) for name, entries in values.items()}
        components = {name: sums[name] for name in ("bond", "angle", "dihedral", "lj")}
        components["coulomb_realspace"] = sums["coulomb_screened"]+sums["coulomb_exclusion"]
        components["coulomb_kspace"] = reciprocal+self_energy
        components["total"] = math.fsum(components.values())
    except OverflowError as exc:
        raise ValueError("periodic energy exceeds the finite numeric range") from exc
    if any(not math.isfinite(v) for v in components.values()) or any(not math.isfinite(v) for row in forces for v in row):
        raise ValueError("periodic energy/force exceeds the finite numeric range")
    return {
        "schema_version": "1", "kind": "periodic_alkane_ewald_energy_force_reference",
        "model_id": data["model_id"], "atom_ids": ids,
        "assignment_schema_version": data["schema_version"],
        "pair_work": {"pair_evaluations": pair_evaluations,
            "max_pair_evaluations": max_pair_evaluations,
            "global_pair_array_materialized_by_kernel": False},
        "evaluated_coordinates_angstrom": [list(p) for p in points],
        "coordinates_override_used": coordinates_angstrom is not None,
        "energy_components_kj_mol": components, "forces_kj_mol_angstrom": forces,
        "energy_units": "kJ/mol", "force_units": "kJ/(mol angstrom)",
        "boundary_conditions": "three_dimensional_periodic_orthorhombic_tinfoil",
        "box_lengths_angstrom": list(box), "cutoff_angstrom": cutoff,
        "ewald_alpha_per_angstrom": alpha, "reciprocal_nmax": reciprocal_nmax,
        "reciprocal_truncation": truncation,
        "coulomb_decomposition_kj_mol": {
            "screened_real": sums["coulomb_screened"], "special_pair_exclusion": sums["coulomb_exclusion"],
            "reciprocal": reciprocal, "self": self_energy,
            "surface": 0.0, "neutralizing_background": 0.0,
        },
        "coulomb_constant_kj_mol_angstrom_e2": COULOMB_CONSTANT,
        "real_space_erfc": "Python math.erfc; no engine-specific polynomial approximation is substituted",
        "dielectric_relative": 1.0, "lj_mixing_rule": "geometric_sigma_and_epsilon",
        "lj_energy_shift": False, "lj_tail_correction": False,
        "pair_counts": {"inside_real_cutoff": pairs_inside, "outside_real_cutoff": pairs_outside,
                        "special_pairs": special_pairs, "intermolecular_inside_cutoff": intermolecular_inside},
        "force_derivation": "analytical Ewald/LJ pair and reciprocal derivatives; project local analytic differentiation for bonded terms",
        "finite_difference_forces": False, "md_integration_steps": 0,
        "equilibrium_status": "not_assessed", "material_property_validation": "not_performed",
        "convergence_status": "not_established; explicit finite real and reciprocal sums",
    }
