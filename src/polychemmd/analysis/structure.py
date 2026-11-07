"""Static chain dimensions and apolar orientational order for validated alkanes.

These are configuration observables, not equilibrium ensemble averages or phase
classification. The implementation uses unwrapped molecular coordinates and
only the Python standard library.
"""

import hashlib
import json
import math


_MAX_ATOMS = 8192
_COORDINATE_TOLERANCE = 1e-8
_DIRECTOR_GAP_TOLERANCE = 1e-12


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be representable") from exc
    if not math.isfinite(value) or positive and value <= 0:
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


def _vector(value, name):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly three finite numbers")
    return [_number(v, name) for v in value]


def _integer(value, name, *, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _eigensystem_symmetric3(matrix):
    """Bounded independently written Jacobi diagonalization of a 3x3 matrix."""
    if not isinstance(matrix, (list, tuple)) or len(matrix) != 3:
        raise ValueError("symmetric matrix must be 3 by 3")
    original = [_vector(row, "symmetric matrix row") for row in matrix]
    scale = max(abs(v) for row in original for v in row)
    tolerance = 1e-12*max(1., scale)
    if any(abs(original[i][j]-original[j][i]) > tolerance for i in range(3) for j in range(3)):
        raise ValueError("matrix must be symmetric")
    if scale == 0:
        return [0., 0., 0.], [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]], 0
    a = [[v/scale for v in row] for row in original]
    normalized = [row[:] for row in a]
    vectors = [[float(i == j) for j in range(3)] for i in range(3)]
    for iteration in range(65):
        p, q = max(((0, 1), (0, 2), (1, 2)), key=lambda pair: abs(a[pair[0]][pair[1]]))
        if abs(a[p][q]) <= 1e-15:
            break
        if iteration == 64:
            raise ValueError("symmetric eigensolver exhausted its 64-rotation budget")
        tau = (a[q][q]-a[p][p])/(2*a[p][q])
        tangent = math.copysign(1., tau)/(abs(tau)+math.hypot(1., tau))
        cosine = 1/math.sqrt(1+tangent*tangent)
        sine = tangent*cosine
        off = a[p][q]
        a[p][p] -= tangent*off
        a[q][q] += tangent*off
        a[p][q] = a[q][p] = 0.
        for index in range(3):
            if index not in (p, q):
                first, second = a[index][p], a[index][q]
                a[index][p] = a[p][index] = cosine*first-sine*second
                a[index][q] = a[q][index] = sine*first+cosine*second
            first, second = vectors[index][p], vectors[index][q]
            vectors[index][p] = cosine*first-sine*second
            vectors[index][q] = sine*first+cosine*second
    ordered = sorted(range(3), key=lambda i: a[i][i], reverse=True)
    values = [a[i][i]*scale for i in ordered]
    axes = [[vectors[row][i] for row in range(3)] for i in ordered]
    if not all(math.isfinite(v) for v in values):
        raise ValueError("matrix eigenvalues exceed the finite numeric range")
    for value, axis in zip([a[i][i] for i in ordered], axes):
        residual = math.hypot(*(math.fsum(normalized[i][j]*axis[j] for j in range(3))-value*axis[i]
                                for i in range(3)))
        if residual > 1e-12:
            raise ValueError("symmetric eigensolver residual check failed")
    if any(abs(math.fsum(x*y for x, y in zip(first, second))-float(i == j)) > 1e-12
           for i, first in enumerate(axes) for j, second in enumerate(axes)):
        raise ValueError("symmetric eigensolver orthogonality check failed")
    return values, axes, iteration


def _gyration(points, masses):
    """Mass-weighted tensor around a center calculated relative to one anchor."""
    if not points or len(points) != len(masses):
        raise ValueError("gyration requires matching nonempty point and mass arrays")
    points = [_vector(point, "gyration coordinate") for point in points]
    masses = [_number(mass, "atomic molar mass", positive=True) for mass in masses]
    try:
        total = math.fsum(masses)
        weights = [mass/total for mass in masses]
        relative = [[value-anchor for value, anchor in zip(point, points[0])] for point in points]
        if not math.isfinite(total) or any(not math.isfinite(v) for row in relative for v in row):
            raise ValueError("gyration centering exceeds the finite numeric range")
        center = [math.fsum(weight*point[i] for weight, point in zip(weights, relative)) for i in range(3)]
        centered = [[value-offset for value, offset in zip(point, center)] for point in relative]
        tensor = [[math.fsum(weight*point[i]*point[j] for weight, point in zip(weights, centered))
                   for j in range(3)] for i in range(3)]
        if not all(math.isfinite(value) for row in tensor for value in row):
            raise ValueError("gyration tensor exceeds the finite numeric range")
        trace = math.fsum(tensor[i][i] for i in range(3))
        if trace == 0 and any(point != relative[0] for point in relative[1:]):
            raise ValueError("nonzero gyration dimensions underflow the representable tensor range")
        eigenvalues, _, iterations = _eigensystem_symmetric3(tensor)
    except (OverflowError, ZeroDivisionError) as exc:
        raise ValueError("gyration statistics exceed the finite numeric range") from exc
    if not math.isfinite(trace) or trace < 0 or min(eigenvalues) < -1e-12*max(1., trace):
        raise ValueError("gyration tensor is not finite positive semidefinite")
    eigenvalues = [max(0., value) for value in eigenvalues]
    anisotropy = (1.-3.*math.fsum((eigenvalues[i]/trace)*(eigenvalues[j]/trace)
                                 for i, j in ((0, 1), (0, 2), (1, 2)))) if trace > 0 else None
    if anisotropy is not None:
        if not -1e-12 <= anisotropy <= 1+1e-12:
            raise ValueError("relative shape anisotropy lies outside its allowed range")
        anisotropy = min(1., max(0., anisotropy))
    return {"radius_of_gyration_angstrom": math.sqrt(trace),
            "gyration_tensor_angstrom2": tensor, "gyration_eigenvalues_angstrom2": eigenvalues,
            "relative_shape_anisotropy": anisotropy, "eigensolver_rotations": iterations}


def _nematic(vectors):
    if not vectors:
        raise ValueError("nematic order requires at least one chain vector")
    normalized = []
    for vector in vectors:
        vector = _vector(vector, "backbone end-to-end vector")
        norm = math.hypot(*vector)
        if not math.isfinite(norm):
            raise ValueError("backbone end-to-end length is unrepresentable")
        if norm == 0:
            return {"orientation_tensor": None, "nematic_order_parameter": None,
                    "nematic_director": None, "orientation_status": "inconclusive_degenerate_endpoints"}
        scale = max(abs(value) for value in vector)
        scaled = [value/scale for value in vector]
        scaled_norm = math.hypot(*scaled)
        normalized.append([value/scaled_norm for value in scaled])
    count = len(normalized)
    tensor = [[1.5*math.fsum(vector[i]*vector[j]/count for vector in normalized)-.5*float(i == j)
               for j in range(3)] for i in range(3)]
    values, axes, iterations = _eigensystem_symmetric3(tensor)
    order = values[0]
    if not -1e-12 <= order <= 1+1e-12 or abs(math.fsum(values)) > 1e-12:
        raise ValueError("nematic tensor eigenvalues violate the normalized traceless convention")
    order = min(1., max(0., order))
    if values[0]-values[1] <= _DIRECTOR_GAP_TOLERANCE:
        director, status = None, "inconclusive_degenerate_director"
    else:
        director = axes[0]
        # Fix a presentation sign only; n and -n denote the same apolar axis.
        index = max(range(3), key=lambda i: abs(director[i]))
        if director[index] < 0:
            director = [-value for value in director]
        status = "computed"
    return {"orientation_tensor": tensor, "orientation_eigenvalues": values,
            "nematic_order_parameter": order, "nematic_director": director,
            "orientation_status": status, "orientation_eigensolver_rotations": iterations}


def _state_coordinates(record, state):
    atoms = record["atoms"]
    if state is None:
        state = {"step": 0, "atom_count": len(atoms),
                 "box": {name: record["box"][name] for name in ("origin_angstrom", "lengths_angstrom", "periodic")},
                 "atoms": [{"engine_id": i+1, "cell_atom_id": atom["cell_atom_id"],
                            "molecule_index": atom["molecule_index"], "charge_e": atom["charge_e"],
                            **{name: atom[name] for name in ("coordinates_angstrom", "wrapped_coordinates_angstrom", "image_flags")}}
                           for i, atom in enumerate(atoms)]}
        source = "assigned_initial_geometry"
    else:
        source = "supplied_engine_state"
    if not isinstance(state, dict) or not {"step", "atom_count", "box", "atoms"}.issubset(state):
        raise ValueError("state requires a full step/box/atom snapshot, not a public summary")
    _integer(state["step"], "state step")
    if type(state["atom_count"]) is not int or state["atom_count"] != len(atoms):
        raise ValueError("state atom_count differs from the frozen assignment")
    try:
        state_hash = hashlib.sha256(_canonical(state)).hexdigest()
    except (ValueError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError("state must contain finite serializable JSON data") from exc
    box = state["box"]
    if (not isinstance(box, dict) or set(box) != {"origin_angstrom", "lengths_angstrom", "periodic"}
            or not isinstance(box["periodic"], list) or any(v is not True for v in box["periodic"])
            or len(box["periodic"]) != 3):
        raise ValueError("state must retain an orthorhombic fully periodic box")
    origin = _vector(box["origin_angstrom"], "box origin")
    lengths = [_number(value, "box length", positive=True) for value in _vector(box["lengths_angstrom"], "box lengths")]
    bounds = [low+length for low, length in zip(origin, lengths)]
    volume = math.prod(lengths)
    if (not math.isfinite(volume) or volume <= 0 or any(not math.isfinite(high) or high <= low
            for low, high in zip(origin, bounds))):
        raise ValueError("state box volume or bounds are unrepresentable")
    if "volume_angstrom3" in state and not math.isclose(_number(state["volume_angstrom3"], "state volume", positive=True),
                                                       volume, rel_tol=1e-10, abs_tol=1e-8):
        raise ValueError("state volume differs from the actual box")
    rows = state["atoms"]
    if not isinstance(rows, list) or len(rows) != len(atoms):
        raise ValueError("state requires one full atom row per source atom")
    lookup = {atom["cell_atom_id"]: (index+1, atom) for index, atom in enumerate(atoms)}
    coordinates = {}
    required = {"engine_id", "cell_atom_id", "molecule_index", "charge_e", "coordinates_angstrom",
                "wrapped_coordinates_angstrom", "image_flags"}
    for row in rows:
        if not isinstance(row, dict) or not required.issubset(row):
            raise ValueError("state atom rows lack required identity/charge/coordinate/image fields")
        identity = row["cell_atom_id"]
        if not isinstance(identity, str) or identity not in lookup or identity in coordinates:
            raise ValueError("state atom identities must uniquely match the assignment")
        index, atom = lookup[identity]
        if (type(row["engine_id"]) is not int or row["engine_id"] != index
                or type(row["molecule_index"]) is not int or row["molecule_index"] != atom["molecule_index"]):
            raise ValueError("state engine or molecule identity changed")
        if not math.isclose(_number(row["charge_e"], "state partial charge"), atom["charge_e"], rel_tol=0, abs_tol=1e-12):
            raise ValueError("state fixed partial charge changed")
        for name in ("element", "isotope", "formal_charge", "type_id"):
            if name in row and (type(row[name]) is not type(atom[name]) or row[name] != atom[name]):
                raise ValueError(f"state {name} differs from its frozen atom")
        if "mass_g_mol" in row and not math.isclose(_number(row["mass_g_mol"], "state mass", positive=True),
                                                  atom["mass_g_mol"], rel_tol=0, abs_tol=1e-12):
            raise ValueError("state mass differs from the frozen parameter mass")
        unwrapped = _vector(row["coordinates_angstrom"], "state unwrapped coordinates")
        wrapped = _vector(row["wrapped_coordinates_angstrom"], "state wrapped coordinates")
        images = row["image_flags"]
        if not isinstance(images, (list, tuple)) or len(images) != 3 or any(type(v) is not int for v in images):
            raise ValueError("state images must be three integer flags")
        for value, point, image, low, high, length in zip(unwrapped, wrapped, images, origin, bounds, lengths):
            if not low <= point < high:
                raise ValueError("state wrapped coordinates lie outside their half-open box")
            try:
                restored = point+image*length
            except OverflowError as exc:
                raise ValueError("state image reconstruction is unrepresentable") from exc
            if not math.isfinite(restored) or abs(value-restored) > _COORDINATE_TOLERANCE:
                raise ValueError("state wrapped/unwrapped/image coordinates disagree")
        for name in ("forces_kj_mol_angstrom", "velocities_angstrom_fs"):
            if name in row:
                _vector(row[name], name)
        coordinates[identity] = unwrapped
    for bond in record["interactions"]["bonds"]:
        first, second = (coordinates[identity] for identity in bond["atom_ids"])
        delta = [a-b for a, b in zip(first, second)]
        if any(not math.isfinite(value) or abs(value) >= length/2 for value, length in zip(delta, lengths)):
            raise ValueError("state bonded image is inconsistent with the molecular unwrapping convention")
        if math.hypot(*delta) == 0:
            raise ValueError("state contains coincident bonded atoms")
    return coordinates, state_hash, source, state["step"]


def analyze_chain_structure(assigned, state=None):
    """Compute snapshot chain dimensions and carbon-endpoint apolar order.

    Only immutable validated periodic assignments are accepted. An optional
    evolved state must retain source atom identity and molecular image continuity.
    Repeated motifs and molecule counts never become independent-cell counts.
    """
    from ..forcefields.periodic import AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell
    if not isinstance(assigned, (AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell)):
        raise ValueError("assigned must be a validated legacy or compact periodic alkane assignment")
    record = assigned.to_dict()
    atoms = record["atoms"]
    if not 0 < len(atoms) <= _MAX_ATOMS:
        raise ValueError(f"chain structure supports 1 to {_MAX_ATOMS} atoms")
    coordinates, state_hash, source, step = _state_coordinates(record, state)
    groups = {index: [] for index in range(len(record["molecules"]))}
    for atom in atoms:
        groups[atom["molecule_index"]].append(atom)
    carbon_graph = {atom["cell_atom_id"]: set() for atom in atoms if atom["element"] == "C"}
    for bond in record["interactions"]["bonds"]:
        first, second = bond["atom_ids"]
        if first in carbon_graph and second in carbon_graph:
            carbon_graph[first].add(second); carbon_graph[second].add(first)
    per_chain, vectors = [], []
    for index, molecule in groups.items():
        carbons = [atom for atom in molecule if atom["element"] == "C"]
        endpoints = [atom for atom in carbons if len(carbon_graph[atom["cell_atom_id"]]) == 1]
        if len(carbons) < 2 or len(endpoints) != 2 or any(len(carbon_graph[atom["cell_atom_id"]]) not in (1, 2) for atom in carbons):
            raise ValueError("chain analysis requires the supported single-path carbon backbone")
        start = min(endpoints, key=lambda atom: atom["atom_map"])["cell_atom_id"]
        path, previous, current = [], None, start
        while current is not None:
            if current in path:
                raise ValueError("carbon backbone contains a cycle")
            path.append(current)
            neighbors = carbon_graph[current]-({previous} if previous else set())
            if len(neighbors) > 1:
                raise ValueError("carbon backbone contains a branch")
            previous, current = current, next(iter(neighbors)) if neighbors else None
        if set(path) != {atom["cell_atom_id"] for atom in carbons}:
            raise ValueError("carbon backbone is disconnected")
        all_atom = _gyration([coordinates[atom["cell_atom_id"]] for atom in molecule],
                             [atom["mass_g_mol"] for atom in molecule])
        backbone = _gyration([coordinates[atom["cell_atom_id"]] for atom in carbons],
                             [atom["mass_g_mol"] for atom in carbons])
        vector = [b-a for a, b in zip(coordinates[path[0]], coordinates[path[-1]])]
        distance = math.hypot(*vector)
        if not math.isfinite(distance):
            raise ValueError("carbon endpoint separation exceeds the finite numeric range")
        vectors.append(vector)
        per_chain.append({"molecule_index": index, "atom_count": len(molecule), "carbon_count": len(carbons),
            "backbone_atom_ids": path, "terminal_carbon_atom_ids": [path[0], path[-1]],
            "mass_weighted_radius_of_gyration_angstrom": all_atom["radius_of_gyration_angstrom"],
            "mass_weighted_gyration_tensor_angstrom2": all_atom["gyration_tensor_angstrom2"],
            "mass_weighted_gyration_eigenvalues_angstrom2": all_atom["gyration_eigenvalues_angstrom2"],
            "mass_weighted_relative_shape_anisotropy": all_atom["relative_shape_anisotropy"],
            "backbone_radius_of_gyration_angstrom": backbone["radius_of_gyration_angstrom"],
            "backbone_gyration_tensor_angstrom2": backbone["gyration_tensor_angstrom2"],
            "end_to_end_distance_angstrom": distance, "end_to_end_vector_angstrom": vector})
    orientation = _nematic(vectors)
    result = {"schema_version": "1", "kind": "chain_structure_snapshot", "status": "computed",
        "assignment_canonical_sha256": hashlib.sha256(_canonical(record)).hexdigest(),
        "state_canonical_sha256": state_hash, "coordinate_source": source, "step": step,
        "atom_count": len(atoms), "molecule_count": len(groups), "per_chain": per_chain,
        "mean_mass_weighted_radius_of_gyration_angstrom": math.fsum(row["mass_weighted_radius_of_gyration_angstrom"]/len(groups) for row in per_chain),
        "mean_backbone_radius_of_gyration_angstrom": math.fsum(row["backbone_radius_of_gyration_angstrom"]/len(groups) for row in per_chain),
        "mean_end_to_end_distance_angstrom": math.fsum(row["end_to_end_distance_angstrom"]/len(groups) for row in per_chain),
        **orientation, "orientation_vector_definition": "unit vector between terminal backbone carbons, apolar under endpoint reversal",
        "orientation_weighting": "equal weight per molecule; no atom-count or mass weighting of chain axes",
        "dimension_weighting": "all-atom molar masses and carbon-only molar masses; equal-chain means",
        "configuration_statistics_are_equilibrium_averages": False,
        "independent_material_cell_count": 1, "molecules_are_independent_material_cells": False,
        "equilibrium_status": "not_assessed", "phase_classification": "not_performed",
        "production_ready": False, "material_property_validation": "not_performed",
        "uncertainty": "not_estimated_from_one_configuration",
        "units": {"length": "angstrom", "gyration_tensor": "angstrom^2", "nematic_order": "dimensionless"}}
    _canonical(result)
    return result
