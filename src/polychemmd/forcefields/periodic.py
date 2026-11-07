"""Identity-preserving assembly of sourced n-alkane assignments in a cell.

The input is a serialized rigid packed starting cell, not a trajectory frame.
Coordinates/images/rigid transforms and the molecular graph are checked anew;
caller-supplied validation flags are retained only as source metadata.
"""

from dataclasses import InitVar, dataclass
import hashlib
from itertools import combinations
import json
import math

from ..analysis.density import density_from_mass_volume
from ..geometry import MolecularGeometry
from .alkanes import MODEL_ID, assign_alkane


_CELL_KIND = "orthorhombic_rigid_packed_cell"
_ASSIGNED_KIND = "assigned_periodic_opls_aa_n_alkane_cell"
_COMPACT_KIND = "compact_assigned_periodic_opls_aa_n_alkane_cell"
_COMPACT_MAX_ATOMS = 8192
_COMPACT_MAX_MOLECULE_ATOMS = 512
_COMPACT_MAX_INTRAMOLECULAR_PAIRS = 2_000_000
_COMPACT_RULES = {
    "representation": "sparse_topological_exceptions",
    "default": {"lj_scale": 1.0, "coulomb_scale": 1.0},
    "excluded_bonded_distances": [1, 2],
    "one_four": {"lj_scale": .5, "coulomb_scale": .5, "handling": "scaled_pair_once"},
}
_COORDINATE_TOLERANCE = 1e-8
_MATRIX_TOLERANCE = 1e-10
_ATOM_FIELDS = {"cell_atom_id", "atom_map", "source_atom_id", "source_atom_map",
                "molecule_index", "element", "isotope", "formal_charge", "mass_g_mol",
                "lineage", "coordinates_angstrom", "wrapped_coordinates_angstrom", "image_flags"}


def _canonical(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _load_json(raw):
    if not isinstance(raw, str):
        raise ValueError("assigned periodic cell must be JSON object text")
    def constant(value):
        raise ValueError("nonfinite JSON values are unsupported")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON fields are unsupported")
            result[key] = value
        return result
    data = json.loads(raw, parse_constant=constant, object_pairs_hook=unique)
    _canonical(data)
    return data


def _integer(value, name, minimum=0):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _number(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be finite numeric data")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is not representable") from exc
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


def _triple(value, name):
    if not isinstance(value, (tuple, list)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly three components")
    return tuple(_number(item, name) for item in value)


def _close(actual, expected, name, tolerance=_COORDINATE_TOLERANCE):
    if not math.isfinite(actual) or not math.isfinite(expected) or abs(actual-expected) > tolerance:
        raise ValueError(f"{name} does not preserve the source cell identity")


def _determinant(matrix):
    a, b, c = matrix
    return (a[0]*(b[1]*c[2]-b[2]*c[1]) - a[1]*(b[0]*c[2]-b[2]*c[0])
            + a[2]*(b[0]*c[1]-b[1]*c[0]))


@dataclass(frozen=True)
class AssignedPeriodicAlkaneCell:
    """Immutable assembly with a strict frozen-model serialized boundary."""

    _payload_json: str
    max_atoms: InitVar[int] = 512

    def __post_init__(self, max_atoms):
        try:
            data = _load_json(self._payload_json)
            if not isinstance(data, dict) or "source_cell" not in data:
                raise ValueError("the serialized assignment requires its source_cell")
            expected = _assignment_payload(data["source_cell"], max_atoms=max_atoms)
            actual_json, expected_json = _canonical(data), _canonical(expected)
        except (ValueError, TypeError, KeyError, OverflowError, RecursionError) as exc:
            raise ValueError(f"invalid assigned periodic n-alkane cell: {exc}") from exc
        if actual_json != expected_json:
            raise ValueError("assigned periodic cell differs from its validated frozen model/topology/source")
        object.__setattr__(self, "_payload_json", expected_json)

    @classmethod
    def _from_validated(cls, payload):
        instance = object.__new__(cls)
        object.__setattr__(instance, "_payload_json", _canonical(payload))
        return instance

    def to_dict(self):
        return json.loads(self._payload_json)


@dataclass(frozen=True)
class CompactAssignedPeriodicAlkaneCell:
    """Strict frozen assignment with implicit unit-scale nonbonded defaults.

    Unlike schema 1, this representation never materializes global unordered
    pairs. All supported chemistry, rigid source and parameter checks are still
    performed. Local molecule assignment and rigidity work remain quadratic in
    molecule size, with explicit per-molecule and aggregate work limits.
    """

    _payload_json: str
    max_atoms: InitVar[int] = _COMPACT_MAX_ATOMS
    max_atoms_per_molecule: InitVar[int] = _COMPACT_MAX_MOLECULE_ATOMS
    max_intramolecular_pairs: InitVar[int] = _COMPACT_MAX_INTRAMOLECULAR_PAIRS

    def __post_init__(self, max_atoms, max_atoms_per_molecule, max_intramolecular_pairs):
        try:
            data = _load_json(self._payload_json)
            if not isinstance(data, dict) or "source_cell" not in data:
                raise ValueError("the serialized compact assignment requires its source_cell")
            expected = _assignment_payload(data["source_cell"], max_atoms=max_atoms,
                compact=True, max_atoms_per_molecule=max_atoms_per_molecule,
                max_intramolecular_pairs=max_intramolecular_pairs)
            actual_json, expected_json = _canonical(data), _canonical(expected)
        except (ValueError, TypeError, KeyError, OverflowError, RecursionError) as exc:
            raise ValueError(f"invalid compact assigned periodic n-alkane cell: {exc}") from exc
        if actual_json != expected_json:
            raise ValueError("compact assignment differs from its validated frozen model/topology/source")
        object.__setattr__(self, "_payload_json", expected_json)

    @classmethod
    def _from_validated(cls, payload):
        instance = object.__new__(cls)
        object.__setattr__(instance, "_payload_json", _canonical(payload))
        return instance

    def to_dict(self):
        return json.loads(self._payload_json)


def _compact_preflight(cell_record, max_atoms, max_atoms_per_molecule,
                       max_intramolecular_pairs):
    for value, name, ceiling in (
        (max_atoms, "max_atoms", _COMPACT_MAX_ATOMS),
        (max_atoms_per_molecule, "max_atoms_per_molecule", _COMPACT_MAX_MOLECULE_ATOMS),
        (max_intramolecular_pairs, "max_intramolecular_pairs", _COMPACT_MAX_INTRAMOLECULAR_PAIRS),
    ):
        _integer(value, name, 1)
        if value > ceiling:
            raise ValueError(f"{name} cannot exceed this compact implementation's limit {ceiling}")
    molecules = cell_record.get("molecules")
    if not isinstance(molecules, list) or not molecules or len(molecules) > max_atoms:
        raise ValueError("source molecule array is outside the compact work budget")
    work, largest = 0, 0
    for molecule in molecules:
        if not isinstance(molecule, dict) or not isinstance(molecule.get("source_geometry"), dict):
            raise ValueError("each source molecule requires its actual geometry")
        atoms = molecule["source_geometry"].get("atoms")
        if not isinstance(atoms, list) or not 1 <= len(atoms) <= max_atoms_per_molecule:
            raise ValueError("source molecule exceeds max_atoms_per_molecule")
        largest = max(largest, len(atoms))
        work += len(atoms)*(len(atoms)-1)//2
        if work > max_intramolecular_pairs:
            raise ValueError("aggregate intramolecular pair-work budget exceeded")
    return work, largest


def _assignment_payload(cell_record, *, max_atoms, compact=False,
                        max_atoms_per_molecule=None, max_intramolecular_pairs=None):
    _integer(max_atoms, "max_atoms", 1)
    if not isinstance(cell_record, dict):
        raise ValueError("cell_record must be a packed-cell mapping, not an envelope")
    if not isinstance(cell_record.get("atoms"), list) or not cell_record["atoms"]:
        raise ValueError("the source cell requires a nonempty atom array")
    if len(cell_record["atoms"]) > max_atoms:
        raise ValueError("source cell atom count exceeds max_atoms")
    if compact:
        local_pair_work, largest_molecule = _compact_preflight(cell_record, max_atoms,
            max_atoms_per_molecule, max_intramolecular_pairs)
    try:
        cell = json.loads(_canonical(cell_record))
    except (ValueError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError("source cell must contain finite JSON data") from exc
    if (cell.get("schema_version") != "1" or cell.get("kind") != _CELL_KIND
            or cell.get("coordinate_units") != "angstrom"):
        raise ValueError("unsupported packed-cell schema, kind or coordinate units")
    box = cell.get("box")
    if not isinstance(box, dict) or set(box) != {"origin_angstrom", "lengths_angstrom", "periodic", "volume_angstrom3"}:
        raise ValueError("source box must be an explicit orthorhombic box")
    if _triple(box["origin_angstrom"], "box origin") != (0., 0., 0.):
        raise ValueError("only zero-origin boxes are supported")
    lengths = _triple(box["lengths_angstrom"], "box lengths")
    if any(length <= 0 for length in lengths):
        raise ValueError("box lengths must be positive")
    if not isinstance(box["periodic"], list) or len(box["periodic"]) != 3 or any(v is not True for v in box["periodic"]):
        raise ValueError("all three box axes must be explicitly periodic")
    volume = math.prod(lengths)
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("box volume is not positive and representable")
    declared_volume = _number(box["volume_angstrom3"], "declared box volume", True)
    if not math.isclose(volume, declared_volume, rel_tol=1e-12, abs_tol=1e-12):
        raise ValueError("declared volume differs from the orthorhombic box")
    atoms, molecules, source_bonds = cell["atoms"], cell.get("molecules"), cell.get("bonds")
    if not isinstance(molecules, list) or not molecules or not isinstance(source_bonds, list):
        raise ValueError("source cell requires molecule and bond arrays")
    for name, actual in (("atom_count", len(atoms)), ("molecule_count", len(molecules)), ("bond_count", len(source_bonds))):
        if _integer(cell.get(name), name, 1) != actual:
            raise ValueError(f"source {name} differs from its records")
    if len(molecules) > len(atoms):
        raise ValueError("more molecules than atoms are unsupported")

    atom_ids, atom_maps = set(), set()
    grouped = {index: {} for index in range(len(molecules))}
    for atom in atoms:
        if not isinstance(atom, dict) or set(atom) != _ATOM_FIELDS:
            raise ValueError("source cell atom fields differ from the coordinate contract")
        identity = atom["cell_atom_id"]
        if not isinstance(identity, str) or not identity or identity in atom_ids:
            raise ValueError("cell atom IDs must be nonempty and unique")
        mapping = _integer(atom["atom_map"], "global atom map", 1)
        if mapping in atom_maps:
            raise ValueError("global atom maps must be unique")
        atom_ids.add(identity); atom_maps.add(mapping)
        molecule_index = _integer(atom["molecule_index"], "molecule index")
        if molecule_index not in grouped:
            raise ValueError("atom molecule index is outside the molecule array")
        source_id = atom["source_atom_id"]
        if not isinstance(source_id, str) or not source_id or source_id in grouped[molecule_index]:
            raise ValueError("source atom IDs must be unique within each molecule")
        if identity != f"molecule:{molecule_index}/{source_id}":
            raise ValueError("global cell atom ID does not preserve its molecule/source ID")
        _integer(atom["source_atom_map"], "source atom map", 1)
        _integer(atom["isotope"], "isotope")
        if type(atom["formal_charge"]) is not int:
            raise ValueError("formal charge must be an integer")
        _number(atom["mass_g_mol"], "atom mass", True)
        position = _triple(atom["coordinates_angstrom"], "unwrapped coordinates")
        wrapped = _triple(atom["wrapped_coordinates_angstrom"], "wrapped coordinates")
        images = atom["image_flags"]
        if not isinstance(images, list) or len(images) != 3 or any(type(v) is not int for v in images):
            raise ValueError("image flags must contain three integers")
        for value, w, image, length in zip(position, wrapped, images, lengths):
            if not 0 <= w < length:
                raise ValueError("wrapped coordinates must lie in the half-open box")
            try:
                restored = w + image*length
            except OverflowError as exc:
                raise ValueError("image reconstruction is not representable") from exc
            _close(value, restored, "wrapped/unwrapped/image relationship")
        grouped[molecule_index][source_id] = atom

    assignments, atom_types, expected_bonds = [], {}, {}
    parameter_set = source_snapshot = None
    assembled = {kind: [] for kind in ("bonds", "angles", "dihedrals", "impropers")}
    intramolecular_pairs = {}
    compact_counts = {"intermolecular": 0, "excluded": 0, "scaled_14": 0,
                      "unscaled_intramolecular": 0}
    for index, molecule in enumerate(molecules):
        if not isinstance(molecule, dict) or _integer(molecule.get("molecule_index"), "molecule index") != index:
            raise ValueError("source molecule records must have unique ordered indices")
        if not isinstance(molecule.get("source_geometry"), dict) or len(molecule["source_geometry"].get("atoms", [])) > max_atoms:
            raise ValueError("source molecule geometry exceeds the finite atom budget")
        geometry = MolecularGeometry.from_mapping(molecule["source_geometry"])
        if molecule.get("source_graph_sha256") != geometry.source_graph_sha256:
            raise ValueError("molecule source graph digest changed")
        current = grouped[index]
        source_by_id = {atom.atom_id: atom for atom in geometry.atoms}
        if set(current) != set(source_by_id):
            raise ValueError("source and packed atom identities differ")
        placement = molecule.get("placement")
        if not isinstance(placement, dict):
            raise ValueError("source rigid placement is required")
        raw_matrix = placement.get("rotation_matrix")
        if not isinstance(raw_matrix, list) or len(raw_matrix) != 3:
            raise ValueError("a three-by-three rotation matrix is required")
        matrix = tuple(_triple(row, "rotation matrix") for row in raw_matrix)
        _close(_determinant(matrix), 1., "proper rotation determinant", _MATRIX_TOLERANCE)
        for i, row in enumerate(matrix):
            for j, other in enumerate(matrix):
                _close(math.fsum(a*b for a,b in zip(row,other)), float(i == j),
                       "rotation orthogonality", _MATRIX_TOLERANCE)
        center = _triple(placement.get("source_centroid_angstrom"), "source centroid")
        translation = _triple(placement.get("translated_centroid_angstrom"), "translated centroid")
        for axis, declared in enumerate(center):
            actual = math.fsum(point[axis]/len(geometry.atoms) for point in geometry.coordinates_angstrom)
            _close(declared, actual, "source geometric centroid")
        for original, source_position in zip(geometry.atoms, geometry.coordinates_angstrom):
            record = current[original.atom_id]
            if record["source_atom_map"] != original.atom_map:
                raise ValueError("source atom maps changed")
            for key in ("element", "isotope", "formal_charge"):
                if type(record[key]) is not type(getattr(original,key)) or record[key] != getattr(original,key):
                    raise ValueError(f"source atomic {key} changed")
            _close(float(record["mass_g_mol"]), original.mass_g_mol, "source atomic mass", 1e-12)
            if _canonical(record["lineage"]) != _canonical(original.lineage):
                raise ValueError("source atomic lineage changed")
            try:
                centered = tuple(value-c for value,c in zip(source_position,center))
                transformed = tuple(math.fsum(a*b for a,b in zip(row,centered))+t
                                    for row,t in zip(matrix,translation))
            except (OverflowError, ValueError) as exc:
                raise ValueError("source rigid transform is not representable") from exc
            for value, expected in zip(record["coordinates_angstrom"], transformed):
                _close(float(value), expected, "source rigid coordinates")
        # A formally proper transform can still lose molecular dimensions when
        # huge translations exhaust floating-point coordinate resolution.
        # Check the actual serialized geometry, not only the transform formula.
        for first, second in combinations(range(len(geometry.atoms)),2):
            before = math.dist(geometry.coordinates_angstrom[first], geometry.coordinates_angstrom[second])
            after = math.dist(current[geometry.atoms[first].atom_id]["coordinates_angstrom"],
                              current[geometry.atoms[second].atom_id]["coordinates_angstrom"])
            if not math.isfinite(before) or before <= 0:
                raise ValueError("source molecular positions are coincident or unrepresentable")
            _close(after, before, "actual intramolecular pair distance")
        assigned = assign_alkane(geometry,
            max_atoms=max_atoms_per_molecule if compact else max_atoms).to_dict()
        if parameter_set is None:
            parameter_set, source_snapshot = assigned["parameter_set"], assigned["source_snapshot"]
        elif (_canonical(parameter_set) != _canonical(assigned["parameter_set"])
              or _canonical(source_snapshot) != _canonical(assigned["source_snapshot"])):
            raise ValueError("molecule parameter sources differ")
        if compact:
            assignments.append({"molecule_index": index, "atom_count": len(geometry.atoms),
                "source_assignment_sha256": hashlib.sha256(_canonical(assigned).encode()).hexdigest()})
        else:
            assignments.append({"molecule_index": index, "source_assignment": assigned})
        for row in assigned["atom_assignments"]:
            atom_types[current[row["atom_id"]]["cell_atom_id"]] = (row["type_id"], row["charge_e"])
        for kind, rows in assigned["interactions"].items():
            for row in rows:
                global_ids = [current[name]["cell_atom_id"] for name in row["atom_ids"]]
                assembled[kind].append({"id": f"molecule:{index}/{row['id']}",
                                        "type_id": row["type_id"], "atom_ids": global_ids})
        for bond in geometry.bonds:
            ids = tuple(current[name]["cell_atom_id"] for name in bond.atom_ids)
            expected_bonds[frozenset(ids)] = {"molecule_index": index, "ids": ids,
                "maps": {name: next(atom["atom_map"] for atom in current.values() if atom["cell_atom_id"] == name)
                         for name in ids}, "order": bond.order, "aromatic": bond.aromatic, "stereo": bond.stereo}
        for pair in assigned["nonbonded_pairs"]:
            if compact:
                compact_counts["excluded" if pair["lj_scale"] == 0 else
                    ("scaled_14" if pair["lj_scale"] == .5 else "unscaled_intramolecular")] += 1
                if pair["topological_distance"] > 3:
                    continue
            ids = [current[name]["cell_atom_id"] for name in pair["atom_ids"]]
            intramolecular_pairs[frozenset(ids)] = pair

    seen_bonds = set()
    for bond in source_bonds:
        if not isinstance(bond,dict) or set(bond) != {"molecule_index", "cell_atom_ids", "atom_maps", "order", "aromatic", "stereo"}:
            raise ValueError("invalid source cell bond fields")
        ids, maps = bond["cell_atom_ids"], bond["atom_maps"]
        if (not isinstance(ids,list) or len(ids) != 2 or not all(isinstance(v,str) for v in ids)
                or len(set(ids)) != 2 or not isinstance(maps,list) or len(maps) != 2
                or any(type(v) is not int or v <= 0 for v in maps)):
            raise ValueError("source cell bond endpoints/maps are invalid")
        key = frozenset(ids)
        if key in seen_bonds or key not in expected_bonds:
            raise ValueError("duplicate, missing or intermolecular source bond")
        seen_bonds.add(key); expected = expected_bonds[key]
        if type(bond["molecule_index"]) is not int or bond["molecule_index"] != expected["molecule_index"]:
            raise ValueError("bond molecule identity changed")
        if dict(zip(ids,maps)) != expected["maps"]:
            raise ValueError("bond global maps changed")
        _number(bond["order"], "bond order", True)
        if (bond["order"] != expected["order"] or type(bond["aromatic"]) is not bool
                or bond["aromatic"] != expected["aromatic"] or bond["stereo"] != expected["stereo"]):
            raise ValueError("source chemical bond attributes changed")
    if seen_bonds != set(expected_bonds):
        raise ValueError("source cell lost a molecular bond")

    output_atoms = [{**atom, "type_id": atom_types[atom["cell_atom_id"]][0],
                    "charge_e": atom_types[atom["cell_atom_id"]][1]} for atom in atoms]
    pairs, counts = [], {"intermolecular": 0, "excluded": 0, "scaled_14": 0, "unscaled_intramolecular": 0}
    if compact:
        counts = compact_counts
        counts["intermolecular"] = len(output_atoms)*(len(output_atoms)-1)//2-local_pair_work
        order = {atom["cell_atom_id"]: index for index, atom in enumerate(output_atoms)}
        for key, original in intramolecular_pairs.items():
            ids = sorted(key, key=order.__getitem__)
            pairs.append({"atom_ids": ids, "same_molecule": True,
                          "topological_distance": original["topological_distance"],
                          "lj_scale": original["lj_scale"], "coulomb_scale": original["coulomb_scale"]})
        pairs.sort(key=lambda pair: tuple(order[identity] for identity in pair["atom_ids"]))
    else:
        for first, second in combinations(output_atoms,2):
            ids = [first["cell_atom_id"], second["cell_atom_id"]]
            same = first["molecule_index"] == second["molecule_index"]
            if same:
                original = intramolecular_pairs[frozenset(ids)]
                distance, lj_scale, q_scale = original["topological_distance"], original["lj_scale"], original["coulomb_scale"]
                counts["excluded" if lj_scale == 0 else ("scaled_14" if lj_scale == .5 else "unscaled_intramolecular")] += 1
            else:
                distance, lj_scale, q_scale = None, 1., 1.
                counts["intermolecular"] += 1
            pairs.append({"atom_ids": ids, "same_molecule": same, "topological_distance": distance,
                          "lj_scale": lj_scale, "coulomb_scale": q_scale})
    charge = math.fsum(atom["charge_e"] for atom in output_atoms)
    mass = math.fsum(atom["mass_g_mol"] for atom in output_atoms)
    density = density_from_mass_volume([atom["mass_g_mol"] for atom in output_atoms], volume)
    if abs(charge) > 1e-12:
        raise ValueError("assembled partial charges are not neutral")
    if type(cell.get("formal_charge")) is not int or cell["formal_charge"] != sum(a["formal_charge"] for a in output_atoms):
        raise ValueError("source total formal charge changed")
    for name, expected in (("total_molar_mass_g_mol", mass), ("construction_density_g_cm3", density)):
        actual = _number(cell.get(name), name, True)
        if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"source {name} differs from the actual cell")
    payload = {
        "schema_version": "1", "kind": _ASSIGNED_KIND, "model_id": MODEL_ID,
        "source_cell": cell, "box": box, "parameter_set": parameter_set,
        "atoms": output_atoms, "interactions": assembled, "nonbonded_pairs": pairs,
        "molecules": assignments, "source_snapshot": source_snapshot,
        "validation": {"source_coordinate_and_image_identity": True,
            "proper_source_rigid_transforms": True, "source_graph_bonds_and_lineage_preserved": True,
            "actual_intramolecular_pair_distances_preserved": True,
            "chemical_typing": "performed_per_supported_neutral_explicit_n_alkane_molecule",
            "molecule_count": len(molecules), "atom_count": len(atoms),
            "topology_counts": {kind: len(rows) for kind,rows in assembled.items()},
            "all_unordered_pairs_once": True, "pair_scaling_counts": counts,
            "total_partial_charge_e": charge, "formal_charge_e": 0,
            "total_molar_mass_g_mol": mass, "construction_density_g_cm3": density,
            "density_interpretation": "mass divided by prescribed box volume, not predicted density",
            "coordinate_absolute_tolerance_angstrom": _COORDINATE_TOLERANCE,
            "rotation_absolute_tolerance": _MATRIX_TOLERANCE,
            "source_validation_flags": "untrusted supplied metadata; no reliance on their values",
            "source_nonbonded_contacts": "not re-certified by parameter assignment",
            "force_field_assigned": True, "production_ready": False,
            "cutoff_and_long_range_model": "not_selected_by_assignment",
            "simulation_performed": False, "equilibrium_status": "not_assessed",
            "material_property_validation": "not_performed"},
    }
    if compact:
        payload.update(schema_version="2", kind=_COMPACT_KIND,
                       special_pairs=payload.pop("nonbonded_pairs"),
                       nonbonded_rules=json.loads(_canonical(_COMPACT_RULES)))
        payload["validation"].pop("all_unordered_pairs_once")
        payload["validation"].update(
            implicit_pair_rule_coverage="all_unordered_pairs_via_unit_default_and_unique_exceptions",
            explicit_global_all_pairs_materialized=False,
            special_pair_count=len(pairs),
            total_unordered_pair_count=len(output_atoms)*(len(output_atoms)-1)//2,
            intramolecular_pair_work=local_pair_work,
            largest_molecule_atom_count=largest_molecule,
            implementation_limits={"max_atoms": _COMPACT_MAX_ATOMS,
                "max_atoms_per_molecule": _COMPACT_MAX_MOLECULE_ATOMS,
                "max_intramolecular_pairs": _COMPACT_MAX_INTRAMOLECULAR_PAIRS},
            molecular_source_assignments="recomputed_and_hash_bound_not_duplicated",
        )
    return payload


def assign_periodic_alkanes(cell_record: dict, *, max_atoms: int = 512) -> AssignedPeriodicAlkaneCell:
    """Assemble sourced assignments without altering atoms, images or positions.

    The finite all-pair topology budget is explicit. This prepares model identity
    and topology only; a periodic electrostatics model and run protocol remain
    separate decisions. Molecular rigidity/source checks are for starting cells,
    not a claim that a later thermal trajectory must remain rigid.
    """
    try:
        payload = _assignment_payload(cell_record, max_atoms=max_atoms)
    except (KeyError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError(f"invalid periodic source cell: {exc}") from exc
    return AssignedPeriodicAlkaneCell._from_validated(payload)


def assign_compact_periodic_alkanes(cell_record: dict, *, max_atoms: int = _COMPACT_MAX_ATOMS,
        max_atoms_per_molecule: int = _COMPACT_MAX_MOLECULE_ATOMS,
        max_intramolecular_pairs: int = _COMPACT_MAX_INTRAMOLECULAR_PAIRS) -> CompactAssignedPeriodicAlkaneCell:
    """Assign the same frozen model using sparse short-topology pair exceptions.

    All omitted unordered pairs have LJ and Coulomb scale one, including every
    intermolecular pair. The limits bound global atoms, local quadratic work and
    aggregate local assignment work before source copying or graph assignment.
    The record is still a starting-cell assignment, not an equilibrium state.
    """
    try:
        payload = _assignment_payload(cell_record, max_atoms=max_atoms, compact=True,
            max_atoms_per_molecule=max_atoms_per_molecule,
            max_intramolecular_pairs=max_intramolecular_pairs)
    except (KeyError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError(f"invalid compact periodic source cell: {exc}") from exc
    return CompactAssignedPeriodicAlkaneCell._from_validated(payload)
