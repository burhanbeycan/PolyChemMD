"""Sourced OPLS-AA assignments for a deliberately restricted n-alkane graph.

The project code independently identifies CH3/CH2/HC environments, enumerates
the complete graph topology and assigns the pinned numerical parameter subset.
It performs no fitting, energy evaluation, relaxation or molecular dynamics.
"""

from collections import Counter, deque
from dataclasses import InitVar, dataclass
from itertools import combinations
import json
import math

from ..geometry import MolecularGeometry
from .alkane_data import load_alkane_parameters
from .contracts import ParameterSet, TopologyRequirements, check_coverage


MODEL_ID = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"
MASS_TOLERANCE_G_MOL = 1e-6


class UnsupportedAlkaneError(ValueError):
    """The supplied chemistry is outside the explicitly supported subset."""


@dataclass(frozen=True)
class AssignedAlkane:
    """Immutable serialized assignment; the returned mapping is a fresh copy.

    Instances returned by :func:`assign_alkane` have passed its graph/type and
    coverage checks. The payload is not a material-validation certificate.
    """

    _payload_json: str
    max_atoms: InitVar[int] = 256

    def __post_init__(self, max_atoms):
        """Revalidate serialized records against the exact frozen assignment.

        The public constructor is a round-trip boundary, not an unchecked way
        to attach this model's name to changed charges or parameter values.
        Geometry provenance remains supplied data; actual assignment is rebuilt
        from its graph. The same explicit atom budget bounds regeneration.
        """
        if not isinstance(self._payload_json, str):
            raise ValueError("assigned n-alkane record must be JSON object text")
        def reject_constant(value):
            raise ValueError("assigned n-alkane JSON must contain finite numbers")
        def unique(pairs):
            data = {}
            for key, value in pairs:
                if key in data:
                    raise ValueError("assigned n-alkane JSON contains duplicate fields")
                data[key] = value
            return data
        try:
            data = json.loads(self._payload_json, parse_constant=reject_constant,
                              object_pairs_hook=unique)
            if not isinstance(data, dict) or "geometry" not in data:
                raise ValueError("assigned n-alkane record requires its geometry")
            geometry = MolecularGeometry.from_mapping(data["geometry"])
            expected = _assignment_payload(geometry, max_atoms=max_atoms)
            # Coordinate numbers are normalized by their geometry contract;
            # changing positions alone does not change an atom-type model.
            data["geometry"] = geometry.to_dict()
            actual_json = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)
            expected_json = json.dumps(expected, sort_keys=True, separators=(",", ":"), allow_nan=False)
        except (ValueError, TypeError, KeyError, OverflowError, RecursionError) as exc:
            raise ValueError(f"invalid assigned n-alkane record: {exc}") from exc
        if actual_json != expected_json:
            raise ValueError("assigned record does not match the frozen n-alkane model, topology and charges")
        object.__setattr__(self, "_payload_json", expected_json)

    @classmethod
    def _from_validated(cls, payload):
        """Internal factory after `_assignment_payload` performed all checks."""
        instance = object.__new__(cls)
        object.__setattr__(instance, "_payload_json",
                           json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False))
        return instance

    def to_dict(self):
        return json.loads(self._payload_json)

    @property
    def geometry(self):
        return MolecularGeometry.from_mapping(self.to_dict()["geometry"])


def _adjacency(geometry):
    adjacent = [set() for _ in geometry.atoms]
    index = {atom.atom_id: i for i, atom in enumerate(geometry.atoms)}
    for bond in geometry.bonds:
        if bond.order != 1.0 or bond.aromatic or bond.stereo != "STEREONONE":
            raise UnsupportedAlkaneError("only nonaromatic single bonds with no bond-stereo assignment are supported")
        first, second = (index[name] for name in bond.atom_ids)
        adjacent[first].add(second)
        adjacent[second].add(first)
    return adjacent


def _validate_chemistry(geometry, adjacent):
    atoms = geometry.atoms
    if any(atom.element not in {"C", "H"} for atom in atoms):
        raise UnsupportedAlkaneError("the supported n-alkane subset contains only explicit C and H atoms")
    if any(atom.formal_charge != 0 for atom in atoms):
        raise UnsupportedAlkaneError("charged atoms are outside the neutral n-alkane subset")
    if any(atom.isotope != 0 for atom in atoms):
        raise UnsupportedAlkaneError("isotope-specific parameter masses are unsupported in this subset")
    carbons = [i for i, atom in enumerate(atoms) if atom.element == "C"]
    if len(carbons) < 2:
        raise UnsupportedAlkaneError("at least two carbon atoms are required; methane is a separate OPLS type")
    carbon_degrees = {}
    for i, atom in enumerate(atoms):
        if atom.element == "H":
            if len(adjacent[i]) != 1 or atoms[next(iter(adjacent[i]))].element != "C":
                raise UnsupportedAlkaneError("every explicit H must have exactly one carbon neighbour")
            continue
        carbon_degree = sum(atoms[j].element == "C" for j in adjacent[i])
        if len(adjacent[i]) != 4 or carbon_degree not in (1, 2):
            raise UnsupportedAlkaneError("each carbon must be saturated CH3/CH2 in an unbranched carbon path")
        carbon_degrees[i] = carbon_degree
    if Counter(carbon_degrees.values()).get(1, 0) != 2:
        raise UnsupportedAlkaneError("an acyclic n-alkane carbon path must have exactly two CH3 endpoints")
    carbon_edges = sum(carbon_degrees.values()) // 2
    if carbon_edges != len(carbons) - 1:
        raise UnsupportedAlkaneError("cycles are outside the supported acyclic n-alkane subset")
    reached, frontier = set(), [carbons[0]]
    while frontier:
        current = frontier.pop()
        if current in reached:
            continue
        reached.add(current)
        frontier.extend(j for j in adjacent[current] if atoms[j].element == "C" and j not in reached)
    if len(reached) != len(carbons):
        raise UnsupportedAlkaneError("the carbon skeleton must be one connected path")
    if len(atoms) != 3 * len(carbons) + 2:
        raise UnsupportedAlkaneError("explicit atom balance does not match C_n H_(2n+2)")
    return carbons, carbon_degrees


def _interaction(kind, rows):
    return [{"id": f"{kind}:{i + 1}", "type_id": type_id, "atom_ids": atom_ids}
            for i, (type_id, atom_ids) in enumerate(rows)]


def _shortest_distances(start, adjacent):
    distances, queue = {start: 0}, deque([start])
    while queue:
        current = queue.popleft()
        for neighbour in sorted(adjacent[current]):
            if neighbour not in distances:
                distances[neighbour] = distances[current] + 1
                queue.append(neighbour)
    return distances


def _assignment_payload(geometry: MolecularGeometry, *, max_atoms: int):
    """Build checked assignment data, shared by construction and round trips."""
    if not isinstance(geometry, MolecularGeometry):
        raise ValueError("geometry must be a MolecularGeometry record")
    if type(max_atoms) is not int or max_atoms <= 0:
        raise ValueError("max_atoms must be a positive integer, not a boolean")
    if len(geometry.atoms) > max_atoms:
        raise ValueError("atom count exceeds max_atoms for explicit topology/pair enumeration")
    adjacent = _adjacency(geometry)
    carbons, carbon_degrees = _validate_chemistry(geometry, adjacent)
    data = load_alkane_parameters()
    if data.get("model_id") != MODEL_ID:
        raise ValueError("loaded parameter data do not match the supported frozen model")
    parameters = ParameterSet.from_mapping(data["parameter_set"])
    if parameters.source.kind != "literature" or any(row.source.kind != "literature" for row in parameters.definitions):
        raise ValueError("n-alkane assignment requires the attributed literature parameter subset")
    atom_definitions = {row.type_id: row for row in parameters.definitions if row.kind == "atoms"}
    charges = data["charge_assignments"]
    selected_types = {"opls_135", "opls_136", "opls_140"}
    if not isinstance(charges, dict) or not selected_types.issubset(charges):
        raise ValueError("source charge assignments are incomplete")
    for type_id in selected_types:
        value = charges[type_id]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError("source partial charges must be finite numbers")
        try:
            numeric = float(value)
        except (ValueError, OverflowError) as exc:
            raise ValueError("source partial charges must be finite numbers") from exc
        if not math.isfinite(numeric):
            raise ValueError("source partial charges must be finite numbers")
        if type_id not in atom_definitions:
            raise ValueError("source atom parameter definitions are incomplete")

    atom_assignments = []
    for i, atom in enumerate(geometry.atoms):
        type_id = "opls_140" if atom.element == "H" else ("opls_135" if carbon_degrees[i] == 1 else "opls_136")
        parameter_mass = atom_definitions[type_id].values["mass_g_mol"]
        if not math.isclose(atom.mass_g_mol, parameter_mass, rel_tol=0.0, abs_tol=MASS_TOLERANCE_G_MOL):
            raise UnsupportedAlkaneError("geometry atomic masses must match the attributed parameter masses; no silent mass replacement")
        atom_assignments.append({"atom_id": atom.atom_id, "atom_map": atom.atom_map,
                                 "type_id": type_id, "charge_e": charges[type_id]})
    charge_sum = math.fsum(row["charge_e"] for row in atom_assignments)
    if not math.isfinite(charge_sum) or abs(charge_sum) > 1e-12:
        raise ValueError("source partial charges failed neutral n-alkane charge balance")

    ids = [atom.atom_id for atom in geometry.atoms]
    elements = [atom.element for atom in geometry.atoms]
    bond_rows, angle_rows, dihedral_rows = [], [], []
    for first in range(len(ids)):
        for second in sorted(j for j in adjacent[first] if j > first):
            type_id = "bond_CC" if elements[first] == elements[second] == "C" else "bond_CH"
            bond_rows.append((type_id, [ids[first], ids[second]]))
    for center in range(len(ids)):
        for first, last in combinations(sorted(adjacent[center]), 2):
            outer_carbons = int(elements[first] == "C") + int(elements[last] == "C")
            type_id = ("angle_HCH", "angle_CCH", "angle_CCC")[outer_carbons]
            angle_rows.append((type_id, [ids[first], ids[center], ids[last]]))
    for center in carbons:
        for next_center in sorted(j for j in adjacent[center] if j > center and elements[j] == "C"):
            for first in sorted(adjacent[center] - {next_center}):
                for last in sorted(adjacent[next_center] - {center}):
                    outer_carbons = int(elements[first] == "C") + int(elements[last] == "C")
                    type_id = ("torsion_HCCH", "torsion_CCCH", "torsion_CCCC")[outer_carbons]
                    dihedral_rows.append((type_id, [ids[first], ids[center], ids[next_center], ids[last]]))
    interactions = {"bonds": _interaction("bond", bond_rows),
                    "angles": _interaction("angle", angle_rows),
                    "dihedrals": _interaction("dihedral", dihedral_rows), "impropers": []}
    requirements = {"atoms": [{"id": row["atom_id"], "type_id": row["type_id"]} for row in atom_assignments],
                    **{kind: [{"id": row["id"], "type_id": row["type_id"]} for row in rows]
                       for kind, rows in interactions.items()}}
    coverage = check_coverage(parameters, TopologyRequirements.from_mapping(requirements)).require_complete()
    pairs, pair_counts = [], Counter()
    for first in range(len(ids)):
        distances = _shortest_distances(first, adjacent)
        for second in range(first + 1, len(ids)):
            distance = distances[second]
            scale = 0.0 if distance <= 2 else (0.5 if distance == 3 else 1.0)
            pairs.append({"atom_ids": [ids[first], ids[second]], "topological_distance": distance,
                          "lj_scale": scale, "coulomb_scale": scale})
            pair_counts["excluded" if scale == 0 else ("scaled_14" if scale == .5 else "unscaled")] += 1
    result = {
        "schema_version": "1", "kind": "assigned_opls_aa_n_alkane", "model_id": MODEL_ID,
        "geometry": geometry.to_dict(), "parameter_set": parameters.to_dict(),
        "atom_assignments": atom_assignments, "interactions": interactions,
        "nonbonded_pairs": pairs, "coverage": coverage.to_dict(),
        "source_snapshot": data["source_snapshot"],
        "validation": {
            "chemical_typing": "performed_for_supported_neutral_n_alkane_graph",
            "explicit_hydrogen_valence": True, "carbon_skeleton": "single_acyclic_path",
            "carbon_count": len(carbons), "molecular_formula": f"C{len(carbons)}H{2*len(carbons)+2}",
            "topology_enumeration": "complete_bonds_angles_and_all_proper_four_atom_paths",
            "improper_policy": "no_improper_for_saturated_alkane_subset",
            "partial_charge_model": "attributed_fixed_OPLS_atom_type_charges_not_QM_fitted_here",
            "total_partial_charge_e": charge_sum, "formal_charge_e": 0,
            "pair_scaling_counts": {key: pair_counts[key] for key in ("excluded", "scaled_14", "unscaled")},
            "mass_policy": "geometry_mass_matches_source_without_override",
            "mass_absolute_tolerance_g_mol": MASS_TOLERANCE_G_MOL,
            "force_field_assigned": True, "parameters_assigned": True,
            "physical_accuracy": "not_assessed",
            "production_ready": False, "equilibrium_status": "not_assessed",
            "simulation_performed": False,
            "coverage_report_scope": "coverage only; chemical typing is recorded by this outer validation object",
        },
    }
    return result


def assign_alkane(geometry: MolecularGeometry, *, max_atoms: int = 256) -> AssignedAlkane:
    """Type a neutral explicit n-alkane and enumerate all assigned interactions.

    Allowed atoms are unlabelled-isotope saturated C/H with source-compatible
    masses. The carbon skeleton must be a single path containing at least two
    carbons. Pair records include *all* unordered atom pairs: graph distances
    one/two are excluded, distance three is scaled by 0.5, and larger distances
    have unit scaling. All proper four-atom paths are retained, including
    hydrogen-ended torsions. Coordinates and input lineage remain unchanged.
    """
    return AssignedAlkane._from_validated(_assignment_payload(geometry, max_atoms=max_atoms))
