"""Immutable coordinate records, independent of cheminformatics and engines.

These checks establish a well-formed coordinate/topology record. They do not
infer a chemical graph, verify parameter masses, or certify physical geometry.
The conformer and packing layers supply their own method-specific checks.
"""

from copy import deepcopy
from dataclasses import dataclass
import json
import math
import re


_ELEMENTS = frozenset((
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni "
    "Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe "
    "Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au "
    "Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf "
    "Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og"
).split())


def _integer(value, name, minimum=None):
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer, not a boolean")
    if minimum is not None and value < minimum:
        raise ValueError(f"{name} must be >= {minimum}")
    return value


def _number(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be a representable finite number") from exc
    if not math.isfinite(value) or (positive and value <= 0):
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return value


def _identifier(value, name):
    if not isinstance(value, str) or not value.strip() or any(ord(c) < 32 for c in value):
        raise ValueError(f"{name} must be a nonempty identifier without control characters")
    return value


def _json_object(raw, name):
    if not isinstance(raw, str):
        raise ValueError(f"{name} must be JSON object text")
    def constant(value):
        raise ValueError(f"{name} contains a nonfinite JSON value")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"{name} contains duplicate JSON fields")
            result[key] = value
        return result
    try:
        data = json.loads(raw, parse_constant=constant, object_pairs_hook=unique)
        # This also rejects overflow such as JSON's 1e400 in nested metadata.
        canonical = json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False)
    except (ValueError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError(f"{name} must be finite, unambiguous JSON") from exc
    if not isinstance(data, dict):
        raise ValueError(f"{name} must encode a JSON object")
    return data, canonical


@dataclass(frozen=True)
class GeometryAtom:
    atom_id: str
    atom_map: int
    element: str
    isotope: int
    formal_charge: int
    mass_g_mol: float
    lineage_json: str

    def __post_init__(self):
        _identifier(self.atom_id, "atom_id")
        _integer(self.atom_map, "atom_map", 1)
        if not isinstance(self.element, str) or self.element not in _ELEMENTS:
            raise ValueError("element must be a chemical element symbol")
        _integer(self.isotope, "isotope", 0)
        _integer(self.formal_charge, "formal_charge")
        object.__setattr__(self, "mass_g_mol", _number(self.mass_g_mol, "mass_g_mol", True))
        lineage, canonical = _json_object(self.lineage_json, "lineage_json")
        for name in ("atom_id", "atom_map", "element", "isotope", "formal_charge"):
            if name in lineage and (type(lineage[name]) is not type(getattr(self, name))
                                    or lineage[name] != getattr(self, name)):
                raise ValueError(f"lineage_json {name} conflicts with its atom")
        object.__setattr__(self, "lineage_json", canonical)

    @property
    def lineage(self):
        return json.loads(self.lineage_json)

    def to_dict(self):
        return {"atom_id": self.atom_id, "atom_map": self.atom_map,
                "element": self.element, "isotope": self.isotope,
                "formal_charge": self.formal_charge, "mass_g_mol": self.mass_g_mol,
                "lineage": self.lineage}


@dataclass(frozen=True)
class GeometryBond:
    atom_ids: tuple[str, str]
    order: float
    aromatic: bool = False
    stereo: str = "STEREONONE"

    def __post_init__(self):
        if not isinstance(self.atom_ids, (tuple, list)) or len(self.atom_ids) != 2:
            raise ValueError("bond atom_ids must contain two endpoints")
        endpoints = tuple(_identifier(a, "bond atom_id") for a in self.atom_ids)
        if endpoints[0] == endpoints[1]:
            raise ValueError("self bonds are unsupported")
        order = _number(self.order, "bond order", True)
        if order not in (1.0, 1.5, 2.0, 3.0):
            raise ValueError("unsupported bond order")
        if not isinstance(self.aromatic, bool) or (self.aromatic and order != 1.5):
            raise ValueError("aromatic must be boolean and aromatic bonds must have order 1.5")
        if not isinstance(self.stereo, str) or self.stereo not in {"STEREONONE", "STEREOANY", "STEREOZ", "STEREOE", "STEREOCIS", "STEREOTRANS"}:
            raise ValueError("unsupported bond stereochemistry")
        object.__setattr__(self, "atom_ids", endpoints)
        object.__setattr__(self, "order", order)

    def to_dict(self):
        return {"atom_ids": list(self.atom_ids), "order": self.order,
                "aromatic": self.aromatic, "stereo": self.stereo}


@dataclass(frozen=True)
class MolecularGeometry:
    atoms: tuple[GeometryAtom, ...]
    bonds: tuple[GeometryBond, ...]
    coordinates_angstrom: tuple[tuple[float, float, float], ...]
    source_graph_sha256: str
    provenance_json: str

    def __post_init__(self):
        if not isinstance(self.atoms, (tuple, list)) or not self.atoms:
            raise ValueError("atoms must be a nonempty sequence")
        if not all(isinstance(a, GeometryAtom) for a in self.atoms):
            raise ValueError("atoms must contain GeometryAtom records")
        atoms = tuple(self.atoms)
        ids, maps = [a.atom_id for a in atoms], [a.atom_map for a in atoms]
        if len(set(ids)) != len(ids) or len(set(maps)) != len(maps):
            raise ValueError("atom IDs and maps must each be unique")
        if not isinstance(self.bonds, (tuple, list)) or not all(isinstance(b, GeometryBond) for b in self.bonds):
            raise ValueError("bonds must contain GeometryBond records")
        bonds, seen = tuple(self.bonds), set()
        adjacent = {atom_id: set() for atom_id in ids}
        for bond in bonds:
            pair = frozenset(bond.atom_ids)
            if pair in seen or not pair.issubset(adjacent):
                raise ValueError("bond endpoints must exist and bond pairs must be unique")
            seen.add(pair)
            first, second = bond.atom_ids
            adjacent[first].add(second)
            adjacent[second].add(first)
        visited, frontier = set(), [ids[0]]
        while frontier:
            current = frontier.pop()
            if current not in visited:
                visited.add(current)
                frontier.extend(adjacent[current] - visited)
        if len(visited) != len(atoms):
            raise ValueError("one MolecularGeometry must describe one connected molecule")
        if not isinstance(self.coordinates_angstrom, (tuple, list)) or len(self.coordinates_angstrom) != len(atoms):
            raise ValueError("coordinates must match the atom sequence")
        coordinates = []
        for row in self.coordinates_angstrom:
            if not isinstance(row, (tuple, list)) or len(row) != 3:
                raise ValueError("each coordinate must have exactly three components")
            coordinates.append(tuple(_number(v, "coordinate") for v in row))
        if not isinstance(self.source_graph_sha256, str) or re.fullmatch(r"[0-9a-f]{64}", self.source_graph_sha256) is None:
            raise ValueError("source_graph_sha256 must be a lowercase SHA-256 digest")
        _, provenance = _json_object(self.provenance_json, "provenance_json")
        object.__setattr__(self, "atoms", atoms)
        object.__setattr__(self, "bonds", bonds)
        object.__setattr__(self, "coordinates_angstrom", tuple(coordinates))
        object.__setattr__(self, "provenance_json", provenance)

    @property
    def provenance(self):
        return json.loads(self.provenance_json)

    def to_dict(self):
        return {"schema_version": "1", "kind": "molecular_geometry",
                "coordinate_units": "angstrom", "atoms": [a.to_dict() for a in self.atoms],
                "bonds": [b.to_dict() for b in self.bonds],
                "coordinates_angstrom": [list(row) for row in self.coordinates_angstrom],
                "source_graph_sha256": self.source_graph_sha256,
                "provenance": self.provenance}

    @classmethod
    def from_mapping(cls, data):
        fields = {"schema_version", "kind", "coordinate_units", "atoms", "bonds",
                  "coordinates_angstrom", "source_graph_sha256", "provenance"}
        if not isinstance(data, dict) or set(data) != fields:
            raise ValueError(f"geometry record requires exactly {sorted(fields)}")
        if data["schema_version"] != "1" or data["kind"] != "molecular_geometry" or data["coordinate_units"] != "angstrom":
            raise ValueError("unsupported geometry schema, kind or coordinate units")
        if not isinstance(data["atoms"], list) or not isinstance(data["bonds"], list):
            raise ValueError("atoms and bonds must be JSON arrays")
        atom_fields = {"atom_id", "atom_map", "element", "isotope", "formal_charge", "mass_g_mol", "lineage"}
        bond_fields = {"atom_ids", "order", "aromatic", "stereo"}
        atoms, bonds = [], []
        for atom in data["atoms"]:
            if not isinstance(atom, dict) or set(atom) != atom_fields:
                raise ValueError("invalid geometry atom fields")
            atoms.append(GeometryAtom(**{k: v for k, v in atom.items() if k != "lineage"},
                                      lineage_json=json.dumps(atom["lineage"], allow_nan=False)))
        for bond in data["bonds"]:
            if not isinstance(bond, dict) or set(bond) != bond_fields:
                raise ValueError("invalid geometry bond fields")
            bonds.append(GeometryBond(**bond))
        try:
            return cls(tuple(atoms), tuple(bonds), deepcopy(data["coordinates_angstrom"]),
                       data["source_graph_sha256"], json.dumps(data["provenance"], allow_nan=False))
        except (TypeError, OverflowError) as exc:
            raise ValueError("invalid geometry record") from exc
