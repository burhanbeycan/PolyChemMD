"""Admission checks for the narrow full-style 12–6 slab pilot format.

LAMMPS data format: https://docs.lammps.org/read_data.html
All registry inputs are operator-supplied numeric data, never LAMMPS programs.
"""

from collections import Counter, defaultdict
from itertools import combinations
import math
import re

ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "F": 18.998,
    "Na": 22.990,
    "Cl": 35.45,
    "S": 32.06,
    "Ar": 39.948,
    "Ti": 47.867,
}

SECTIONS = {
    "Masses": 2,
    "Pair Coeffs": 3,
    "PairIJ Coeffs": 4,
    "Bond Coeffs": 3,
    "Angle Coeffs": 3,
    "Dihedral Coeffs": 5,
    "Improper Coeffs": 3,
    "Atoms": 7,
    "Bonds": 4,
    "Angles": 5,
    "Dihedrals": 6,
    "Impropers": 6,
}
TERMS = {
    "Bonds": ("bond", 2),
    "Angles": ("angle", 3),
    "Dihedrals": ("dihedral", 4),
    "Impropers": ("improper", 4),
}


def _uint(s, positive=True):
    if not re.fullmatch(r"\d+", s):
        raise ValueError(f"Invalid integer: {s}")
    n = int(s)
    if positive and not n:
        raise ValueError("Expected positive integer ID")
    return n


def _float(s):
    x = float(s)
    if not math.isfinite(x):
        raise ValueError("Nonfinite numeric topology parameter")
    return x


def _entries(rows, n, name):
    result = {}
    for tokens in rows:
        i = _uint(tokens[0])
        if i in result:
            raise ValueError(f"Duplicate {name} ID")
        result[i] = tokens[1:]
    if set(result) != set(range(1, n + 1)):
        raise ValueError(f"{name} IDs differ from header counts")
    return result


def parse_data(text):
    rows = defaultdict(list)
    counts = {}
    bounds = {}
    section = None
    lines = text.splitlines()
    if not lines:
        raise ValueError("Empty LAMMPS data file")
    for raw in lines[1:]:
        line = raw.partition("#")[0].strip()
        if not line:
            continue
        if line in SECTIONS:
            if line in rows:
                raise ValueError(f"Duplicate data section: {line}")
            section = line
            rows[section] = []
            continue
        t = line.split()
        if section is None:
            if len(t) in (2, 3) and " ".join(t[1:]) in {
                "atoms",
                "bonds",
                "angles",
                "dihedrals",
                "impropers",
                "atom types",
                "bond types",
                "angle types",
                "dihedral types",
                "improper types",
            }:
                key = " ".join(t[1:])
                if key in counts:
                    raise ValueError(f"Duplicate count: {key}")
                counts[key] = _uint(t[0], False)
                continue
            if len(t) == 4 and (t[2], t[3]) in [(k + "lo", k + "hi") for k in "xyz"]:
                dim = t[2][0]
                if dim in bounds:
                    raise ValueError("Duplicate box extent")
                bounds[dim] = (_float(t[0]), _float(t[1]))
                continue
            raise ValueError(
                f"Unsupported header (orthorhombic cell required): {line[:100]}"
            )
        if section == "Atoms" and len(t) == 10 and t[-3:] == ["0", "0", "0"]:
            t = t[:7]
        if len(t) != SECTIONS[section]:
            raise ValueError(
                f"Unexpected row width in {section}. Class-II, hybrid styles and velocities require another adapter."
            )
        for word in t:
            _float(word)
        rows[section].append(t)
    return rows, counts, bounds


def validate(packed, data, profile):
    supported = (
        "label",
        "parameter_source",
        "applicability",
        "water_model",
        "type_elements",
        "bond_orders",
        "mixing",
        "lj14",
        "coul14",
        "cutoff_A",
        "dt_fs",
        "packed_sha256",
        "schema",
        "water_constraints",
    )
    if set(profile) != set(supported):
        raise ValueError(f"Manifest keys must match: {supported}")
    if profile["schema"] != "studio-full-slab-v1":
        raise ValueError("Unsupported topology manifest schema")
    water = any(a.get("role") == "water" for a in packed.get("atoms", []))
    if profile["water_constraints"] != ("flexible" if water else "no_water"):
        raise ValueError(
            "Specify flexible water with explicit bonded parameters, or no_water for a water-free cell"
        )
    for k in ("label", "parameter_source", "applicability", "water_model"):
        if not isinstance(profile[k], str) or not 1 <= len(profile[k].strip()) <= 500:
            raise ValueError(f"Missing or invalid profile provenance: {k}")
    if profile["mixing"] not in ("arithmetic", "geometric", "explicit"):
        raise ValueError("Unrecognized Lennard–Jones mixing rule")

    def parameter(key, lo, hi):
        x = profile[key]
        if type(x) not in (int, float) or not math.isfinite(x) or not lo <= x <= hi:
            raise ValueError(f"Invalid profile parameter {key}")

    for k in ("lj14", "coul14"):
        parameter(k, 0, 1)
    parameter("dt_fs", 0.01, 1.0)
    parameter("cutoff_A", 2, 20)
    if packed.get("units") != "angstrom" or packed.get("pbc") != [True, True, False]:
        raise ValueError(
            "Only Ångström orthorhombic XY-periodic, Z-nonperiodic packing is supported"
        )
    atoms = packed["atoms"]
    n = len(atoms)
    if not 2 <= n <= 50_000:
        raise ValueError("Atom count outside the local pilot budget")
    rows, counts, bounds = parse_data(data)
    if counts.get("atoms") != n:
        raise ValueError("Atom count differs from packed job")
    box = packed["box_A"]
    if len(box) != 3 or any(not math.isfinite(v) or v <= 0 for v in box):
        raise ValueError("Invalid packed box")
    if min(box[:2]) <= 2 * profile["cutoff_A"]:
        raise ValueError("Nonbonded cutoff must be below half either periodic edge")
    for axis, L in zip("xyz", box):
        lo, hi = bounds[axis]
        if abs(lo) > 2e-5 or abs(hi - L) > 2e-5:
            raise ValueError("Topology box differs from packed job")
    nt = counts["atom types"]
    if not 1 <= nt <= 128:
        raise ValueError("Atom type count outside supported range")
    masses = _entries(rows["Masses"], nt, "Masses")
    elem = profile["type_elements"]
    if not isinstance(elem, dict) or set(elem) != {str(i) for i in masses}:
        raise ValueError("Missing element mapping for atom types")
    for i, v in masses.items():
        element = elem[str(i)]
        mass = ATOMIC_MASSES.get(element) if isinstance(element, str) else None
        if mass is None or abs(_float(v[0]) / mass - 1) > 0.02:
            raise ValueError(
                "Mass and explicit element identity disagree, or isotope/element is unsupported"
            )
    if profile["mixing"] == "explicit":
        if rows["Pair Coeffs"]:
            raise ValueError(
                "Explicit pair matrix cannot be mixed with diagonal coefficients"
            )
        pairs = {}
        for a, b, ep, sg in rows["PairIJ Coeffs"]:
            key = (_uint(a), _uint(b))
            if key in pairs:
                raise ValueError("Repeated PairIJ coefficient")
            pairs[key] = [ep, sg]
        if set(pairs) != {(i, j) for i in masses for j in range(i, nt + 1)}:
            raise ValueError("Incomplete explicit PairIJ matrix")
    else:
        if rows["PairIJ Coeffs"]:
            raise ValueError("Mixed coefficients require diagonal Pair Coeffs")
        pairs = _entries(rows["Pair Coeffs"], nt, "Pair Coeffs")
    for e, s in pairs.values():
        if _float(e) < 0 or _float(s) < 0 or (_float(e) > 0 and _float(s) == 0):
            raise ValueError("Invalid 12–6 parameters")
    typed = _entries(rows["Atoms"], n, "Atoms")
    by_molecule = defaultdict(list)
    mol_ids = {}
    elements_by_id = {}
    slab = []
    mobile = []
    charges = []
    for i, a in enumerate(atoms, 1):
        mol, type_, q, *xyz = typed[i]
        mol = _uint(mol, False)
        type_ = _uint(type_)
        if a.get("id", i - 1) != i - 1 or elem.get(str(type_)) != a.get("element"):
            raise ValueError(f"Element/atom-order mismatch at id {i}")
        elements_by_id[i] = a["element"]
        if max(abs(_float(v) - a[k]) for k, v in zip("xyz", xyz)) > 0.00201:
            raise ValueError(f"Coordinates differ from packed atom {i}")
        q = _float(q)
        charges.append(q)
        if a["role"] == "slab":
            slab.append(i)
        else:
            mobile.append(i)
            # Each input mobile molecule must retain its membership under a distinct positive numeric ID.
            label = str(a["molecule_id"])
            if mol <= 0 or (label in mol_ids and mol != mol_ids[label]):
                raise ValueError("Changed mobile molecule membership")
            mol_ids[label] = mol
            by_molecule[label].append((q, a.get("formal_charge", 0)))
    if not slab or len(mobile) < 2:
        raise ValueError("A fixed slab and at least two mobile atoms are required")
    if len(set(mol_ids.values())) != len(mol_ids):
        raise ValueError("Different mobile molecules share one numeric ID")
    if abs(math.fsum(charges)) > 1e-5:
        raise ValueError("Topology must be charge neutral within 1e-5 e")
    for label, members in by_molecule.items():
        if abs(math.fsum(x - y for x, y in members)) > 1e-4:
            raise ValueError(
                f"Mobile molecular charge differs from explicit formal charge: {label}"
            )
    seen_bonds = set()
    orders = profile["bond_orders"]
    higher = defaultdict(list)
    for section, (kind, arity) in TERMS.items():
        num = counts.get(kind + " types", 0)
        coeff = _entries(rows[kind.title() + " Coeffs"], num, kind + " coeffs")
        assigned = _entries(rows[section], counts.get(section.lower(), 0), section)
        if section == "Bonds" and (
            not isinstance(orders, dict) or set(orders) != {str(i) for i in coeff}
        ):
            raise ValueError("Bond type to chemical order map is incomplete")
        for constants in coeff.values():
            if kind in ("bond", "angle", "improper") and _float(constants[0]) < 0:
                raise ValueError("Negative harmonic force constant")
            if kind == "bond" and _float(constants[1]) <= 0:
                raise ValueError("Nonpositive harmonic bond equilibrium length")
            if kind == "angle" and not 0 < _float(constants[1]) <= 180:
                raise ValueError("Invalid equilibrium angle")
        for tr in assigned.values():
            ty = _uint(tr[0])
            indices = tuple(_uint(x) for x in tr[1:])
            if (
                ty not in coeff
                or not all(x in typed for x in indices)
                or len(set(indices)) != arity
            ):
                raise ValueError(f"Invalid {kind} indices or missing coefficient")
            if section == "Bonds":
                order = orders[str(ty)]
                if type(order) not in (float, int) or order not in (0.5, 1, 1.5, 2, 3):
                    raise ValueError("Invalid source bond order")
                edge = (*sorted(indices), float(order))
                if edge in seen_bonds:
                    raise ValueError("Duplicate typed bond")
                seen_bonds.add(edge)
            else:
                higher[section].append(indices)
    expected = {(min(i, j) + 1, max(i, j) + 1, float(o)) for i, j, o in packed["bonds"]}
    if seen_bonds != expected:
        raise ValueError(
            "Bond connectivity or chemical order differs from the packing graph"
        )
    edges = {tuple(sorted((i, j))) for i, j, _ in seen_bonds}
    for seq in higher["Angles"] + higher["Dihedrals"]:
        if any(tuple(sorted(pair)) not in edges for pair in zip(seq, seq[1:])):
            raise ValueError("Angle/dihedral does not follow bonded connectivity")
    angles = {tuple(sorted((a, c))) + (b,) for a, b, c in higher["Angles"]}
    neighbors = defaultdict(set)
    for a, b in edges:
        neighbors[a].add(b)
        neighbors[b].add(a)
    torsion_centers = {tuple(sorted((b, c))) for a, b, c, d in higher["Dihedrals"]}
    for b in mobile:
        if atoms[b - 1]["role"] != "polymer":
            continue
        for a, c in combinations(sorted(neighbors[b]), 2):
            if (a, c, b) not in angles:
                raise ValueError("Polymer angle is missing from the typed topology")
        if elements_by_id[b] == "H":
            continue
        for c in neighbors[b]:
            if c <= b or elements_by_id[c] == "H":
                continue
            sides_b = {a for a in neighbors[b] - {c} if elements_by_id[a] != "H"}
            sides_c = {d for d in neighbors[c] - {b} if elements_by_id[d] != "H"}
            if (
                any(a != d for a in sides_b for d in sides_c)
                and (b, c) not in torsion_centers
            ):
                raise ValueError(
                    "No proper torsion defined about an internal polymer heavy-atom bond"
                )
    if water:
        water_sites = defaultdict(list)
        for i, a in enumerate(atoms, 1):
            if a["role"] == "water":
                water_sites[a["molecule_id"]].append(i)
        for ids in water_sites.values():
            oxygen = [i for i in ids if elements_by_id[i] == "O"]
            hydrogens = sorted(i for i in ids if elements_by_id[i] == "H")
            if (
                len(ids) != 3
                or len(oxygen) != 1
                or len(hydrogens) != 2
                or tuple(hydrogens + oxygen) not in angles
            ):
                raise ValueError(
                    "Explicit flexible water requires a defined bonded H–O–H angle"
                )
    return dict(
        atoms=n,
        mobile_atoms=len(mobile),
        fixed_slab_atoms=len(slab),
        slab_ids=slab,
        bonds=len(seen_bonds),
        charge_e=math.fsum(charges),
        has_partial_charges=any(q != 0 for q in charges),
        elements=dict(Counter(elements_by_id.values())),
    )
