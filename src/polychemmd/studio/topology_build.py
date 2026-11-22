"""Compile a packed OH/F Ti3C2 film into a registered fixed-field topology.

Native nonbonded parameters and water are read from a verified IFF source.
Organic valence and nonbonded terms come from a local OPLS/Moltemplate build.
Only the mobile coordinates evolve. Length is Å, energy is kcal/mol in LAMMPS.
"""

from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
import math
from pathlib import Path

from . import iff_sources
from .register_topology import register
from .topology_file import TERMS


LAYER_TYPES = {
    ("middle_Ti", "Ti"): 6,
    ("Ti", "Ti"): 1,
    ("C", "C"): 2,  # The four native carbon environments have identical LJ and charge.
    ("F", "F"): 7,
    ("OH", "O"): 8,
    ("OH", "H"): 9,
}
WATER_TYPES = {"O": 10, "H": 11}
NATIVE_ELEMENTS = {
    1: "Ti",
    2: "C",
    6: "Ti",
    7: "F",
    8: "O",
    9: "H",
    **{10: "O", 11: "H"},
}


def _digest(path):
    return sha256(path.read_bytes()).hexdigest()


def _native_parameters():
    data, _ = iff_sources.load("iff-water")
    # The valence tables of the fixed crystal are not used. Read only the section
    # prefix through the end of Atoms instead of retaining all solid torsions.
    text = data.decode("ascii").split("\nBonds", 1)[0]
    sections = iff_sources.sections(text)

    def table(name):
        return {
            int(t[0]): list(map(float, t[1:]))
            for s in sections[name]
            if (t := iff_sources.fields(s))
        }

    masses = table("Masses")
    pair = table("Pair Coeffs")
    charges = defaultdict(set)
    for row in sections["Atoms"]:
        _, _, atom_type, charge, *_ = iff_sources.fields(row)
        charges[int(atom_type)].add(float(charge))
    if any(len(q) != 1 for q in charges.values()):
        raise ValueError("Nonunique native type charges")
    for k in (3, 4, 5):
        if (pair[k], masses[k], charges[k]) != (pair[2], masses[2], charges[2]):
            raise ValueError(
                "Rigid-field carbon type reduction is invalid for this source"
            )
    types = {
        k: dict(
            source="IFF",
            source_type=k,
            element=el,
            mass=masses[k][0],
            epsilon=pair[k][0],
            sigma=pair[k][1],
            charge=next(iter(charges[k])),
        )
        for k, el in NATIVE_ELEMENTS.items()
    }
    return types, table("Bond Coeffs")[10], table("Angle Coeffs")[29]


def _native_key(a):
    if a["role"] == "water":
        return WATER_TYPES[a["element"]]
    layer = a.get("layer", "")
    if layer.startswith(("top_", "bottom_")):
        layer = layer.split("_", 1)[1]
    try:
        return LAYER_TYPES[(layer, a["element"])]
    except KeyError as exc:
        raise ValueError(
            f"Unrecognised fixed-slab site or chemical element: {layer}"
        ) from exc


def _audit_sheet(atoms, keys):
    sheets = defaultdict(list)
    for a, k in zip(atoms, keys, strict=True):
        if a["role"] == "slab":
            sheets[a["sheet_index"]].append((a, k))
    if not sheets:
        raise ValueError("Expected one or more chemically identified Ti3C2 sheets")
    for pairs in sheets.values():
        counts = Counter(k for _, k in pairs)
        n = counts[6]
        if not n or counts[1] != 2 * n or counts[2] != 2 * n:
            raise ValueError("Expected Ti3C2 basal stoichiometry")
        for face in ("top", "bottom"):
            sites = defaultdict(list)
            for a, k in pairs:
                if k in (7, 8, 9) and a["layer"].startswith(face + "_"):
                    sites[a["site"]].append(k)
            if len(sites) != n or any(
                sorted(s) not in ([7], [8, 9]) for s in sites.values()
            ):
                raise ValueError(
                    "Each face needs one F or complete OH at every basal site"
                )


def _molecule_bonds(packed):
    atoms = packed["atoms"]
    grouped = defaultdict(list)
    for a, b, order in packed["bonds"]:
        m = atoms[a]["molecule_id"]
        if m != atoms[b]["molecule_id"]:
            raise ValueError("Intermolecular covalent links are unsupported")
        grouped[m].append((a, b, order))
    return grouped


def _assemble(packed, shape, directory):
    native, oh, hoh = _native_parameters()
    atoms = packed["atoms"]
    bonds_by_mol = _molecule_bonds(packed)
    family = shape["polymer"]["name"]
    groups = defaultdict(list)
    labels = []
    for i, a in enumerate(atoms):
        if a["role"] not in ("slab", "polymer", "water"):
            raise ValueError(
                "Salt/ions or other species have no matched parameters in this model. Pack at zero salt."
            )
        groups[a["molecule_id"]].append(i)
        labels.append(_native_key(a) if a["role"] != "polymer" else None)
    _audit_sheet(atoms, labels)
    all_types = {("IFF", i): t for i, t in native.items()}
    labels = [("IFF", k) if k is not None else None for k in labels]
    coefficients = {k: {} for k in TERMS}
    terms = {k: [] for k in TERMS}
    bond_orders = {}
    had_polymer = None
    chain_count = 0
    water_count = 0
    for name, ids in groups.items():
        role = atoms[ids[0]]["role"]
        if role == "slab":
            if bonds_by_mol[name]:
                raise ValueError(
                    "Expected a fixed substrate without internal bonded terms"
                )
            continue
        if role == "water":
            oxygen = [i for i in ids if atoms[i]["element"] == "O"]
            hydrogens = [i for i in ids if atoms[i]["element"] == "H"]
            if len(ids) != 3 or len(oxygen) != 1 or len(hydrogens) != 2:
                raise ValueError(
                    "Only complete three-site native flexible waters are supported"
                )
            o = oxygen[0]
            h1, h2 = hydrogens
            coefficients["Bonds"]["water"] = oh
            coefficients["Angles"]["water"] = hoh
            bond_orders["water"] = 1
            terms["Bonds"] += [("water", o, h1), ("water", o, h2)]
            terms["Angles"].append(("water", h1, o, h2))
            water_count += 1
            continue
        chain_count += 1
        template = [atoms[i] for i in ids]
        local = {i: j for j, i in enumerate(ids)}
        bonds = [(local[a], local[b], v) for a, b, v in bonds_by_mol[name]]
        identity = [
            (a["element"], a["atom_map"], a["atom_id"], a.get("formal_charge", 0))
            for a in template
        ]
        if identity != [
            (a["element"], a["atom_map"], a["atom_id"], a.get("formal_charge", 0))
            for a in shape["templates"]["polymer"]
        ] or (
            {(*sorted((i, j)), o) for i, j, o in bonds}
            != {(*sorted((i, j)), o) for i, j, o in shape["polymer_bonds"]}
        ):
            raise ValueError(
                "Packed chain does not reproduce the source template atom identity and bonds"
            )
        if had_polymer is None:
            from .polymer_opls import compile_polymer

            had_polymer = compile_polymer(template, bonds, family, directory)
            for i, t in had_polymer["types"].items():
                all_types[("OPLS", i)] = dict(t, source="OPLS")
            for section in TERMS:
                coefficients[section].update(had_polymer["coefficients"][section])
        for i, k in zip(ids, had_polymer["atom_types"], strict=True):
            labels[i] = ("OPLS", k)
        order_by_edge = {tuple(sorted((a, b))): v for a, b, v in bonds}
        for section in TERMS:
            for kind, *vs in had_polymer["terms"][section]:
                if section == "Bonds":
                    value = order_by_edge[tuple(sorted(vs))]
                    if kind in bond_orders and bond_orders[kind] != value:
                        raise ValueError(
                            "Inconsistent bond orders under a single OPLS coefficient"
                        )
                    bond_orders[kind] = value
                terms[section].append((kind, *[ids[i] for i in vs]))

    if chain_count != shape["requested_chains"] or not (water_count or chain_count):
        raise ValueError(
            "Chain count differs from geometry or cell contains no mobile molecules"
        )
    return all_types, labels, coefficients, terms, bond_orders, chain_count, water_count


def _write_data(path, packed, all_types, labels, coefficients, terms, orders):
    keys = sorted(set(labels))
    type_ids = {k: i for i, k in enumerate(keys, 1)}
    term_ids = {
        section: {k: i for i, k in enumerate(spec, 1)}
        for section, spec in coefficients.items()
    }

    def line(*tokens):
        return " ".join(
            format(x, ".13g") if isinstance(x, float) else str(x) for x in tokens
        )

    content = [
        "Rigid OH/F Ti3C2T2 surface field; explicit fluid coordinates from Studio packing",
        "",
        f"{len(labels)} atoms",
    ]
    content += [f"{len(terms[s])} {s.lower()}" for s in TERMS]
    content += ["", f"{len(keys)} atom types"]
    content += [
        f"{len(coefficients[s])} {kind} types" for s, (kind, _) in TERMS.items()
    ]
    content += [
        "",
        *[line(0.0, v, d + "lo", d + "hi") for d, v in zip("xyz", packed["box_A"])],
    ]
    content += [
        "",
        "Masses",
        "",
        *[line(type_ids[k], all_types[k]["mass"]) for k in keys],
    ]
    content += ["", "PairIJ Coeffs", ""]
    for i, k in enumerate(keys):
        t = all_types[k]
        for h in keys[i:]:
            u = all_types[h]
            fixed_pair = k[0] == h[0] == "IFF" and k[1] < 10 and h[1] < 10
            # Only the external LJ field of fixed sites acts on the fluid.
            epsilon = 0.0 if fixed_pair else math.sqrt(t["epsilon"] * u["epsilon"])
            sigma = math.sqrt(t["sigma"] * u["sigma"])
            content.append(line(type_ids[k], type_ids[h], epsilon, sigma))
    for s, (kind, _) in TERMS.items():
        if coefficients[s]:
            content.extend(["", kind.title() + " Coeffs", ""])
            for k, coeff in coefficients[s].items():
                content.append(line(term_ids[s][k], *coeff))
    content.extend(["", "Atoms # full", ""])
    molecule_id = {}
    for i, (a, k) in enumerate(zip(packed["atoms"], labels, strict=True), 1):
        m = molecule_id.setdefault(a["molecule_id"], len(molecule_id) + 1)
        content.append(
            line(i, m, type_ids[k], all_types[k]["charge"], a["x"], a["y"], a["z"])
        )
    for s in TERMS:
        if terms[s]:
            content.extend(["", s, ""])
            content.extend(
                line(i, term_ids[s][k], *[a + 1 for a in ids])
                for i, (k, *ids) in enumerate(terms[s], 1)
            )
    path.write_text("\n".join(content) + "\n", encoding="ascii")
    return (
        {str(i): all_types[k]["element"] for k, i in type_ids.items()},
        {str(term_ids["Bonds"][k]): o for k, o in orders.items()},
        type_ids,
    )


def run_topology(config, out, *, source_dir):
    """Compile and register one exact packing job. No physical solver is run here."""
    unexpected = set(config) - {"source_job", "source_job_id", "cutoff_A"}
    if unexpected:
        raise ValueError(f"Unknown topology inputs: {sorted(unexpected)}")
    source_dir, out = Path(source_dir), Path(out)
    cutoff = config.get("cutoff_A", 8.0)
    if (
        type(cutoff) not in (float, int)
        or not math.isfinite(cutoff)
        or not 6 <= cutoff <= 12
    ):
        raise ValueError("Choose a finite real-space cutoff between 6 and 12 Å")
    pp = source_dir / "packed.json"
    sp = source_dir / "structure.json"
    shape = json.loads(sp.read_bytes())
    packing = json.loads((source_dir / "packing_report.json").read_bytes())
    if packing["geometry_input_sha256"] != _digest(sp):
        raise ValueError("Packing provenance does not match the source geometry")
    packed = json.loads(pp.read_bytes())
    if packed["pbc"] != [True, True, False] or packed["units"] != "angstrom":
        raise ValueError(
            "Only XY-periodic, orthorhombic, open-z Ångström cells are supported"
        )
    if min(packed["box_A"][:2]) <= 2 * cutoff:
        raise ValueError(
            "Both lateral cell lengths must exceed twice the nonbonded cutoff"
        )
    if (
        packing["salt_M_requested"] != 0
        or packing["salt_M_realized"] != 0
        or packing["ion_pairs"] != 0
        or any(a["role"] in ("sodium", "chloride") for a in packed["atoms"])
    ):
        raise ValueError(
            "This native flexible water model has no qualified NaCl parameters. Pack with salt_M=0."
        )
    family = shape["polymer"]["name"]
    if family not in ("bare", "PVA", "PEO", "PEtOx"):
        raise ValueError(f"No documented organic parameter assignment for {family}")
    out.mkdir(parents=True, exist_ok=True)
    types, labels, coefficients, terms, orders, chains, waters = _assemble(
        packed, shape, out / "oligomer"
    )
    typed_path = out / "typed.data"
    elements, bond_order, type_ids = _write_data(
        typed_path, packed, types, labels, coefficients, terms, orders
    )
    packed_sha = _digest(pp)
    profile = dict(
        schema="studio-fixed-iff-opls-v2",
        label=f"Generated Ti3C2T2 / {family} fixed-field candidate",
        parameter_source="Armstrong et al. ACS Nano 10.1021/acsnano.5c20046; polymer OPLS/2020 10.1021/acs.jpcb.3c06602 via Moltemplate 2.23.0 (for polymer-containing cells).",
        applicability="Frozen, fully OH/F-terminated Ti3C2 sheets with explicit fixed charges; mobile neutral fluid. Geometric cross LJ and slab-PPPM electrostatics. Energies contain a constant substrate self contribution.",
        water_model="Flexible deposited IFF liquid-interface water"
        if waters
        else "none",
        type_elements=elements,
        bond_orders=bond_order,
        mixing="explicit",
        lj14=0.5,
        coul14=0.5,
        cutoff_A=cutoff,
        dt_fs=0.25,
        packed_sha256=packed_sha,
        water_constraints="flexible" if waters else "no_water",
    )
    profile_path = out / "profile.json"
    profile_path.write_text(json.dumps(profile, indent=2) + "\n")
    mapping_path = out / "atom-types.csv"
    with mapping_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "atom_id_1based",
                "molecule",
                "atom_key",
                "element",
                "role",
                "layer",
                "atom_type",
                "parameter_source",
                "source_type",
                "charge_e",
            ]
        )
        for i, (a, k) in enumerate(zip(packed["atoms"], labels, strict=True), 1):
            t = types[k]
            w.writerow(
                [
                    i,
                    a["molecule_id"],
                    a.get("atom_id", ""),
                    a["element"],
                    a["role"],
                    a.get("layer", ""),
                    type_ids[k],
                    t["source"],
                    t["source_type"],
                    t["charge"],
                ]
            )
    record = register(pp, typed_path, profile_path)
    result = dict(
        profile="generated-fixed-iff-opls",
        topology_id=record["id"],
        source_job=config.get("source_job", config.get("source_job_id")),
        atoms=record["atom_count"],
        charge_e=record["charge_e"],
        polymer=family,
        chain_count=chains,
        water_count=waters,
        salt_M_requested=packing["salt_M_requested"],
        salt_M_realized=packing["salt_M_realized"],
        ion_pairs=packing["ion_pairs"],
        packed_sha256=packed_sha,
        geometry_sha256=_digest(sp),
        topology_sha256=_digest(typed_path),
        map_sha256=_digest(mapping_path),
        term_counts={s: len(v) for s, v in terms.items()},
        cutoff_A=cutoff,
        parameter_provenance={
            "native_sha256": iff_sources.HASHES["iff-water"],
            "native_doi": iff_sources.DOI,
            "opls_sha256": iff_sources.OPLS_SHA if chains else None,
            "polymer": shape["polymer"],
            "face_occupancy": shape["face_occupancy"],
        },
        applicability=profile["applicability"],
        validation="Packed coordinates, molecular graphs, partial-charge balances and coefficient coverage verified; adsorption remains uncalibrated.",
    )
    (out / "topology-build.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n"
    )
    return result
