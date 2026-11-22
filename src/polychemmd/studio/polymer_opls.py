"""Type the Studio's neutral explicit-H oligomers and compile bonded OPLS terms.

The IDs below refer to the pinned molecule environments in oplsaa2024.lt, not
new fitted force-field parameters. Moltemplate reads the installed parameter
source and raises an error on missing bonded parameters. Numerical parameters
are not bundled in this module.
"""

from __future__ import annotations

from collections import Counter
from importlib import metadata, util
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sysconfig
import tempfile

from .iff_sources import OPLS_SHA, check, fields, sections

FAMILIES = ("PVA", "PEO", "PEtOx")
TERM_STYLES = {
    "Bonds": "bond",
    "Angles": "angle",
    "Dihedrals": "dihedral",
    "Impropers": "improper",
}


def _graph(atoms, bonds):
    elements = [a["element"] for a in atoms]
    neighbours = [dict() for _ in atoms]
    if not atoms or any(a.get("formal_charge", 0) != 0 for a in atoms):
        raise ValueError("OPLS template requires one neutral, explicit-H chain")
    for i, j, order in bonds:
        if (
            type(i) is not int
            or type(j) is not int
            or not 0 <= i < len(atoms)
            or not 0 <= j < len(atoms)
            or i == j
            or j in neighbours[i]
            or type(order) not in (int, float)
            or order not in (1, 2)
        ):
            raise ValueError(
                "Invalid polymer bond index, chemical order or duplicate edge"
            )
        neighbours[i][j] = neighbours[j][i] = int(order)
    valence = {"H": 1, "C": 4, "N": 3, "O": 2}
    if any(sum(ns.values()) != valence.get(e) for e, ns in zip(elements, neighbours)):
        raise ValueError(
            "Incomplete explicit hydrogens, unsupported element or valence"
        )
    seen = {0}
    pending = [0]
    while pending:
        for nb in neighbours[pending.pop()]:
            if nb not in seen:
                seen.add(nb)
                pending.append(nb)
    if len(seen) != len(atoms):
        raise ValueError("Polymer template must be a single connected chain")
    if len(bonds) != len(atoms) - 1:
        raise ValueError("These templates support acyclic polymer chains")
    return elements, neighbours


def source_atom_types(atoms, bonds, family):
    """Assign upstream type IDs using explicit bond orders and local environments."""
    if family not in FAMILIES:
        raise ValueError(f"Unsupported OPLS polymer family {family!r}")
    el, ns = _graph(atoms, bonds)
    signatures = [tuple(sorted(el[j] for j in n)) for n in ns]
    carbonyl = {
        i
        for i, n in enumerate(ns)
        if el[i] == "C" and any(el[j] == "O" and order == 2 for j, order in n.items())
    }
    typed = [None] * len(atoms)
    for i, element in enumerate(el):
        sig = "".join(signatures[i])
        if element == "H":
            continue
        if family == "PVA":
            if element == "O":
                if sig == "C" and next(iter(ns[i])) in carbonyl:
                    typed[i] = 466
                elif sig == "CC" and any(j in carbonyl for j in ns[i]):
                    typed[i] = 467
                elif sig == "CH":
                    typed[i] = 154
            elif element == "C":
                if i in carbonyl and sig == "COO":
                    typed[i] = 465
                elif sig == "CCHO":
                    oxygen = next(j for j in ns[i] if el[j] == "O")
                    typed[i] = 158 if "H" in signatures[oxygen] else 491
                else:
                    typed[i] = {"CHHH": 135, "CCHH": 136}.get(sig)
        elif family == "PEO":
            if element == "O":
                typed[i] = {"CH": 169, "CC": 180}.get(sig)
            elif element == "C" and sig == "CHHO":
                oxygen = next(j for j in ns[i] if el[j] == "O")
                typed[i] = 173 if "H" in signatures[oxygen] else 182
        elif family == "PEtOx":
            if (
                element == "N"
                and sig == "CCC"
                and sum(j in carbonyl for j in ns[i]) == 1
            ):
                typed[i] = 1037
            elif element == "O" and sig == "C" and next(iter(ns[i])) in carbonyl:
                typed[i] = 1036
            elif element == "C":
                if i in carbonyl and sig == "CNO":
                    typed[i] = 1035
                else:
                    typed[i] = {
                        "HHHN": 1038,
                        "CHHN": 1039,
                        "CHHH": 135,
                        "CCHH": 136,
                    }.get(sig)
        if typed[i] is None:
            raise ValueError(
                f"No {family} OPLS type for {element} at index {i} with neighbors {sig}"
            )
    for i, element in enumerate(el):
        if element != "H":
            continue
        parent = next(iter(ns[i]))
        parent_type = typed[parent]
        if (
            family == "PVA"
            and parent_type == 135
            and any(j in carbonyl for j in ns[parent])
        ):
            typed[i] = 282
        else:
            typed[i] = {
                135: 140,
                136: 140,
                154: 155,
                158: 140,
                491: 469,
                169: 170,
                173: 176,
                182: 185,
                1038: 1042,
                1039: 1042,
            }.get(parent_type)
        if typed[i] is None:
            raise ValueError(f"No {family} OPLS hydrogen type on parent {parent_type}")
    count = Counter(typed)
    if family == "PVA" and not (
        count[154] + count[466] > 0
        and count[154] == count[158]
        and count[465] == count[466] == count[467] == count[491]
    ):
        raise ValueError("Expected matched vinyl-alcohol and vinyl-acetate groups")
    if family == "PEO" and not (
        count[169] == 2 and count[173] == 2 and count[182] == 2 * count[180]
    ):
        raise ValueError("Expected a linear PEG diol with two hydroxyl ends")
    if family == "PEtOx" and not (
        count[1037] > 0 and count[1035] == count[1036] == count[1037]
    ):
        raise ValueError("Expected one tertiary propionamide group per PEtOx repeat")
    return typed


def _template(atoms, bonds, type_ids):
    xyz = [tuple(float(a[k]) for k in ("x", "y", "z")) for a in atoms]
    if not all(math.isfinite(x) for row in xyz for x in row):
        raise ValueError("Nonfinite polymer coordinates")
    lines = [
        'import "oplsaa2024.lt"',
        "Chain inherits OPLSAA {",
        ' write("Data Atoms") {',
    ]
    for i, (coords, t) in enumerate(zip(xyz, type_ids), 1):
        lines.append(
            f"  $atom:a{i} $mol:. @atom:{t} 0.0 " + " ".join(f"{v:.9f}" for v in coords)
        )
    lines += [" }", ' write("Data Bond List") {']
    for k, (i, j, _) in enumerate(bonds, 1):
        lines.append(f"  $bond:b{k} $atom:a{i + 1} $atom:a{j + 1}")
    lines += [" }", "}", "chain = new Chain", 'write_once("Data Boundary") {']
    for axis, coords in zip("xyz", zip(*xyz)):
        lines.append(f" {min(coords) - 2:.9f} {max(coords) + 2:.9f} {axis}lo {axis}hi")
    return "\n".join(lines + ["}", ""])


def _compiled(directory, atoms, bonds, source_ids):
    data = sections((directory / "system.data").read_text())
    types = {}
    for row in data["Masses"]:
        fields_, _, comment = row.partition("#")
        idx, mass = fields_.split()
        match = re.match(r"(\d+)_", comment.strip())
        if not match:
            raise ValueError("Moltemplate output lacks OPLS type identity")
        types[int(idx)] = {"source_type": int(match[1]), "mass": float(mass)}
    atom_types = [None] * len(atoms)
    for row in data["Atoms"]:
        idx, _, atom_type = map(int, fields(row)[:3])
        if (
            not 1 <= idx <= len(atoms)
            or atom_types[idx - 1] is not None
            or types[atom_type]["source_type"] != source_ids[idx - 1]
        ):
            raise ValueError("Moltemplate atom identity/order differs from polymer")
        atom_types[idx - 1] = atom_type
        types[atom_type]["element"] = atoms[idx - 1]["element"]
    if None in atom_types:
        raise ValueError("Incomplete Moltemplate atom table")

    by_command = {style + "_coeff": name for name, style in TERM_STYLES.items()}
    coefficients = {kind: {} for kind in TERM_STYLES}
    for row in (directory / "system.in.settings").read_text().splitlines():
        args = fields(row)
        if not args:
            continue
        if args[0] == "pair_coeff" and args[1] == args[2] and len(args) == 5:
            types[int(args[1])].update(epsilon=float(args[3]), sigma=float(args[4]))
        elif (
            args[0] == "set"
            and args[1] == "type"
            and args[3] == "charge"
            and len(args) == 5
        ):
            types[int(args[2])]["charge"] = float(args[4])
        elif args[0] in by_command:
            values = list(map(float, args[2:]))
            if args[0] == "improper_coeff":
                values = [float(args[2]), int(args[3]), int(args[4])]
            coefficients[by_command[args[0]]][int(args[1])] = values
        else:
            raise ValueError("Unsupported Moltemplate setting: " + " ".join(args[:3]))
    for t in types.values():
        if not all(
            key in t and math.isfinite(t[key])
            for key in ("mass", "epsilon", "sigma", "charge")
        ):
            raise ValueError("Incomplete or nonfinite OPLS coefficient")
    if abs(math.fsum(types[t]["charge"] for t in atom_types)) > 1e-7:
        raise ValueError("Assigned OPLS chain partial charges do not sum to zero")
    terms = {}
    for name in TERM_STYLES:
        rows = [list(map(int, fields(line))) for line in data[name]]
        terms[name] = [[t, *[j - 1 for j in ids]] for _, t, *ids in rows]
        if {term[0] for term in terms[name]} != set(coefficients[name]):
            raise ValueError(f"Missing or unused coefficients for {name}")
    # A second check independent of Moltemplate's -checkff. The Studio
    # families are acyclic. Enumerate all expected paths including H ends.
    _, neighbours = _graph(atoms, bonds)
    expected_bonds = {tuple(sorted((i, j))) for i, j, _ in bonds}
    expected_angles = {
        (a, i, b) for i, n in enumerate(neighbours) for a in n for b in n if a < b
    }
    expected_dihedrals = set()
    for i, j in expected_bonds:
        for a in neighbours[i].keys() - {j}:
            for b in neighbours[j].keys() - {i}:
                path = (a, i, j, b)
                expected_dihedrals.add(min(path, path[::-1]))
    for name, expected in (
        ("Bonds", expected_bonds),
        ("Angles", expected_angles),
        ("Dihedrals", expected_dihedrals),
    ):
        paths = [tuple(row[1:]) for row in terms[name]]
        if len(paths) != len(expected) or {min(p, p[::-1]) for p in paths} != expected:
            raise ValueError(f"Compiled {name} differ from complete polymer bond graph")
    return {
        "types": types,
        "atom_types": atom_types,
        "coefficients": coefficients,
        "terms": terms,
    }


def compile_polymer(atoms, bonds, family, directory):
    """Compile one mapped chain; return clean numbering and bonded tables in real units."""
    type_ids = source_atom_types(atoms, bonds, family)
    if family == "PEtOx" and type_ids.count(1037) > 1:
        raise ValueError(
            "Missing sourced NM–CT–CT–NM inter-repeat torsion for PEtOx in "
            "oplsaa2024.lt. Supply a complete independently parameterized topology."
        )
    spec = util.find_spec("moltemplate")
    if spec is None:
        raise RuntimeError("Install moltemplate==2.23.0 (Studio full or studio-models)")
    table = Path(spec.origin).parent / "force_fields" / "oplsaa2024.lt"
    check(table.read_bytes(), OPLS_SHA)
    script_dir = sysconfig.get_path("scripts")
    env = dict(os.environ)
    env["PATH"] = os.pathsep.join([script_dir, env.get("PATH", "")])
    executables = [
        shutil.which(x, path=env["PATH"])
        for x in ("moltemplate.sh", "cleanup_moltemplate.sh")
    ]
    if not all(executables):
        raise RuntimeError("Both Moltemplate POSIX tools are required")
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    template = _template(atoms, bonds, type_ids)
    (directory / "polymer.lt").write_text(template)
    with tempfile.TemporaryDirectory(dir=directory) as tmp:
        work = Path(tmp)
        # Resolve the import to the exact parameter file we fingerprinted.
        (work / "oplsaa2024.lt").symlink_to(table.resolve())
        (work / "system.lt").write_text(template)
        for args in ([executables[0], "-checkff", "system.lt"], [executables[1]]):
            run = subprocess.run(
                args, cwd=work, env=env, capture_output=True, text=True, timeout=180
            )
            if run.returncode:
                raise ValueError(
                    "Moltemplate failed to assign the polymer topology: "
                    + run.stderr[-1800:]
                )
        compiled = _compiled(work, atoms, bonds, type_ids)
    compiled["source"] = {
        "distribution": "moltemplate",
        "version": metadata.version("moltemplate"),
        "table": "oplsaa2024.lt",
        "sha256": OPLS_SHA,
        "doi": "10.1021/acs.jpcb.3c06602",
        "scope": "Explicit neutral linear oligomer; candidate local-group transfer to hydrated interface.",
    }
    return compiled
