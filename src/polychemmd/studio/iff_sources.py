"""Prepare source-pinned explicit-water Ti3C2Tx experiments in a local cache.

The studio never downloads force-field data while executing an MD job. Use
this module once to prepare published inputs, or import the existing Atlas.
"""

import argparse
from collections import defaultdict
import gzip
import hashlib
from importlib import resources
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request
from zipfile import ZipFile

DOI = "10.1021/acsnano.5c20046"
URL = "https://pmc-oa-opendata.s3.amazonaws.com/PMC13450432.1/nn5c20046_si_002.zip"
SI_SHA = "ca1e4c50ee2b154a58cd827f2f31e43c284aabd2cb39672cdc92d1621b0827da"
NATIVE_MEMBER = (
    "SUPPORTING FILES/Run Scripts and Examples in LAMMPS and NAMD/"
    "LAMMPS MXene liquid interfaces/Stack3x3_H2Osys_OH_47.data"
)
HASHES = {
    "iff-water": "ec194edb8b6a5d1923269a1dee4bdf2ba047964ad930a5d52a887baa112fa873",
    "iff-opls-clustered": "c86fc1f129fe8ccc8171a1e5f15639872a56adc6c5aa9884a222bc03fec999b8",
    "iff-opls-dispersed": "62f4e41d8b2a8b1069c7cff4e5dc828cf833212d1e89e0e9bbfa4dde18a13873",
    "iff-charmm-clustered": "9649fbfdc72d05064e6e8ce70efc88a573823e028a44180b6cbcf9b7a76666da",
    "iff-charmm-dispersed": "b324480ab52e2d0af7aeb1b784ce8d8559c5effb13b10d73359917752b0b34cb",
}
PAIR_SHA = "b572e243a7287850887f67bbc0f46da1efb58986b2abb75530ff52dfcb114b09"
OPLS_SHA = "77387af74ff10844d89b1d69785b7d89d90a9a5ffe5c25492f9827b9f211225d"
SEQUENCES = ("clustered", "dispersed")
ASSETS = resources.files(__package__).joinpath("md_assets")


def cache_directory():
    override = os.getenv("POLYCHEMMD_IFF_CACHE")
    return (
        Path(override).expanduser()
        if override
        else Path.home() / ".cache" / "polychemmd" / "iff"
    )


def check(content, expected):
    if hashlib.sha256(content).hexdigest() != expected:
        raise ValueError(
            "SHA-256 mismatch for published/calibrated input; refusing to run"
        )


def install(cache, name, data, digest):
    check(data, digest)
    cache.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(dir=cache, delete=False) as handle:
        tmp = Path(handle.name)
        try:
            handle.write(data)
            handle.flush()
        except BaseException:
            tmp.unlink(missing_ok=True)
            raise
    try:
        tmp.replace(cache / name)
    finally:
        tmp.unlink(missing_ok=True)


def load(profile, sequence="clustered", cache=None):
    key = profile if profile == "iff-water" else f"{profile}-{sequence}"
    if key not in HASHES:
        raise ValueError("Unknown IFF profile or sequence")
    cache = cache or cache_directory()
    try:
        data = (cache / (key + ".data")).read_bytes()
        check(data, HASHES[key])
        pairs = None
        if profile == "iff-charmm":
            pairs = (cache / "fixed_iff_pairs.in").read_bytes()
            check(pairs, PAIR_SHA)
    except FileNotFoundError as exc:
        raise RuntimeError(
            f"{key} is not prepared. Run python -m polychemmd.studio.iff_sources "
            "--fetch --opls or import your existing Atlas with --atlas FILE."
        ) from exc
    return data, pairs


def status(cache=None):
    cache = cache or cache_directory()
    return {
        key: (cache / (key + ".data")).is_file()
        and ("charmm" not in key or (cache / "fixed_iff_pairs.in").is_file())
        for key in HASHES
    }


def import_si(path, cache):
    h = hashlib.sha256()
    with Path(path).open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            h.update(block)
    if h.hexdigest() != SI_SHA:
        raise ValueError("Unexpected Supporting Information archive SHA-256")
    with ZipFile(path) as source:
        install(
            cache, "iff-water.data", source.read(NATIVE_MEMBER), HASHES["iff-water"]
        )


def import_atlas(path, cache):
    def mixed(source):
        for profile, folder in (
            ("iff-opls", "models"),
            ("iff-charmm", "models_charmm"),
        ):
            for sequence in SEQUENCES:
                prefix = f"study/interfaces/{folder}/{sequence}/"
                name = f"{profile}-{sequence}"
                compressed = source.read(prefix + "interface.data.gz")
                with gzip.GzipFile(fileobj=io.BytesIO(compressed)) as gz:
                    data = gz.read(12 << 20)
                install(cache, name + ".data", data, HASHES[name])
                if profile == "iff-charmm":
                    install(
                        cache,
                        "fixed_iff_pairs.in",
                        source.read(prefix + "fixed_iff_pairs.in"),
                        PAIR_SHA,
                    )

    with ZipFile(path) as source:
        if "previous/PVA_case_study.zip" in source.namelist():
            if source.getinfo("previous/PVA_case_study.zip").file_size > 100 << 20:
                raise ValueError("Unexpected nested PVA archive size")
            with ZipFile(
                io.BytesIO(source.read("previous/PVA_case_study.zip"))
            ) as nested:
                mixed(nested)
        else:
            mixed(source)


PARAMS = [
    "Masses",
    "Pair Coeffs",
    "Bond Coeffs",
    "Angle Coeffs",
    "Dihedral Coeffs",
    "Improper Coeffs",
]
TABLES = ["Atoms", "Bonds", "Angles", "Dihedrals", "Impropers"]


def sections(text):
    parts = defaultdict(list)
    section = "header"
    for line in text.splitlines():
        prefix = line.partition("#")[0].strip()
        if prefix in PARAMS + TABLES:
            section = prefix
        elif prefix:
            parts[section].append(line.strip())
    return parts


def fields(row):
    return row.partition("#")[0].split()


def assemble(native, directory, sequence):
    """Merge typed oligomer and unchanged native terms with a pinned placement."""
    s = sections(native.decode())
    p = sections((directory / "system.data").read_text())
    plan = json.loads(ASSETS.joinpath(f"placement_{sequence}.json").read_text())
    settings = [
        fields(r)
        for r in (directory / "system.in.settings").read_text().splitlines()
        if r.strip()
    ]
    offsets = {k: len(s[k]) for k in PARAMS}
    charges = {int(r[2]): float(r[4]) for r in settings if r[0] == "set"}
    for row in settings:
        if row[0].endswith("_coeff"):
            kind = row[0].split("_")[0].capitalize() + " Coeffs"
            tokens = row[3:] if kind == "Pair Coeffs" else row[2:]
            if kind == "Dihedral Coeffs":
                tokens = ["opls", *tokens]
            p[kind].append(" ".join([str(int(row[1]) + offsets[kind]), *tokens]))
    merged = {}
    merged["Masses"] = s["Masses"] + [
        f"{int(fields(row)[0]) + offsets['Masses']} " + " ".join(row.split()[1:])
        for row in p["Masses"]
    ]
    for kind in PARAMS[1:]:
        src = s[kind]
        if kind == "Dihedral Coeffs":
            src = [r.split()[0] + " harmonic " + " ".join(r.split()[1:]) for r in src]
        merged[kind] = src + p[kind]
    oxygen = set(plan["removed_water_oxygen_ids"])
    removed = set(oxygen)
    for row in s["Bonds"]:
        _, kind, a, b = map(int, fields(row))
        if kind == 10 and (a in oxygen or b in oxygen):
            removed.update([a, b])
    if len(oxygen) != 70 or len(removed) != 210:
        raise ValueError("Expected 70 complete three-site waters")
    native_rows = [fields(r) for r in s["Atoms"] if int(fields(r)[0]) not in removed]
    atom_map = {int(a[0]): i for i, a in enumerate(native_rows, 1)}
    merged["Atoms"] = [" ".join([str(i), *a[1:]]) for i, a in enumerate(native_rows, 1)]
    for row, xyz in zip(p["Atoms"], plan["placed_polymer_xyz_A"], strict=True):
        a = fields(row)
        merged["Atoms"].append(
            " ".join(
                [
                    str(len(native_rows) + int(a[0])),
                    "99000",
                    str(int(a[2]) + offsets["Masses"]),
                    str(charges[int(a[2])]),
                    *[f"{v:.8f}" for v in xyz],
                    "0",
                    "0",
                    "0",
                ]
            )
        )
    for table, param in zip(TABLES[1:], PARAMS[2:]):
        rows = []
        for row in s[table]:
            _, kind, *ids = map(int, fields(row))
            if not any(a in removed for a in ids):
                rows.append(
                    " ".join(
                        map(str, [len(rows) + 1, kind, *[atom_map[a] for a in ids]])
                    )
                )
        for row in p[table]:
            _, kind, *ids = map(int, fields(row))
            rows.append(
                " ".join(
                    map(
                        str,
                        [
                            len(rows) + 1,
                            kind + offsets[param],
                            *[a + len(native_rows) for a in ids],
                        ],
                    )
                )
            )
        merged[table] = rows
    lines = [
        "Fixed substrate, flexible native water, OPLS-AA 2024 PVA/VAc; candidate assembly",
        "",
    ]
    lines += [f"{len(merged[k])} {k.lower()}" for k in TABLES] + [""]
    for kind, param in zip(
        ("atom", "bond", "angle", "dihedral", "improper"), ["Masses", *PARAMS[2:]]
    ):
        lines.append(f"{len(merged[param])} {kind} types")
    lines += [
        "",
        *[
            s
            for s in s["header"]
            if s.endswith(("xlo xhi", "ylo yhi", "zlo zhi", "xy xz yz"))
        ],
    ]
    for kind in PARAMS + TABLES:
        lines += ["", kind, "", *merged[kind]]
    return ("\n".join(lines) + "\n").encode()


def rebuild_opls(cache):
    import importlib.util

    spec = importlib.util.find_spec("moltemplate")
    if spec is None:
        raise RuntimeError(
            "Install moltemplate==2.23.0 and activate the environment containing its executables"
        )
    source = Path(spec.origin).parent / "force_fields" / "oplsaa2024.lt"
    check(source.read_bytes(), OPLS_SHA)
    native, _ = load("iff-water", cache=cache)
    commands = [shutil.which(x) for x in ("moltemplate.sh", "cleanup_moltemplate.sh")]
    if not all(commands):
        raise RuntimeError(
            "Activate moltemplate==2.23.0 with both shell utilities on PATH"
        )
    for sequence in SEQUENCES:
        with tempfile.TemporaryDirectory() as td:
            directory = Path(td)
            (directory / "system.lt").write_bytes(
                ASSETS.joinpath(f"pva_{sequence}.lt").read_bytes()
            )
            for args in ([commands[0], "system.lt"], [commands[1]]):
                subprocess.run(
                    args, cwd=directory, check=True, timeout=180, capture_output=True
                )
            key = "iff-opls-" + sequence
            install(
                cache, key + ".data", assemble(native, directory, sequence), HASHES[key]
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--fetch",
        action="store_true",
        help="Download the 114 MB public Supporting Information",
    )
    group.add_argument(
        "--si", type=Path, help="Use the original local Supporting Information zip"
    )
    parser.add_argument(
        "--atlas",
        type=Path,
        help="Import pinned polymer profiles from the existing two Atlas ZIP layouts",
    )
    parser.add_argument(
        "--opls",
        action="store_true",
        help="Rebuild both OPLS polymer sequences from the reference IFF source",
    )
    parser.add_argument(
        "--status",
        action="store_true",
        help="Print cached profiles (digests rechecked before use)",
    )
    args = parser.parse_args()
    cache = cache_directory()
    if args.fetch:
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "supporting.zip"
            with urllib.request.urlopen(URL, timeout=90) as r, path.open("wb") as f:
                shutil.copyfileobj(r, f)
            import_si(path, cache)
    if args.si:
        import_si(args.si, cache)
    if args.opls:
        rebuild_opls(cache)
    if args.atlas:
        import_atlas(args.atlas, cache)
    print(json.dumps({"cache": str(cache), "profiles": status(cache)}, indent=2))


if __name__ == "__main__":
    main()
