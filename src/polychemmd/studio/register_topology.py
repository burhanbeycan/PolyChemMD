"""Operator-only admission of externally parameterised packed cells for bounded MD."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import tempfile

from polychemmd.engines.lammps import LammpsRunSpec, run_lammps
from .topology_file import validate


def _registry():
    return (
        Path(
            os.getenv(
                "POLYCHEMMD_TOPOLOGY_REGISTRY",
                Path.home() / ".polychemmd" / "studio-topologies",
            )
        )
        .expanduser()
        .absolute()
    )


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _encode(d):
    return json.dumps(
        d, sort_keys=True, allow_nan=False, ensure_ascii=False, separators=(",", ":")
    ).encode()


def _read(path, limit=48_000_000):
    if not path.is_file() or path.is_symlink() or path.stat().st_size > limit:
        raise ValueError(f"Missing, nonregular or oversized input: {path.name}")
    return path.read_bytes()


def _id(data, manifest):
    return _sha(_encode({"parameter_data_sha256": _sha(data), "profile": manifest}))


def register(packed_path, data_path, manifest_path):
    packed_bytes = _read(Path(packed_path))
    data = _read(Path(data_path))
    manifest = json.loads(_read(Path(manifest_path), 100_000))
    if not isinstance(manifest, dict) or manifest.get("packed_sha256") != _sha(
        packed_bytes
    ):
        raise ValueError("Packed file checksum differs from the operator manifest")
    validation = validate(json.loads(packed_bytes), data.decode("ascii"), manifest)
    key = _id(data, manifest)
    root = _registry()
    root.mkdir(exist_ok=True, parents=True)
    location = root / key
    if not location.exists():
        with tempfile.TemporaryDirectory(dir=root) as tmp:
            stage = Path(tmp) / "entry"
            stage.mkdir()
            (stage / "profile.json").write_bytes(_encode(manifest))
            (stage / "typed.data").write_bytes(data)
            try:
                stage.rename(location)
            except OSError:
                if not location.is_dir():
                    raise
    if (
        location.is_symlink()
        or _id(
            _read(location / "typed.data"),
            json.loads(_read(location / "profile.json", 100_000)),
        )
        != key
    ):
        raise ValueError("Stored registry entry failed its content hash")
    return {
        "id": key,
        "label": manifest["label"],
        "packed_sha256": manifest["packed_sha256"],
        "atom_count": validation["atoms"],
        "charge_e": validation["charge_e"],
    }


def profiles():
    root = _registry()
    if not root.is_dir():
        return []
    result = []
    for loc in sorted(root.iterdir()):
        if (
            not re.fullmatch("[0-9a-f]{64}", loc.name)
            or loc.is_symlink()
            or not loc.is_dir()
        ):
            continue
        try:
            m = json.loads(_read(loc / "profile.json", 100_000))
            result.append(
                {
                    k: m[k]
                    for k in (
                        "label",
                        "packed_sha256",
                        "parameter_source",
                        "applicability",
                    )
                }
                | {"id": loc.name}
            )
        except (KeyError, OSError, ValueError):
            continue
    return result


def _int(config, k, lo, hi, default):
    v = config.get(k, default)
    if type(v) is not int or not lo <= v <= hi:
        raise ValueError(f"{k} must be an integer from {lo} to {hi}")
    return v


def _ranges(ids):
    start = prev = ids[0]
    groups = []
    for n in ids[1:]:
        if n != prev + 1:
            groups.append(str(start) if start == prev else f"{start}:{prev}")
            start = n
        prev = n
    groups.append(str(start) if start == prev else f"{start}:{prev}")
    return " ".join(groups)


def _script(profile, slab, steps, temp, seed, charged, minimize_steps=0):
    # Numeric model parameters have passed admission. No manifest field is a command.
    mix = profile["mixing"] if profile["mixing"] != "explicit" else "arithmetic"
    lines = [
        "units real",
        "atom_style full",
        "boundary p p f",
        f"pair_style {'lj/cut/coul/long' if charged else 'lj/cut'} {profile['cutoff_A']}",
        "bond_style harmonic",
        "angle_style harmonic",
        "dihedral_style opls",
        "improper_style cvff"
        if profile["schema"] == "studio-fixed-iff-opls-v2"
        else "improper_style harmonic",
        f"special_bonds lj 0 0 {profile['lj14']} coul 0 0 {profile['coul14']}",
        f"pair_modify mix {mix}",
        "read_data typed.data",
        *(["kspace_style pppm 1.0e-5", "kspace_modify slab 3.0"] if charged else []),
        f"timestep {profile['dt_fs']}",
        f"group substrate id {_ranges(slab)}",
        "group mobile subtract all substrate",
        "neighbor 1.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        "variable force_magnitude atom sqrt(fx*fx+fy*fy+fz*fz)",
        "compute fmax mobile reduce max v_force_magnitude",
        f"thermo {max(1, steps // 10)}",
        "thermo_style custom step atoms pe ke etotal c_fmax",
        "thermo_modify lost error flush yes format float %.12g",
    ]
    if minimize_steps:
        # Keep the substrate and periodic cell fixed. The max force reports mobile atoms only.
        lines += [
            "fix immobilize substrate setforce 0.0 0.0 0.0",
            "min_style cg",
            "min_modify dmax 0.05",
            f"minimize 0.0 0.0 {minimize_steps} {minimize_steps * 10}",
            "unfix immobilize",
            "reset_timestep 0",
        ]
    lines += [
        f"dump coords all custom {max(1, steps // 10)} trajectory.dump id type xs ys zs fx fy fz",
        "dump_modify coords sort id format float %.10g",
    ]
    if steps:
        lines.extend(
            [
                f"velocity mobile create {temp} {seed} mom yes rot no dist gaussian",
                "fix integrate mobile nve",
            ]
        )
    lines.extend(
        [
            f"run {steps}",
            "write_dump all custom final.dump id x y z modify sort id format float %.14g",
            'print "POLYCHEMMD_RUN_COMPLETE"',
        ]
    )
    return "\n".join(lines) + "\n"


def _final_structure(packed, dump, step):
    """Join actual end-step positions by LAMMPS ID, independently of trajectory stride."""
    lines = dump.read_text(encoding="ascii").splitlines()
    n = len(packed["atoms"])
    if (
        len(lines) != n + 9
        or lines[:5]
        != [
            "ITEM: TIMESTEP",
            str(step),
            "ITEM: NUMBER OF ATOMS",
            str(n),
            "ITEM: BOX BOUNDS pp pp ff",
        ]
        or lines[8] != "ITEM: ATOMS id x y z"
    ):
        raise ValueError("Incomplete or incompatible final-coordinate snapshot")
    for line, length in zip(lines[5:8], packed["box_A"], strict=True):
        lo, hi = map(float, line.split())
        if not math.isclose(lo, 0.0, abs_tol=2e-5) or not math.isclose(
            hi, length, abs_tol=2e-5
        ):
            raise ValueError("Final orthorhombic box differs from the admitted cell")
    atoms = [None] * n
    for line in lines[9:]:
        atom_id, *coords = line.split()
        i = int(atom_id) - 1
        xyz = list(map(float, coords))
        if (
            not 0 <= i < n
            or atoms[i] is not None
            or len(xyz) != 3
            or not all(map(math.isfinite, xyz))
        ):
            raise ValueError("Invalid or repeated atom in final LAMMPS snapshot")
        if not -1e-8 <= xyz[2] <= packed["box_A"][2] + 1e-8:
            raise ValueError(
                "An atom crossed the nonperiodic z boundary; increase vacuum headspace or revise the pilot"
            )
        atoms[i] = packed["atoms"][i] | dict(zip("xyz", xyz, strict=True))
    return packed | {
        "atoms": atoms,
        "snapshot_step": step,
        "coordinate_source": "final LAMMPS write_dump",
    }


def _groups(packed):
    atoms = packed["atoms"]
    layers = {}
    for i, a in enumerate(atoms, 1):
        if a["role"] == "slab" and a.get("layer") == "middle_Ti":
            layers.setdefault(a["sheet_index"], []).append(i)
    expected = {"single": 1, "slit": 2}.get(packed.get("geometry"))
    if not layers or expected != len(layers):
        return None
    # Select the solvent-exposed side(s) of the actual GUI packing region by z.
    ordered = sorted(
        layers, key=lambda k: sum(atoms[i - 1]["z"] for i in layers[k]) / len(layers[k])
    )
    exposed = {ordered[0]: "top"}
    if len(ordered) == 2:
        exposed[ordered[1]] = "bottom"
    sites = [
        i
        for i, a in enumerate(atoms, 1)
        if a["role"] == "slab"
        and a.get("sheet_index") in exposed
        and a["element"] in ("O", "F")
        and a.get("layer", "")
        in (exposed[a["sheet_index"]] + "_OH", exposed[a["sheet_index"]] + "_F")
    ]
    if not sites:
        return None
    return dict(
        reference_plane=layers[ordered[0]],
        comparison_plane=layers[ordered[1]] if len(ordered) == 2 else [],
        surface_sites=sites,
        water_oxygen=[
            i
            for i, a in enumerate(atoms, 1)
            if a["role"] == "water" and a["element"] == "O"
        ],
        polymer_heavy=[
            i
            for i, a in enumerate(atoms, 1)
            if a["role"] == "polymer" and a["element"] != "H"
        ],
        plane_restraint="fixed",
    )


def run_registered(config, out, source_dir):
    allowed = {
        "source_job",
        "source_job_id",
        "topology_id",
        "profile",
        "md_profile",
        "temperature_k",
        "minimize_steps",
        "seed",
        "steps",
    }
    if set(config) - allowed:
        raise ValueError(
            "Unknown registered pilot input fields: "
            + ", ".join(sorted(set(config) - allowed))
        )
    if any(config[k] != "registered" for k in ("profile", "md_profile") if k in config):
        raise ValueError("Wrong profile for registered topology")
    source_dir = Path(source_dir)
    out = Path(out)
    key = config.get("topology_id", "")
    if not isinstance(key, str) or not re.fullmatch("[0-9a-f]{64}", key):
        raise ValueError("Select a local registered topology ID")
    loc = _registry() / key
    if loc.is_symlink():
        raise ValueError("Registry entries must be local regular directories")
    data = _read(loc / "typed.data")
    p = json.loads(_read(loc / "profile.json", 100_000))
    packed_bytes = _read(source_dir / "packed.json")
    packed = json.loads(packed_bytes)
    if _id(data, p) != key or _sha(packed_bytes) != p.get("packed_sha256"):
        raise ValueError("Registry content or selected source-job checksum mismatch")
    audited = validate(packed, data.decode("ascii"), p)
    steps = _int(config, "steps", 0, 200, 0)
    minimize_steps = _int(config, "minimize_steps", 0, 100, 0)
    seed = _int(config, "seed", 1, 900_000_000, 17)
    temp = config.get("temperature_k", 298.15)
    if (
        type(temp) not in (int, float)
        or not math.isfinite(temp)
        or not 1 <= temp <= 1000
    ):
        raise ValueError("Initial mobile temperature must be within [1,1000] K")
    out.mkdir(parents=True, exist_ok=True)
    script = _script(
        p,
        audited["slab_ids"],
        steps,
        temp,
        seed,
        audited["has_partial_charges"],
        minimize_steps,
    )
    result = run_lammps(
        LammpsRunSpec(
            input_script=script,
            run_directory=out / "engine",
            required_thermo_columns=("Step", "Atoms", "PotEng", "KinEng", "TotEng"),
            timeout_seconds=120,
            input_files=(("typed.data", data.decode("ascii")),),
            units_style="real",
            expected_artifacts=("log.lammps", "trajectory.dump", "final.dump"),
        )
    )
    record = result.to_dict()
    thermo = record["thermo_blocks"][-1]
    idx = thermo["columns"].index("PotEng")
    force = thermo["columns"].index("c_fmax")
    first = record["thermo_blocks"][0]["rows"][0]
    # Nonperiodic z and fixed basal layers are part of the recipe, including run 0.
    end_dump = out / "engine" / "final.dump"
    final = _final_structure(packed, end_dump, steps)
    embedded = {
        "profile": "registered",
        "final_snapshot_sha256": _sha(end_dump.read_bytes()),
        "final_snapshot_step": steps,
        "topology_id": key,
        "label": p["label"],
        "source_doi_or_reference": p["parameter_source"],
        "source_job": config.get("source_job", config.get("source_job_id")),
        "packed_sha256": p["packed_sha256"],
        "parameter_data_sha256": _sha(data),
        "atoms": audited["atoms"],
        "charge_e": audited["charge_e"],
        "steps": steps,
        "dt_fs": p["dt_fs"],
        "initial_kcal_mol": first[idx],
        "minimize_steps_requested": minimize_steps,
        "minimize_iterations": int(record["thermo_blocks"][0]["rows"][-1][0])
        if minimize_steps
        else 0,
        "post_minimize_kcal_mol": thermo["rows"][0][idx],
        "initial_mobile_fmax_kcal_mol_A": first[force],
        "post_minimize_mobile_fmax_kcal_mol_A": thermo["rows"][0][force],
        "final_mobile_fmax_kcal_mol_A": thermo["rows"][-1][force],
        "final_kcal_mol": thermo["rows"][-1][idx],
        "dump_force_unit": "kcal mol^-1 angstrom^-1",
        "ensemble": "fixed cell and fixed substrate; initialized mobile NVE pilot"
        if steps
        else "static potential energy",
        "scope": "Source-matched topology and engine execution. Equilibrium properties and adsorption remain unvalidated.",
    }
    (out / "registered.json").write_bytes(_encode(embedded))
    (out / "structure.json").write_bytes(_encode(final))
    shutil.copyfile(out / "engine" / "log.lammps", out / "registered-lammps.log")
    shutil.copyfile(out / "engine" / "run_record.json", out / "registered-engine.json")
    native_dump = out / "engine" / "trajectory.dump"
    # Keep the native trace in the engine run and expose the same bytes as a GUI download.
    try:
        os.link(native_dump, out / "registered-trajectory.dump")
    except OSError:
        shutil.copyfile(native_dump, out / "registered-trajectory.dump")
    (out / "registered.in").write_text(script, encoding="ascii")
    (out / "registered-profile.json").write_bytes(_encode(p))
    groups = _groups(packed)
    if groups is not None:
        (out / "trajectory-groups.json").write_bytes(_encode(groups))
    return embedded


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    sub = cli.add_subparsers(dest="action", required=True)
    add = sub.add_parser("add")
    add.add_argument("--packed", type=Path, required=True)
    add.add_argument("--data", type=Path, required=True)
    add.add_argument("--manifest", type=Path, required=True)
    sub.add_parser("list")
    pilot = sub.add_parser(
        "pilot", help="Local bounded verification of an already registered topology"
    )
    pilot.add_argument("--source-dir", type=Path, required=True)
    pilot.add_argument("--id", required=True)
    pilot.add_argument("--out", type=Path, required=True)
    pilot.add_argument("--steps", type=int, default=0)
    pilot.add_argument("--seed", type=int, default=17)
    pilot.add_argument("--minimize-steps", type=int, default=0)
    pilot.add_argument("--temperature-k", type=float, default=298.15)
    a = cli.parse_args()
    if a.action == "add":
        result = register(a.packed, a.data, a.manifest)
    elif a.action == "pilot":
        result = run_registered(
            dict(
                topology_id=a.id,
                steps=a.steps,
                minimize_steps=a.minimize_steps,
                seed=a.seed,
                temperature_k=a.temperature_k,
            ),
            a.out,
            a.source_dir,
        )
    else:
        result = profiles()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
