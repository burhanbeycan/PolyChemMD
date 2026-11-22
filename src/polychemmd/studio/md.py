"""Source-pinned, finite MD of water and two PVA/VAc sequences on Ti3C2(OH,F).

This explicit preset bypasses ideal packing geometries: they do not contain
interfacial force-field assignments. Results describe this short trajectory.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

from ..engines.lammps import LammpsRunSpec, run_lammps
from .iff_sources import DOI, HASHES, OPLS_SHA, load, sections, fields


def _integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an integer between {low} and {high}")
    return value


def recipe(profile, seed, steps, *, continuation=False):
    charmm = profile == "iff-charmm"
    polymer = profile != "iff-water"
    dihedral = (
        "hybrid harmonic fourier"
        if charmm
        else ("hybrid harmonic opls" if polymer else "harmonic")
    )
    improper = "hybrid cvff harmonic" if charmm else "cvff"
    mixing = "arithmetic" if charmm else "geometric"
    scale = 0.5 if profile == "iff-opls" else 1
    interval = max(1, (steps + 9) // 10)
    lines = [
        "units real",
        "atom_style full",
        "boundary p p p",
        "bond_style harmonic",
        "angle_style harmonic",
        "dihedral_style " + dihedral,
        "improper_style " + improper,
        "pair_style lj/cut/coul/long 12.0",
        "pair_modify mix " + mixing,
        "kspace_style pppm 1e-5",
        f"special_bonds lj/coul 0 0 {scale}",
        "read_data interface.data",
        "kspace_style pppm 1e-5",
    ]
    if continuation:
        # Hybrid dihedral/improper and mixed-pair coefficients are explicitly
        # re-applied from the same hash-checked scientific data after restart.
        lines = [
            "units real",
            "read_restart parent.restart",
            "kspace_style pppm 1e-5",
            "pair_modify mix " + mixing,
            "include restore-coeffs.in",
        ]
    if charmm:
        lines.append("include fixed_iff_pairs.in")
    lines += [
        "neighbor 2.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        "group solid type 1:9",
        "group mobile subtract all solid",
        "fix slab solid setforce 0 0 0",
        "velocity solid set 0 0 0",
        "compute Tmob mobile temp",
        f"thermo {interval}",
        "thermo_style custom step atoms c_Tmob pe ke etotal",
        "thermo_modify format float %.12g flush yes",
        "min_style cg",
        "min_modify dmax 0.1",
        "minimize 0.0 1e-4 25 100",
        "reset_timestep 0",
        f"velocity mobile create 298.15 {seed} mom yes rot yes dist gaussian",
        "fix propagate mobile nvt temp 298.15 298.15 100.0",
        "fix_modify propagate temp Tmob",
        "timestep 0.5",
        f"dump trace all custom {interval} trajectory.dump id type xs ys zs",
        "dump_modify trace first yes sort id format float %.9g",
        f"run {steps}",
        "write_restart state.restart",
        'print "POLYCHEMMD_RUN_COMPLETE"',
    ]
    if continuation:
        omitted = (
            "group solid ",
            "group mobile ",
            "velocity ",
            "min_style ",
            "min_modify ",
            "minimize ",
            "reset_timestep ",
        )
        lines = [line for line in lines if not line.startswith(omitted)]
    return "\n".join(lines) + "\n"


def restart_coefficients(data):
    tables = sections(data.decode())
    records = []
    for table, command in (
        ("Pair Coeffs", "pair_coeff"),
        ("Bond Coeffs", "bond_coeff"),
        ("Angle Coeffs", "angle_coeff"),
        ("Dihedral Coeffs", "dihedral_coeff"),
        ("Improper Coeffs", "improper_coeff"),
    ):
        for row in tables[table]:
            tokens = fields(row)
            if table == "Pair Coeffs":
                tokens.insert(0, tokens[0])
            records.append(" ".join([command, *tokens]))
    return ("\n".join(records) + "\n").encode()


def _checkpoint(source_dir):
    root = Path(source_dir).resolve(strict=True)

    def trusted(name):
        path = root / name
        if path.is_symlink() or not path.resolve(strict=True).is_relative_to(root):
            raise ValueError("Checkpoint must be a source-job artifact")
        return path

    summary = json.loads(trusted("md_summary.json").read_text())
    provenance = json.loads(trusted("engine/run_record.json").read_text())
    checkpoint = trusted("engine/state.restart").read_bytes()
    if (
        provenance["status"] != "completed"
        or not checkpoint
        or hashlib.sha256(checkpoint).hexdigest()
        != provenance["artifact_sha256"]["state.restart"]
        or hashlib.sha256(trusted("engine/input.lammps").read_bytes()).hexdigest()
        != summary["input_script_sha256"]
    ):
        raise ValueError("Parent checkpoint/provenance mismatch")
    return summary, provenance["engine_identity"], checkpoint


def _dump_json(path, data):
    path.write_text(
        json.dumps(data, ensure_ascii=False, allow_nan=False, separators=(",", ":")),
        encoding="utf-8",
    )


def _measure(directory, profile, out):
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
    import numpy as np

    universe = mda.Universe(
        str(directory / "interface.data"),
        str(directory / "trajectory.dump"),
        topology_format="DATA",
        format="LAMMPSDUMP",
        atom_style="id resid type charge x y z",
        lammps_coordinate_convention="scaled",
        dt=0.0005,
    )
    atoms = universe.atoms
    types = atoms.types.astype(int)
    ids = atoms.ids
    central = np.flatnonzero(types == 6)
    z = atoms.positions[:, 2].copy()
    # The source slab sits inside the z cell. Water wraps in z; the fixed stack does not.
    if len(central) != 648:
        raise RuntimeError("Unknown central titanium sublattice")
    by_height = central[np.argsort(z[central])]
    bottom, top = np.split(by_height, 2)
    # In the deposited stack, the two 3x3 supercell repeats are ~9.6 Å apart;
    # only the two outside faces are hydrated. Sort by z around central Ti planes.
    mid = (z[bottom].mean() + z[top].mean()) / 2
    dz = z - mid
    period = float(universe.trajectory.ts.triclinic_dimensions[2, 2])
    dz -= np.rint(dz / period) * period
    exposed = np.isin(types, [7, 8]) & (np.abs(dz) > 6.0)
    surface_f = np.flatnonzero(exposed & (types == 7))
    surface_oh = np.flatnonzero(exposed & (types == 8))
    if (len(surface_f), len(surface_oh)) != (342, 306):
        raise RuntimeError("Reference external-face IDs do not match the pinned source")
    water = np.flatnonzero(types == 10)
    heavy = np.flatnonzero((types >= 12) & (atoms.masses > 2))
    groups = {
        "reference_plane": ids[bottom].tolist(),
        "comparison_plane": ids[top].tolist(),
        "surface_sites": ids[exposed].tolist(),
        "water_oxygen": ids[water].tolist(),
        "polymer_heavy": ids[heavy].tolist(),
        "plane_restraint": "fixed",
        "selection_note": "Central Ti(type 6) in each fixed sheet, external F/OH(O) sites on both water-facing surfaces; all water O(type 10).",
    }
    _dump_json(out / "trajectory-groups.json", groups)
    rows = []
    for ts in universe.trajectory:
        row = {"step": int(ts.data["step"])}
        for label, sites in (("F", surface_f), ("OH", surface_oh)):
            pairs = capped_distance(
                atoms.positions[sites],
                atoms.positions[water],
                3.5,
                box=ts.dimensions,
                method="nsgrid",
                return_distances=False,
            )
            row[label] = len(pairs) / len(sites)
        rows.append(row)
    if not rows:
        raise RuntimeError("No hydration frames")
    with (out / "hydration.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["step", "F", "OH"])
        writer.writeheader()
        writer.writerows(rows)
    # Full last sampled frame in Cartesian Å; IDs and graph edges are retained.
    elements = {1: "H", 12: "C", 16: "O", 19: "F", 48: "Ti"}
    display = [
        {
            "id": int(a.id),
            "element": elements[round(float(a.mass))],
            "role": ("surface" if t <= 9 else "water" if t <= 11 else "polymer"),
            "x": round(float(pos[0]), 5),
            "y": round(float(pos[1]), 5),
            "z": round(float(pos[2]), 5),
        }
        for a, t, pos in zip(atoms, types, atoms.positions, strict=True)
    ]
    wrapped = atoms.positions
    all_bonds = universe.bonds.indices
    visible = (
        np.linalg.norm(wrapped[all_bonds[:, 0]] - wrapped[all_bonds[:, 1]], axis=1)
        < 4.0
    )
    _dump_json(
        out / "structure.json",
        {
            "label": "last sampled MD frame",
            "step": rows[-1]["step"],
            "units": "angstrom",
            "atoms": display,
            "bonds": all_bonds[visible].tolist(),
            "cell_vectors_A": universe.trajectory.ts.triclinic_dimensions.tolist(),
            "bond_display": "Cross-boundary lines omitted from viewer only; topology is complete in interface.data.",
        },
    )
    return {
        "sampled_frames": len(rows),
        "last_sampled_step": rows[-1]["step"],
        "MDAnalysis_version": mda.__version__,
        "water_contacts": {
            label: round(sum(r[label] for r in rows) / len(rows), 7)
            for label in ("F", "OH")
        },
        "contact_definition": "Per-site number of water O within 3.5 Å, triclinic minimum-image distances, mean over sampled frames including step zero. Geometric proximity, no H-bond angle or convergence claim.",
        "atom_count": len(atoms),
        "water_count": len(water),
        "exposed_F_sites": 342,
        "exposed_OH_sites": 306,
    }


def run_md(config: dict, out: Path, *, source_dir=None) -> dict:
    allowed = (
        {"md_steps", "restart_job"}
        if source_dir is not None
        else {"md_steps", "md_profile", "md_sequence", "seed"}
    )
    if set(config) - allowed:
        raise ValueError("Use only fields supported by this MD preset or continuation")
    if source_dir is not None:
        if not config.get("restart_job"):
            raise ValueError("Missing parent job identifier")
        previous, identity, checkpoint = _checkpoint(source_dir)
    else:
        previous, identity, checkpoint = {}, None, None
    profile = previous.get("profile", config.get("md_profile", "iff-water"))
    sequence = previous.get("sequence") or config.get("md_sequence", "clustered")
    if profile not in {"iff-water", "iff-opls", "iff-charmm"} or sequence not in {
        "clustered",
        "dispersed",
    }:
        raise ValueError("Unknown calibrated MD preset or PVA/VAc sequence")
    steps = _integer(config.get("md_steps", 20), "md_steps", 1, 2000)
    seed = _integer(
        previous.get("seed", config.get("seed", 15991)), "seed", 1, 900_000_000
    )
    start = previous.get("end_step", 0)
    # Fail before consuming engine time if the observable analysis dependency is absent.
    try:
        import MDAnalysis  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "Install the studio MD extra with MDAnalysis before commissioning"
        ) from exc
    source, pairs = load(profile, sequence)
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    directory = out / "engine"
    staged = [("interface.data", source)]
    if pairs:
        staged.append(("fixed_iff_pairs.in", pairs))
    if checkpoint is not None:
        if previous.get("source_data_sha256") != hashlib.sha256(source).hexdigest():
            raise ValueError("Parent science input differs from the pinned profile")
        staged.extend(
            [
                ("parent.restart", checkpoint),
                ("restore-coeffs.in", restart_coefficients(source)),
            ]
        )
    spec = LammpsRunSpec(
        input_script=recipe(profile, seed, steps, continuation=checkpoint is not None),
        run_directory=directory,
        expected_engine_identity=identity,
        units_style="real",
        timeout_seconds=1200,
        expected_artifacts=("trajectory.dump", "state.restart"),
        required_thermo_columns=("Step", "c_Tmob", "PotEng"),
        input_binary_files=tuple(staged),
    )
    report = run_lammps(spec).to_dict()
    dynamics = report["thermo_blocks"][-1]
    columns = dynamics["columns"]
    final = dict(zip(columns, dynamics["rows"][-1], strict=True))
    if final["Step"] != start + steps:
        raise RuntimeError("LAMMPS did not finish the requested MD time")
    key = profile if profile == "iff-water" else f"{profile}-{sequence}"
    summary = {
        "profile": profile,
        "sequence": None if profile == "iff-water" else sequence,
        "engine": "LAMMPS",
        "engine_version": report["engine_version_banner"],
        "seed": seed,
        "steps": steps,
        "timestep_fs": 0.5,
        "duration_ps": steps * 0.0005,
        "start_step": start,
        "end_step": start + steps,
        "restart_job": config.get("restart_job"),
        "parent_checkpoint_sha256": hashlib.sha256(checkpoint).hexdigest()
        if checkpoint
        else None,
        "target_temperature_K": 298.15,
        "final_mobile_temperature_K": final["c_Tmob"],
        "source_doi": DOI,
        "source_data_sha256": HASHES[key],
        "scope": "Fixed deposited two-sheet Ti3C2(OH,F) stack, water on outer faces, no added salt. Finite NVT engine commissioning. Distances of frozen sheets are prescribed. No equilibrium, binding free energy, dispersion or 24-hour optical prediction.",
        "input_script_sha256": hashlib.sha256(spec.input_script.encode()).hexdigest(),
        "charge_model": "Native IFF liquid-interface parameters; flexible water.",
        "polymer_source": {
            "iff-water": None,
            "iff-opls": "OPLS/2020, DOI 10.1021/acs.jpcb.3c06602, Moltemplate 2.23.0 source SHA256 "
            + OPLS_SHA,
            "iff-charmm": "Atlas MATCH 0.01 CGenFF v2a molecular assignment",
        }[profile],
        "artifacts": [
            "hydration.csv",
            "structure.json",
            "trajectory-groups.json",
            "engine/trajectory.dump",
            "engine/input.lammps",
            "engine/interface.data",
            "engine/log.lammps",
            "engine/run_record.json",
            "engine/state.restart",
        ],
    }
    summary.update(_measure(directory, profile, out))
    _dump_json(out / "md_summary.json", summary)
    return summary
