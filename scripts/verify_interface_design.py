"""Replay the controlled PVA interface end to end and audit native engine outputs.

Use a new output directory with the Studio dependencies and LAMMPS installed.
The adjacent report records an installation and force-consistency check at one
constructed geometry. No thermodynamic or experimental property is estimated.
"""

import argparse
import csv
import json
import math
import os
from hashlib import sha256
from itertools import pairwise
from pathlib import Path

import numpy as np

from polychemmd.engines.lammps import LammpsRunSpec, run_lammps
from polychemmd.studio.jobs import JobStore, execute_next
from polychemmd.studio.operations import versions

RECIPE = {
    "polymer": "PVA",
    "dp": 8,
    "hydrolysis": 0.75,
    "acetate_pattern": "dyads",
    "acetate_dyads": 1,
    "pattern": "correlated",
    "oh_top": 0.5,
    "oh_bottom": 0.5,
    "surface_correlation_top": 0.22,
    "surface_correlation_bottom": -0.10,
    "geometry": "single",
    "nx": 7,
    "ny": 5,
    "water_height_A": 22.0,
    "chains": 1,
    "water_count": 32,
    "salt_M": 0.0,
    "min_distance_A": 2.0,
    "packing_engine": "slit",
    "seed": 17,
    "tacticity_seed": 101,
    "sequence_seed": 109,
    "surface_seed": 113,
    "conformer_seed": 127,
}


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def queued(store, kind, recipe):
    job = store.submit(kind, recipe)
    execute_next(store, job["id"])
    job = store.get(job["id"])
    if job["status"] != "succeeded":
        raise RuntimeError(f"{kind}: {job['error']}")
    return job, store.output_dir(job["id"])


def count_acetates(packed):
    # Count acetate carbonyl C from its actual covalent C=O bond, independent
    # of the recipe, dyad sampler and the published polymer summary.
    atoms = packed["atoms"]
    positions = []
    for a, b, order in packed["bonds"]:
        if order != 2:
            continue
        pair = [atoms[a], atoms[b]]
        if all(v["role"] == "polymer" for v in pair):
            carbon = next(v for v in pair if v["element"] == "C")
            oxygen = next(v for v in pair if v["element"] == "O")
            require(
                carbon["repeat_index"] == oxygen["repeat_index"],
                "Cross-repeat carbonyl",
            )
            positions.append(carbon["repeat_index"] + 1)
    positions.sort()
    dyads = sum(b == a + 1 for a, b in pairwise(positions))
    require(
        len(positions) == 2 and dyads == RECIPE["acetate_dyads"],
        "PVA chemistry mismatch",
    )
    return {"carbonyl_positions_1based": positions, "adjacent_acetate_dyads": dyads}


def surface_check(packed):
    dims = np.array(packed["box_A"][:2])
    results = {}
    for face in ("top", "bottom"):
        sites = [
            v
            for v in packed["atoms"]
            if v.get("layer") in {face + "_OH", face + "_F"}
            and v["element"] in {"O", "F"}
        ]
        xy = np.array([[v["x"], v["y"]] for v in sites])
        diff = xy[:, None, :] - xy[None, :, :]
        diff -= np.rint(diff / dims) * dims
        distances = np.linalg.norm(diff, axis=2)
        nearest = np.min(distances[distances > 1.0e-6])
        edges = np.isclose(distances, nearest, atol=1e-6, rtol=0.0)
        require(np.all(edges.sum(axis=1) == 6), "Expected a periodic triangular net")
        labels = np.array([v["element"] == "O" for v in sites], dtype=float)
        fraction = labels.mean()
        target = RECIPE["surface_correlation_" + face]
        rho = float(
            ((labels - fraction)[:, None] * (labels - fraction)[None, :])[edges].mean()
            / (fraction * (1.0 - fraction))
        )
        require(
            abs(fraction - RECIPE["oh_" + face]) < 1.0e-12, "Incorrect face composition"
        )
        require(
            abs(rho - target) < 0.05, "Target correlation missed for bounded fixture"
        )
        results[face] = {
            "oh_sites": int(labels.sum()),
            "total_sites": len(sites),
            "undirected_neighbor_bonds": int(edges.sum() // 2),
            "requested_correlation": target,
            "coordinate_based_correlation": rho,
        }
    return results


def frames_with_forces(path):
    with path.open(encoding="ascii") as inp:
        while inp.readline().strip() == "ITEM: TIMESTEP":
            step = int(inp.readline())
            require(
                inp.readline().startswith("ITEM: NUMBER OF ATOMS"),
                "Malformed atom count",
            )
            count = int(inp.readline())
            require(
                inp.readline().startswith("ITEM: BOX BOUNDS pp pp ff"),
                "Unexpected slab boundary",
            )
            bounds = np.array([inp.readline().split() for _ in range(3)], float)
            columns = inp.readline().split()[2:]
            values = np.array([inp.readline().split() for _ in range(count)], float)
            require(np.isfinite(values).all(), "Nonfinite native dump")
            values = values[values[:, columns.index("id")].argsort()]
            require(
                np.array_equal(values[:, columns.index("id")], np.arange(1, count + 1)),
                "Lost atom IDs",
            )
            xyz = values[:, [columns.index(k) for k in ("xs", "ys", "zs")]]
            force = values[:, [columns.index(k) for k in ("fx", "fy", "fz")]]
            yield step, xyz * (bounds[:, 1] - bounds[:, 0]) + bounds[:, 0], force


def finite_difference(source, atom, axis, shift, destination):
    # Same native coefficients, cell, exclusions and PPPM settings. Change only
    # one Cartesian coordinate before run zero; no MD trajectory enters dE/dx.
    original = (source / "registered.in").read_text(encoding="ascii")
    upstream, command, tail = original.partition("run 0\n")
    require(bool(command) and "run 0" not in tail, "Expected a static registered run")
    direction = [shift if i == axis else 0.0 for i in range(3)]
    displacement = " ".join(f"{value:.12g}" for value in direction)
    script = (
        upstream
        + f"group testatom id {atom}\n"
        + f"displace_atoms testatom move {displacement} units box\n"
        + command
        + tail
    )
    data = (source / "engine" / "typed.data").read_text(encoding="ascii")
    calculation = run_lammps(
        LammpsRunSpec(
            input_script=script,
            run_directory=destination,
            units_style="real",
            required_thermo_columns=("Step", "PotEng"),
            input_files=(("typed.data", data),),
            timeout_seconds=120,
        )
    )
    thermo = calculation.to_dict()["thermo_blocks"][-1]
    return thermo["rows"][-1][thermo["columns"].index("PotEng")]


def audit_motion(initial, final, packed):
    box = np.array(packed["box_A"])
    delta = final - initial
    delta[:, :2] -= np.rint(delta[:, :2] / box[:2]) * box[:2]
    distance = np.linalg.norm(delta, axis=1)
    solid = np.array([a["role"] == "slab" for a in packed["atoms"]])
    worst_fixed = float(distance[solid].max())
    largest_mobile = float(distance[~solid].max())
    require(
        worst_fixed < 1e-8 and largest_mobile > 1e-7,
        "Frozen group or mobile integration failed",
    )
    return {
        "maximum_fixed_displacement_A": worst_fixed,
        "maximum_mobile_displacement_A": largest_mobile,
    }


def verify(output):
    if output.exists() and any(output.iterdir()):
        raise ValueError("Use an empty output directory")
    output.mkdir(exist_ok=True, parents=True)
    # Keep this verification's generated registry and SQLite store self-contained.
    os.environ["POLYCHEMMD_TOPOLOGY_REGISTRY"] = str((output / "registry").resolve())
    store = JobStore(output / "queue")
    pack, pack_dir = queued(store, "packing", RECIPE)
    packed = read(pack_dir / "packed.json")
    audit = {
        "recipe": RECIPE,
        "polymer": count_acetates(packed),
        "faces": surface_check(packed),
    }
    geometry = read(pack_dir / "structure.json")
    for face, value in audit["faces"].items():
        recorded = geometry["face_occupancy"][0][face]["neighbor_statistics"]
        require(
            abs(recorded["oh_correlation"] - value["coordinate_based_correlation"])
            < 1e-12,
            "Recorded surface correlation differs from Cartesian pair enumeration",
        )
    require(
        audit["polymer"]["carbonyl_positions_1based"]
        == geometry["polymer"]["acetate_residues_1based"],
        "Reported sequence positions differ from the packed bond graph",
    )
    topo, topo_dir = queued(
        store, "topology", {"source_job": pack["id"], "cutoff_A": 8.0}
    )
    key = topo["result"]["topology_id"]
    pilot = {
        "profile": "registered",
        "source_job": pack["id"],
        "topology_id": key,
        "seed": 7001,
    }
    static, static_dir = queued(store, "md", {**pilot, "steps": 0})
    dynamic, dynamic_dir = queued(store, "md", {**pilot, "steps": 4})
    _, analysis_dir = queued(store, "trajectory", {"source_job": dynamic["id"]})
    (static_frame,) = list(
        frames_with_forces(static_dir / "registered-trajectory.dump")
    )
    frames = list(frames_with_forces(dynamic_dir / "registered-trajectory.dump"))
    require([f[0] for f in frames] == list(range(5)), "Unexpected sampled timesteps")
    raw_atoms = packed["atoms"]
    n = len(raw_atoms)
    with (topo_dir / "atom-types.csv").open(encoding="utf-8") as infile:
        assignments = list(csv.DictReader(infile))
    require(len(assignments) == n, "Incomplete per-atom parameter trace")
    by_molecule = {}
    for i, (at, row) in enumerate(zip(raw_atoms, assignments, strict=True), 1):
        require(
            int(row["atom_id_1based"]) == i
            and row["element"] == at["element"]
            and row["molecule"] == at["molecule_id"],
            "Changed atom in source parameter map",
        )
        by_molecule.setdefault(row["molecule"], []).append(float(row["charge_e"]))
    largest_molecule_charge = max(abs(math.fsum(q)) for q in by_molecule.values())
    require(largest_molecule_charge < 1e-8, "Component charge imbalance")
    require(
        topo["result"]["packed_sha256"] == digest(pack_dir / "packed.json")
        and static["result"]["parameter_data_sha256"]
        == digest(topo_dir / "typed.data"),
        "Lost source-to-engine file identity",
    )
    raw = np.array([[at[k] for k in "xyz"] for at in raw_atoms])
    require(static_frame[1].shape == raw.shape, "Native atom count changed")
    diff = static_frame[1] - raw
    dims = np.array(packed["box_A"])
    diff[:, :2] -= np.rint(diff[:, :2] / dims[:2]) * dims[:2]
    position_error = float(np.linalg.norm(diff, axis=1).max())
    require(position_error < 2e-7, "Packed coordinates differ from engine step zero")
    require(
        np.array_equal(frames[0][1], static_frame[1]),
        "Static and dynamic launches differ",
    )
    force = static_frame[2]
    mobile = np.array([at["role"] != "slab" for at in raw_atoms])
    fmax = float(np.linalg.norm(force[mobile], axis=1).max())
    require(
        np.isclose(fmax, static["result"]["initial_mobile_fmax_kcal_mol_A"], rtol=1e-8),
        "Computed mobile maximum force differs from LAMMPS thermo",
    )
    polymer_oxygens = np.array(
        [
            i
            for i, at in enumerate(raw_atoms)
            if at["role"] == "polymer" and at["element"] == "O"
        ]
    )
    at, axis = np.unravel_index(
        np.abs(force[polymer_oxygens]).argmax(), (len(polymer_oxygens), 3)
    )
    index = int(polymer_oxygens[at])
    h = 0.001
    energies = [
        finite_difference(static_dir, index + 1, axis, sign * h, output / label)
        for sign, label in ((-1, "minus"), (1, "plus"))
    ]
    predicted = (energies[0] - energies[1]) / (2 * h)
    actual = float(force[index, axis])
    require(
        math.isclose(predicted, actual, rel_tol=0.01, abs_tol=0.06),
        f"Force consistency failed: engine {actual}, finite difference {predicted}",
    )
    with (analysis_dir / "density.csv").open(encoding="utf-8") as infile:
        water_count = sum(
            float(row["mean_bin_volume_A3"]) * float(row["water_O_number_density_A-3"])
            for row in csv.DictReader(infile)
        )
    require(
        abs(water_count - RECIPE["water_count"]) < 1e-6,
        "Water density normalization changed",
    )
    result = {
        "scope": "Finite construction, topology, energy/force consistency and integration at one sparse aqueous test composition.",
        "software": versions(),
        **audit,
        "packing": {
            "atoms": n,
            "water_molecules": pack["result"]["water_count"],
            "salt_pairs": pack["result"]["ion_pairs"],
            "closest_mobile_pair_or_slab_contact_A": pack["result"]["audit"][
                "minimum_mobile_to_other_A"
            ],
            "coordinates_sha256": digest(pack_dir / "packed.json"),
        },
        "topology": {
            "file_sha256": digest(topo_dir / "typed.data"),
            "charge_e": static["result"]["charge_e"],
            "max_abs_component_charge_e": largest_molecule_charge,
            "source_IFF_sha256": topo["result"]["parameter_provenance"][
                "native_sha256"
            ],
            "source_OPLS_table_sha256": topo["result"]["parameter_provenance"][
                "opls_sha256"
            ],
            "cutoff_A": 8.0,
            "packed_to_engine_max_coordinate_error_A": position_error,
        },
        "native_energy_force": {
            "static_potential_kcal_mol": static["result"]["final_kcal_mol"],
            "polymer_oxygen_atom_id": index + 1,
            "cartesian_component": "xyz"[axis],
            "initial_mobile_max_force_kcal_mol_A": fmax,
            "delta_A": h,
            "engine_force_kcal_mol_A": actual,
            "central_difference_kcal_mol_A": predicted,
            "abs_difference_kcal_mol_A": abs(actual - predicted),
            "acceptance_relative_tolerance": 0.01,
            "acceptance_absolute_tolerance_kcal_mol_A": 0.06,
        },
        "trajectory": {
            "nve_steps": 4,
            "timestep_fs": dynamic["result"]["dt_fs"],
            "frames": len(frames),
            "water_count_recovered_from_density": water_count,
            **audit_motion(frames[0][1], frames[-1][1], packed),
        },
    }
    out = output / "verification.json"
    out.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return out


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--output-dir", type=Path, required=True)
    print(verify(cli.parse_args().output_dir))
