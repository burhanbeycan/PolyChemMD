"""Run the three canonical MD profiles with matched split/continuous controls."""

import argparse
import json
from pathlib import Path

import numpy as np

from polychemmd.studio.md import run_md
from polychemmd.studio.trajectory import frames

PROFILES = (
    ("iff-water", "clustered"),
    ("iff-opls", "dispersed"),
    ("iff-charmm", "clustered"),
)
SCOPE = (
    "Finite commissioning on Linux CPU, LAMMPS real units. The comparison uses printed fractional "
    "coordinates (9 significant digits) with triclinic minimum images and thermo printed to "
    "12 significant digits. OPLS and MATCH have different sequences in these commissioning "
    "fixtures. Each restart is compared to its own same-profile unsplit trajectory."
)


def summarize(folder):
    def get(stage, name):
        return json.loads((folder / stage / name).read_text())

    def last_frame(stage):
        last = None
        for frame in frames(folder / stage / "engine/trajectory.dump"):
            last = frame
        return last

    a, b = last_frame("second"), last_frame("direct")
    if (
        not np.array_equal(a.ids, b.ids)
        or not np.allclose(a.cell, b.cell)
        or a.step != b.step
    ):
        raise ValueError("Cannot pair simulation end points")
    delta = a.fractional - b.fractional
    delta -= np.rint(delta)
    maximum = float(np.linalg.norm(delta @ a.cell, axis=1).max())
    direct = get("direct", "md_summary.json")
    second = get("second", "md_summary.json")
    first = get("first", "md_summary.json")

    def potential(stage):
        block = get(stage, "engine/run_record.json")["thermo_blocks"][-1]
        return block["rows"][-1][block["columns"].index("PotEng")]

    energy = abs(potential("second") - potential("direct"))
    temperature = abs(
        second["final_mobile_temperature_K"] - direct["final_mobile_temperature_K"]
    )
    if maximum > 2e-6 or energy > 1e-5 or temperature > 1e-5:
        raise RuntimeError(
            f"Continuation tolerance exceeded: {maximum} Å, {energy} kcal/mol, {temperature} K"
        )
    return {
        "profile": direct["profile"],
        "sequence": direct["sequence"],
        "seed": direct["seed"],
        "engine": direct["engine_version"],
        "analysis_engine": "MDAnalysis " + direct["MDAnalysis_version"],
        "atom_count": direct["atom_count"],
        "water_count": direct["water_count"],
        "source_sha256": direct["source_data_sha256"],
        "input_script_sha256": direct["input_script_sha256"],
        "end_step": a.step,
        "continuous_steps": direct["steps"],
        "split_steps": [first["steps"], second["steps"]],
        "timestep_fs": direct["timestep_fs"],
        "duration_ps": direct["duration_ps"],
        "observation_frames_continuous": direct["sampled_frames"],
        "continuous_final_mobile_temperature_K": direct["final_mobile_temperature_K"],
        "mean_water_O_within_3p5A_per_external_surface_site": direct["water_contacts"],
        "maximum_final_coordinate_difference_A_printed": maximum,
        "absolute_final_temperature_difference_K": temperature,
        "absolute_final_box_potential_energy_difference_kcal_per_mol": energy,
        "source_doi": direct["source_doi"],
        "trajectory_sha256": get("direct", "engine/run_record.json")["artifact_sha256"][
            "trajectory.dump"
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    args.output_dir.mkdir(exist_ok=False, parents=True)
    results = []
    for profile, sequence in PROFILES:
        folder = args.output_dir / profile
        config = {"md_profile": profile, "md_sequence": sequence, "seed": 15991}
        run_md({**config, "md_steps": 10}, folder / "first")
        run_md(
            {"restart_job": "first", "md_steps": 10},
            folder / "second",
            source_dir=folder / "first",
        )
        run_md({**config, "md_steps": 20}, folder / "direct")
        results.append(summarize(folder))
    output = args.output_dir / "iff_engine_checks.json"
    output.write_text(json.dumps({"scope": SCOPE, "runs": results}, indent=2) + "\n")
    print(output)


if __name__ == "__main__":
    main()
