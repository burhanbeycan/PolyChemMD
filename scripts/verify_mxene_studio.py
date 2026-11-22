"""Replay the registered Studio tasks through its queue. No web service is needed.

Run in a new output directory after installing the engines and model data.
Only a small machine-readable measurement ledger is intended for version control.
"""

import argparse
import json
from pathlib import Path

from polychemmd.studio.jobs import JobStore, execute_next
from polychemmd.studio.operations import versions


def verify(directory, with_md):
    if directory.exists() and any(directory.iterdir()):
        raise ValueError("Choose a new or empty verification directory.")
    store = JobStore(directory)
    reports = []

    def run(kind, config, fields):
        job = store.submit(kind, config)
        execute_next(store, job["id"])
        actual = store.get(job["id"])
        if actual["status"] != "succeeded":
            raise RuntimeError(f"{kind} failed: {actual['error']}")
        result = actual["result"]
        reports.append(
            {
                "kind": kind,
                "config": {k: v for k, v in config.items() if k != "source_job"},
                "output": {k: result[k] for k in fields if k in result},
            }
        )
        return actual

    for backend in ("builtin", "molpack", "packmol"):
        j = run(
            "packing",
            {"polymer": "PEO", "dp": 8, "packing_engine": backend, "seed": 17},
            (
                "engine",
                "atom_count",
                "water_count",
                "ion_pairs",
                "salt_M_requested",
                "salt_M_realized",
                "audit",
            ),
        )
        if j["result"]["ion_pairs"] != 1:
            raise AssertionError("Unexpected ion pair count")
    run(
        "radonpy",
        {"polymer": "PEO", "dp": 8, "seed": 17},
        (
            "engine",
            "n_atoms",
            "heavy_atoms",
            "formula",
            "software",
            "gaff2_typed",
            "charge_model",
        ),
    )
    run(
        "experiment",
        {"assay": "xrd", "reference": "xing-pei"},
        ("evidence", "summary", "source_doi"),
    )
    if with_md:
        for profile in ("iff-water", "iff-opls", "iff-charmm"):
            j = run(
                "md",
                {
                    "md_profile": profile,
                    "md_sequence": "clustered",
                    "md_steps": 20,
                    "seed": 15991,
                },
                (
                    "profile",
                    "last_sampled_step",
                    "sampled_frames",
                    "atom_count",
                    "water_contacts",
                    "MDAnalysis_version",
                ),
            )
            if j["result"]["last_sampled_step"] != 20:
                raise AssertionError("MD did not reach the requested step")
            tr = run(
                "trajectory",
                {"source_job": j["id"]},
                (
                    "n_frames",
                    "mean_water_O_contact_fraction",
                    "mean_plane_separation_nm",
                    "plane_restraint",
                    "bin_width_A",
                    "contact_cutoff_A",
                ),
            )
            if tr["result"]["n_frames"] != 11:
                raise AssertionError("Expected eleven saved frames")
    ledger = {
        "software": versions(),
        "scope": "Short engineering and matched-readout checks. No equilibrium or experimental prediction inferred.",
        "jobs": reports,
    }
    (directory / "verification.json").write_text(
        json.dumps(ledger, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(directory / "verification.json")


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--with-md", action="store_true")
    args = p.parse_args()
    verify(args.output_dir, args.with_md)
