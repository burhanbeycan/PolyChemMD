"""Paired packing-engine benchmark on fixed molecular templates and measured wall time."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
from importlib import metadata
import json
import os
from pathlib import Path
import platform
import shutil
import signal
from statistics import median
import subprocess
import sys
import tempfile
from time import perf_counter

from .packing import ENGINE_BORDER_A, ENGINE_DISTANCE_BUFFER_A, pack_system


CASES = {
    "aqueous": dict(
        geometry="single",
        polymer="PEO",
        dp=8,
        nx=6,
        ny=4,
        water_height_A=22,
        vacuum_padding_A=12,
        chains=1,
        salt_M=0,
        water_density=0.997,
        min_distance_A=1.8,
        seed=17,
    ),
    "confined": dict(
        geometry="slit",
        polymer="PVA",
        dp=12,
        hydrolysis=0.75,
        acetate_pattern="clustered",
        nx=6,
        ny=4,
        gap_A=9,
        chains=1,
        salt_M=0,
        water_density=0.85,
        min_distance_A=1.8,
        seed=17,
    ),
}
ENGINES = ("builtin", "slit", "molpack", "packmol")
TIME_BUDGET_S = 60
FORCE_MIN_STEPS = 60
FORCE_TIME_BUDGET_S = 150


def _list(config, key, default, options=None):
    values = config.get(key, list(default))
    if not isinstance(values, list) or not values:
        raise ValueError(f"{key} must be a nonempty list")
    if options is None:
        if len(values) > 8 or any(
            type(v) is not int or not 1 <= v <= 2147483646 for v in values
        ):
            raise ValueError("Choose 1 to 8 positive 32-bit packing seeds")
    elif any(not isinstance(v, str) or v not in options for v in values):
        raise ValueError(f"Invalid {key}; allowed values are {', '.join(options)}")
    if len(values) != len(set(values)):
        raise ValueError(f"Repeated {key} would duplicate paired observations")
    return values


def _worker(directory):
    request = json.loads((directory / "pack-request.json").read_text())
    started = perf_counter()
    try:
        packed = pack_system(request, directory, matched_constraints=True)
        measurement = {
            "status": "success",
            "reason": "",
            "minimum_mobile_to_other_A": packed["audit"]["minimum_mobile_to_other_A"],
            "rigid_molecule_identity": packed["audit"]["rigid_molecule_identity"],
            "fixed_surface_invariant": packed["audit"]["fixed_surface_invariant"],
            "proper_rigid_transforms": packed["audit"]["proper_rigid_transforms"],
            "self_images_checked": packed["audit"]["mobile_self_image_check"],
            "atom_count": packed["atom_count"],
            "water_count": packed["water_count"],
            "packed_xyz_sha256": packed["packed_xyz_sha256"],
            "placement_attempts": packed["engine_info"].get("total_placement_attempts"),
            "rejected_rigid_starts": len(packed["engine_info"]["rejected_starts"])
            if "rejected_starts" in packed["engine_info"]
            else None,
            "polymer_orientation_proposals": packed["engine_info"].get(
                "polymer_orientation_proposals"
            ),
        }
    except (
        RuntimeError,
        ValueError,
        OSError,
        subprocess.SubprocessError,
        ImportError,
    ) as exc:
        reason = f"{type(exc).__name__}: {exc}"
        unavailable = isinstance(exc, ImportError) or " unavailable." in reason
        measurement = {
            "status": "unavailable" if unavailable else "failed",
            "reason": reason,
        }
    measurement["wall_seconds"] = perf_counter() - started
    (directory / "trial.json").write_text(json.dumps(measurement, allow_nan=False))


def _subprocess(directory, flag, budget):
    """Keep native child processes within the same wall budget on Linux/WSL."""
    env = os.environ.copy()
    env.update(
        {
            k: "1"
            for k in (
                "OMP_NUM_THREADS",
                "MKL_NUM_THREADS",
                "OPENBLAS_NUM_THREADS",
                "RAYON_NUM_THREADS",
            )
        }
    )
    env["POLYCHEMMD_TOPOLOGY_REGISTRY"] = str((directory / "registry").resolve())
    with subprocess.Popen(
        [
            sys.executable,
            "-m",
            "polychemmd.studio.packing_benchmark",
            flag,
            str(directory.resolve()),
        ],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        start_new_session=os.name == "posix",
    ) as proc:
        try:
            _, errors = proc.communicate(timeout=budget)
        except subprocess.TimeoutExpired:
            if os.name == "posix":
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                proc.kill()
            proc.communicate()
            raise
        if proc.returncode:
            raise RuntimeError(f"Worker exit {proc.returncode}: {errors[-500:]}")


def _timed_pack(directory, config):
    (directory / "pack-request.json").write_text(json.dumps(config, allow_nan=False))
    started = perf_counter()
    try:
        _subprocess(directory, "--worker", TIME_BUDGET_S)
        measurement = json.loads((directory / "trial.json").read_text())
    except subprocess.TimeoutExpired:
        measurement = {
            "status": "timeout",
            "reason": f"Exceeded shared {TIME_BUDGET_S} s subprocess wall limit",
        }
    except (OSError, RuntimeError) as exc:
        measurement = {"status": "failed", "reason": str(exc)}
    measurement["process_wall_seconds"] = perf_counter() - started
    return measurement


def _measure_forces(directory):
    from .register_topology import run_registered
    from .topology_build import run_topology

    top = directory / "topology"
    md = directory / "force"
    top.mkdir()
    md.mkdir()
    typed = run_topology({"cutoff_A": 8.0}, top, source_dir=directory)
    result = run_registered(
        {
            "topology_id": typed["topology_id"],
            "profile": "registered",
            "steps": 0,
            "minimize_steps": FORCE_MIN_STEPS,
        },
        md,
        directory,
    )
    return {
        "force_status": "success",
        "force_topology_id": typed["topology_id"],
        "force_topology_sha256": typed["topology_sha256"],
        "force_cutoff_A": typed["cutoff_A"],
        "minimize_iterations": result["minimize_iterations"],
        "initial_potential_kcal_mol": result["initial_kcal_mol"],
        "final_potential_kcal_mol": result["post_minimize_kcal_mol"],
        "initial_mobile_fmax_kcal_mol_A": result["initial_mobile_fmax_kcal_mol_A"],
        "final_mobile_fmax_kcal_mol_A": result["final_mobile_fmax_kcal_mol_A"],
    }


def _summary(rows, cases, engines):
    summary = []
    for case in cases:
        for engine in engines:
            subset = [
                row for row in rows if row["case"] == case and row["engine"] == engine
            ]
            passed = [row for row in subset if row["status"] == "success"]
            counts = Counter(row["status"] for row in subset)
            summary.append(
                {
                    "case": case,
                    "engine": engine,
                    "n_requested": len(subset),
                    "success_count": counts["success"],
                    "failure_count": counts["failed"],
                    "unavailable_count": counts["unavailable"],
                    "timeout_count": counts["timeout"],
                    "median_success_s": median(row["wall_seconds"] for row in passed)
                    if passed
                    else None,
                    "minimum_distance_A": min(
                        (row["minimum_mobile_to_other_A"] for row in passed),
                        default=None,
                    ),
                    "force_success_count": sum(
                        row["force_status"] == "success" for row in subset
                    ),
                }
            )
    return summary


def _paired_comparisons(rows, cases, engines):
    """Keep unsuccessful paired starts visible and compare timings only for joint successes."""
    reference = "builtin" if "builtin" in engines else engines[0]
    comparisons = []
    for case in cases:
        groups = {
            engine: {
                r["seed"]: r
                for r in rows
                if r["engine"] == engine and r["case"] == case
            }
            for engine in engines
        }
        for engine in engines:
            if engine == reference:
                continue
            paired = [
                (groups[reference][seed], row) for seed, row in groups[engine].items()
            ]
            both = [
                (a, b) for a, b in paired if a["status"] == b["status"] == "success"
            ]
            comparisons.append(
                {
                    "case": case,
                    "reference": reference,
                    "engine": engine,
                    "joint_success_count": len(both),
                    "reference_only_success_count": sum(
                        a["status"] == "success" and b["status"] != "success"
                        for a, b in paired
                    ),
                    "engine_only_success_count": sum(
                        a["status"] != "success" and b["status"] == "success"
                        for a, b in paired
                    ),
                    "neither_success_count": sum(
                        a["status"] != "success" and b["status"] != "success"
                        for a, b in paired
                    ),
                    "median_paired_seconds_engine_minus_reference": median(
                        b["wall_seconds"] - a["wall_seconds"] for a, b in both
                    )
                    if both
                    else None,
                }
            )
    return comparisons


def _versions():
    versions = {"python": platform.python_version(), "platform": platform.platform()}
    for pkg in (
        "numpy",
        "scipy",
        "ase",
        "rdkit",
        "molcrafts-molpack",
        "packmol",
        "lammps",
        "moltemplate",
    ):
        try:
            versions[pkg] = metadata.version(pkg)
        except metadata.PackageNotFoundError:
            versions[pkg] = None
    return versions


def run_benchmark(config, out):
    """All selected engines receive the same prebuilt templates for each seed and case."""
    from .geometry import build_geometry

    allowed = {"cases", "engines", "seeds", "evaluate_forces"}
    if not isinstance(config, dict) or config.keys() - allowed:
        raise ValueError("Unknown packing benchmark option")
    cases = _list(config, "cases", CASES, CASES)
    engines = _list(config, "engines", ENGINES, ENGINES)
    seeds = _list(config, "seeds", (17, 31, 53))
    forces = config.get("evaluate_forces", False)
    if not isinstance(forces, bool):
        raise ValueError("evaluate_forces must be boolean")
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    for case in cases:
        with tempfile.TemporaryDirectory(prefix="packing-bench-", dir=out) as tmp:
            work = Path(tmp)
            geometry_dir = work / "template"
            build_geometry(CASES[case], geometry_dir)
            source = geometry_dir / "structure.json"
            geometry_hash = hashlib.sha256(source.read_bytes()).hexdigest()
            # Rotate execution order across seeds to distribute machine warm-up effects.
            for i, seed in enumerate(seeds):
                first = i % len(engines)
                for engine in engines[first:] + engines[:first]:
                    sub = work / f"{seed}-{engine}"
                    sub.mkdir()
                    shutil.copy2(source, sub / "structure.json")
                    recipe = {
                        **CASES[case],
                        "packing_engine": engine,
                        "packing_seed": seed,
                    }
                    row = {
                        "case": case,
                        "engine": engine,
                        "seed": seed,
                        "geometry_sha256": geometry_hash,
                    }
                    row.update(_timed_pack(sub, recipe))
                    row["force_status"] = (
                        "not_requested" if not forces else "packing_incomplete"
                    )
                    if forces and row["status"] == "success":
                        try:
                            _subprocess(sub, "--force-worker", FORCE_TIME_BUDGET_S)
                            row.update(
                                json.loads((sub / "force-trial.json").read_text())
                            )
                        except (
                            RuntimeError,
                            ValueError,
                            OSError,
                            subprocess.SubprocessError,
                            ImportError,
                        ) as exc:
                            row.update(
                                force_status="failed",
                                force_reason=f"{type(exc).__name__}: {exc}",
                            )
                    rows.append(row)
                    shutil.rmtree(sub)
    summary = _summary(rows, cases, engines)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with (out / "benchmark_trials.csv").open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "stage": "paired_rigid_packing_benchmark",
        "rows": rows,
        "summary": summary,
        "paired_comparisons": _paired_comparisons(rows, cases, engines),
        "cases": {case: CASES[case] for case in cases},
        "seeds": seeds,
        "software": _versions(),
        "pack_wall_limit_s_per_trial": TIME_BUDGET_S,
        "common_generation_distance_buffer_A": ENGINE_DISTANCE_BUFFER_A,
        "common_generation_boundary_margin_A": ENGINE_BORDER_A,
        "thread_caps": 1,
        "evaluate_forces": forces,
        "force_min_iterations_requested": FORCE_MIN_STEPS if forces else None,
        "protocol": (
            "Fixed chain/surface template within each case, matched placement seeds across engines. "
            "Wall seconds include packing and independent audits inside the isolated worker; process time also includes interpreter start. "
            "Force diagnostics, when requested, use the same generated Hamiltonian and iteration limit for a case. "
            "Each case is an initial geometry. Sampling and physical stability require separate equilibrated runs."
        ),
    }
    (out / "packing_benchmark.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--force-worker", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--config", type=Path)
    args = parser.parse_args()
    if args.worker:
        _worker(args.worker)
    elif args.force_worker:
        try:
            values = _measure_forces(args.force_worker)
        except (
            RuntimeError,
            ValueError,
            OSError,
            subprocess.SubprocessError,
            ImportError,
        ) as exc:
            detail = f"{type(exc).__name__}: {exc}"
            log = args.force_worker / "force" / "engine" / "log.lammps"
            if log.is_file():
                errors = [
                    ln.strip()
                    for ln in log.read_text(errors="replace").splitlines()
                    if ln.startswith("ERROR")
                ]
                if errors:
                    detail += f" ({errors[-1][:300]})"
            values = {"force_status": "failed", "force_reason": detail}
        (args.force_worker / "force-trial.json").write_text(
            json.dumps(values, allow_nan=False)
        )
    elif args.output:
        settings = json.loads(args.config.read_text()) if args.config else {}
        print(json.dumps(run_benchmark(settings, args.output)["summary"], indent=2))
    else:
        parser.error("--output is required to run a benchmark.")
