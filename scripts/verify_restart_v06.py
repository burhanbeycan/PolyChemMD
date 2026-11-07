#!/usr/bin/env python3
"""Execute a bounded same-cell native NPT restart reproducibility audit.

Four engine invocations are predeclared: preparation, one 40-step continuation,
and two 20-step continuations from the identical preparation checkpoint. Native
binary checkpoints remain private; hashes and actual text outputs are public.
This is an integration check, not an equilibrated-density calculation.
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import time

from polychemmd import __version__
from polychemmd.analysis.density import density_from_mass_volume
from polychemmd.analysis.sampling import DensitySamplingSpec, analyze_density_sampling
from polychemmd.forcefields.periodic import assign_periodic_alkanes
from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, run_periodic_alkane_protocol


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE_CELL_PATH = "evidence/geometry_v03/cell.json"
TOTAL_WALL_BUDGET_SECONDS = 240.0
PER_STAGE_TIMEOUT_SECONDS = 60.0
FULL_STEPS = 40
SPLIT_STEPS = 20
# Fixed before the engine is invoked; failures do not adjust these margins.
TOLERANCES = {
    "position_component_angstrom": 1e-9,
    "velocity_component_angstrom_fs": 1e-10,
    "box_component_angstrom": 1e-9,
    "energy_kj_mol": 1e-7,
    "force_component_kj_mol_angstrom": 1e-6,
}
PREPARATION = {
    "temperature_k": 300.0, "pressure_bar": 1.0, "timestep_fs": 0.25,
    "minimize_max_iterations": 200, "minimize_max_evaluations": 1000,
    "nvt_steps": 100, "npt_steps": 100, "velocity_seed": 49217,
    "timeout_seconds": PER_STAGE_TIMEOUT_SECONDS, "cutoff_angstrom": 10.0,
    "ewald_alpha_per_angstrom": 0.30, "reciprocal_nmax": 16,
}
SOURCE_MODULES = (
    "polychemmd", "polychemmd.engines.lammps", "polychemmd.engines.periodic_alkane",
    "polychemmd.forcefields.periodic", "polychemmd.protocols.alkane",
    "polychemmd.protocols.restart", "polychemmd.analysis.density", "polychemmd.analysis.sampling",
)
TEXT_ARTIFACTS = (
    "input.lammps", "system.data", "log.lammps", "prepared.data", "final.dump",
    "checkpoint.dump", "resumed_initial.data", "resumed_initial.dump",
)


class RestartAuditError(RuntimeError):
    """A predeclared source, resource or reproducibility check failed."""


def require(condition, message):
    if not condition:
        raise RestartAuditError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical_bytes(record):
    return json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def public_record(raw):
    """Retain scientific output/hash identity while removing host-only paths."""
    record = deepcopy(raw)
    changes = []
    def visit(value, location="record"):
        if isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, f"{location}[{index}]")
        elif isinstance(value, dict):
            if "run_directory" in value:
                value["run_directory"] = "private_engine_run_not_published"
                changes.append(location+".run_directory replaced")
            if isinstance(value.get("executable"), str):
                value["executable"] = Path(value["executable"]).name
                changes.append(location+".executable reduced to basename; identity hashes retained")
            if isinstance(value.get("command"), list) and value["command"]:
                value["command"][0] = Path(value["command"][0]).name
            line = value.get("engine_help_first_line")
            if isinstance(line, str) and any(word in line for word in ("UCX", "hwloc")):
                value["engine_help_first_line"] = "Incidental host diagnostic omitted; see engine_version_banner"
                changes.append(location+".engine_help_first_line omitted host diagnostic")
            for key, item in value.items():
                visit(item, location+"."+key)
    visit(record)
    record["public_record_sanitization"] = {
        "changes": changes,
        "scope": "Only host-only path/diagnostic text changes; scientific values, native checkpoint identities and raw artifact hashes remain explicit.",
    }
    return record


def source_records():
    records = {}
    for name in SOURCE_MODULES:
        actual = Path(importlib.import_module(name).__file__).resolve()
        relative = Path("src")/Path(*name.split("."))
        expected = REPOSITORY/(relative/"__init__.py" if name == "polychemmd" else relative.with_suffix(".py"))
        require(actual.is_file() and expected.is_file(), "Run this auditor from an editable source checkout")
        actual_hash = digest(actual.read_bytes())
        require(actual_hash == digest(expected.read_bytes()), f"Imported {name} differs from checkout")
        records[expected.relative_to(REPOSITORY).as_posix()] = {"sha256": actual_hash}
    records["scripts/verify_restart_v06.py"] = {"sha256": digest(Path(__file__).read_bytes())}
    return records


def hardware_record():
    model = None
    cpuinfo = Path("/proc/cpuinfo")
    if cpuinfo.is_file():
        for line in cpuinfo.read_text(encoding="utf-8").splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {"platform": platform.platform(), "machine": platform.machine(),
            "processor": platform.processor() or model, "logical_cpu_count": os.cpu_count(),
            "engine_parallelism": "one MPI rank, one OpenMP thread; no GPU requested",
            "timing_interpretation": "Measured on this executing host; not a benchmark of the project owner's desktop or GPU"}


def audit_source_cell(source):
    atoms = source["atoms"]
    molecules = Counter(atom["molecule_index"] for atom in atoms)
    require(len(atoms) == 128 and len(molecules) == 4 and set(molecules.values()) == {32},
            "Audit fixture must remain four explicit C10H22 molecules (128 atoms)")
    for index in molecules:
        elements = Counter(atom["element"] for atom in atoms if atom["molecule_index"] == index)
        require(elements == {"C": 10, "H": 22}, "Audit fixture chemistry differs from decane")
    require(source["box"]["lengths_angstrom"] == [24.0, 24.0, 24.0], "Audit initial box changed")
    return {"molecules": 4, "atoms": 128, "formula_per_molecule": "C10H22",
            "initial_box_lengths_angstrom": [24.0]*3,
            "initial_density_interpretation": "Dilute construction fixture; no liquid or polymer-melt phase validation",
            "independently_prepared_cells": 1,
            "restart_branches_interpretation": "Alternative executions of one cell and identical checkpoint; not independent samples"}


def _maximum_difference(first, second, label):
    require(len(first) == len(second) and len(first) > 0, f"{label} dimensions differ")
    values = []
    for left, right in zip(first, second):
        require(type(left) in (int, float) and type(right) in (int, float)
                and math.isfinite(left) and math.isfinite(right), f"{label} contains nonfinite output")
        values.append(abs(left-right))
    return max(values)


def compare_endpoints(full, split):
    """Compare actual final records without minimum-image hiding of drift."""
    require(full["step"] == split["step"], "Continuous and split endpoint steps differ")
    require(full["atom_count"] == split["atom_count"], "Restart endpoint atom counts differ")
    atoms_a = full["atoms"]
    atoms_b = split["atoms"]
    require(len(atoms_a) == len(atoms_b) == full["atom_count"], "Endpoint atom arrays differ")
    errors = {"position_component_angstrom": 0.0, "velocity_component_angstrom_fs": 0.0,
              "force_component_kj_mol_angstrom": 0.0}
    for left, right in zip(atoms_a, atoms_b):
        for key in ("engine_id", "cell_atom_id", "molecule_index", "charge_e", "image_flags"):
            require(left[key] == right[key], f"Endpoint {key} changed across restart")
        for key, result in (("coordinates_angstrom", "position_component_angstrom"),
                            ("wrapped_coordinates_angstrom", "position_component_angstrom"),
                            ("velocities_angstrom_fs", "velocity_component_angstrom_fs"),
                            ("forces_kj_mol_angstrom", "force_component_kj_mol_angstrom")):
            errors[result] = max(errors[result], _maximum_difference(left[key], right[key], key))
    errors["box_component_angstrom"] = max(
        _maximum_difference(full["box"][key], split["box"][key], "box "+key)
        for key in ("origin_angstrom", "lengths_angstrom"))
    errors["energy_kj_mol"] = max(
        _maximum_difference([full[key]], [split[key]], key)
        for key in ("potential_energy_kj_mol", "total_energy_kj_mol"))
    failures = [key for key, value in errors.items() if value > TOLERANCES[key]]
    require(not failures, "Predeclared split-restart endpoint parity failed: "+", ".join(failures))
    return {"status": "passed", "endpoint_step": full["step"],
            "maximum_absolute_errors": errors, "predeclared_absolute_tolerances": dict(TOLERANCES),
            "exact_atom_molecule_charge_image_identity": True,
            "interpretation": "Same-cell numerical restart continuity under this exact engine/configuration; not independent sampling or physical model validation"}


def copy_text_artifacts(run, run_directory, target):
    available = {}
    hashes = {"input.lammps": run["input_sha256"], **run.get("staged_input_sha256", {}),
              **run.get("artifact_sha256", {})}
    for name in TEXT_ARTIFACTS:
        source = run_directory/name
        if not source.is_file():
            continue
        raw = source.read_bytes()
        require(bool(raw), "Scientific text artifact is empty: "+name)
        raw.decode("utf-8")
        require(name in hashes and digest(raw) == hashes[name], "Scientific artifact hash changed: "+name)
        shutil.copyfile(source, target/name)
        available[name] = {"sha256": digest(raw), "bytes": len(raw)}
    require("input.lammps" in available and "log.lammps" in available and "checkpoint.dump" in available,
            "Actual engine input, log or velocity-bearing checkpoint dump missing")
    return available


def assess_continuation_sampling(result, assignment):
    """Analyze canonical NPT thermo rows; do not concatenate restart branches."""
    blocks = result["engine_run"]["thermo_blocks"]
    require(len(blocks) == 2 and len(blocks[0]["rows"]) == 1,
            "Continuation sampling requires one restore run-zero block and one NPT integration block")
    block = blocks[1]
    require(set(("Step", "Temp", "Press", "Volume")).issubset(block["columns"]),
            "Continuation block lacks physical sampling inputs")
    masses = [atom["mass_g_mol"] for atom in assignment["atoms"]]
    rows = []
    for values in block["rows"]:
        row = dict(zip(block["columns"], values))
        require(row["Step"] == int(row["Step"]), "Noninteger engine sampling step")
        rows.append({"step": int(row["Step"]), "temperature_k": row["Temp"],
                     "pressure_bar": row["Press"]*1.01325,
                     "density_g_cm3": density_from_mass_volume(masses, row["Volume"])})
    spec = DensitySamplingSpec(timestep_fs=PREPARATION["timestep_fs"])
    analysis = analyze_density_sampling(rows, spec)
    require(analysis["sampling_status"] == "insufficient_sampling" and analysis["production_ready"] is False,
            "The short continuation cannot pass the predeclared sampling review")
    return {"source_stage": "continuous_40", "canonical_thermo_rows": rows, "analysis": analysis,
            "interpretation": "One 10 fs continuation window from one prepared cell; split branches are excluded from independent-cell statistics"}


def main(argv=None):
    from polychemmd.protocols.restart import (
        PeriodicAlkaneContinuationSpec, run_periodic_alkane_continuation,
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/restart-v06")
    parser.add_argument("--executable")
    args = parser.parse_args(argv)
    modules = source_records()
    raw = (REPOSITORY/SOURCE_CELL_PATH).read_bytes()
    record = json.loads(raw)
    source = record["payload"] if "payload" in record else record
    system = audit_source_cell(source)
    assigned = assign_periodic_alkanes(source)
    assignment = assigned.to_dict()
    assignment_hash = digest(canonical_bytes(assignment))
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    private = output/"private_engine_runs"
    private.mkdir()
    os.environ["OMP_NUM_THREADS"], os.environ["UCX_TLS"] = "1", "self"
    started = time.monotonic()
    report = {
        "schema_version": "1", "evidence_kind": "executed_same_cell_native_npt_restart_audit",
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"polychemmd": __version__, "python": platform.python_version()},
        "system": system, "hardware": hardware_record(), "source_module_hashes": modules,
        "source_cell": {"path": SOURCE_CELL_PATH, "sha256": digest(raw)},
        "assignment_canonical_sha256": assignment_hash,
        "predeclared": {"preparation": dict(PREPARATION), "continuous_steps": FULL_STEPS,
                         "split_steps": [SPLIT_STEPS, SPLIT_STEPS],
                         "per_engine_timeout_seconds": PER_STAGE_TIMEOUT_SECONDS,
                         "total_engine_and_audit_wall_budget_seconds": TOTAL_WALL_BUDGET_SECONDS,
                         "maximum_engine_invocations": 4, "absolute_tolerances": dict(TOLERANCES)},
        "runtime_environment": {"OMP_NUM_THREADS": "1", "UCX_TLS": "self"},
        "stages": [], "equilibrium_status": "insufficient_sampling", "production_ready": False,
        "density_prediction_validation": "not_performed", "external_comparison": "not_performed",
        "scientific_limits": [
            "One dilute, four-decane starting cell is reused for every branch; branches and frames are correlated and are not independent preparation replicates.",
            "The longest physical path has only 60 fs integration including preparation; minimization iterations add no trajectory duration.",
            "Native checkpoint portability and trajectory continuity are restricted to the exact executable, linked runtime and configuration checked by the continuation protocol.",
            "Restart parity does not validate a liquid phase, polymer melt, force-field predictivity, sampling adequacy or any experimental density.",
            "No quantum calculation, parameter fit or external material-property comparison is performed.",
        ],
        "public_record_policy": "Publish report.json and stage text artifacts only. Binary checkpoints remain in private_engine_runs; their hashes, sizes and lineage are in public manifests. Regenerate checkpoints with this script on the supported matching engine/configuration.",
    }
    write_json(output/"report.json", report)
    def budget():
        remaining = TOTAL_WALL_BUDGET_SECONDS-(time.monotonic()-started)
        require(remaining > 0, "Predeclared total audit wall budget exhausted")
        return min(PER_STAGE_TIMEOUT_SECONDS, remaining)
    def retain(name, result, elapsed):
        directory = output/name
        directory.mkdir()
        public = public_record(result)
        write_json(directory/"protocol.json", public)
        artifacts = copy_text_artifacts(result["engine_run"], private/name, directory)
        report["stages"].append({"name": name, "status": result["status"],
                                 "measured_stage_wall_seconds": elapsed,
                                 "engine_elapsed_seconds": result["engine_run"]["elapsed_seconds"],
                                 "engine_version_banner": result["engine_run"]["engine_version_banner"],
                                 "engine_executable_sha256": result["engine_run"]["executable_sha256"],
                                 "checkpoint": result["checkpoint"], "text_artifacts": artifacts})
        write_json(output/"report.json", report)
    stage_name = "preparation"
    try:
        before = time.monotonic()
        spec = PeriodicAlkaneProtocolSpec.from_mapping({**PREPARATION, "timeout_seconds": budget()})
        prepared = run_periodic_alkane_protocol(assigned, spec, run_directory=private/stage_name,
                                                executable=args.executable)
        retain(stage_name, prepared, time.monotonic()-before)
        results = {}
        for name, parent, steps in (("continuous_40", "preparation", FULL_STEPS),
                                     ("split_20a", "preparation", SPLIT_STEPS),
                                     ("split_20b", "split_20a", SPLIT_STEPS)):
            stage_name = name
            previous = prepared if parent == "preparation" else results[parent]
            before = time.monotonic()
            continuation = PeriodicAlkaneContinuationSpec(steps=steps, timeout_seconds=budget(),
                                                         thermo_stride_steps=10)
            result = run_periodic_alkane_continuation(
                assigned, previous["checkpoint"], continuation,
                checkpoint_directory=private/parent, run_directory=private/name,
                executable=args.executable)
            results[name] = result
            retain(name, result, time.monotonic()-before)
        report["endpoint_comparison"] = compare_endpoints(results["continuous_40"]["final_state"],
                                                           results["split_20b"]["final_state"])
        report["density_sampling_assessment"] = assess_continuation_sampling(results["continuous_40"], assignment)
        report.update(status="passed", actual_engine_invocations=4,
                      integrated_steps_actually_executed_including_redundant_branches=280,
                      longest_physical_path_integration_steps=240,
                      longest_physical_path_duration_fs=240*spec.timestep_fs,
                      unique_continuation_path_duration_fs=40*spec.timestep_fs,
                      minimization_iterations_are_md_steps=False)
        require(time.monotonic()-started <= TOTAL_WALL_BUDGET_SECONDS,
                "Actual audit exceeded the predeclared total wall budget")
    except Exception as exc:
        report.update(status="failed", failure_stage=stage_name, failure_type=type(exc).__name__,
                      failure_reason=str(exc))
        if isinstance(getattr(exc, "record", None), dict):
            report["failed_protocol_record"] = public_record(exc.record)
    report.update(finished_utc=datetime.now(timezone.utc).isoformat(),
                  measured_total_wall_seconds=time.monotonic()-started)
    write_json(output/"report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output/"report.json"),
                      "measured_total_wall_seconds": report["measured_total_wall_seconds"],
                      "equilibrium_status": report["equilibrium_status"],
                      "endpoint_comparison": report.get("endpoint_comparison"),
                      "failure_reason": report.get("failure_reason")}, indent=2))
    return 0 if report["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
