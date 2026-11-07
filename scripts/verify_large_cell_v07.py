#!/usr/bin/env python3
"""Execute a bounded dense-construction, sparse-assignment and engine audit.

This uses one ordered, tiled decane starting configuration. Its input density is
not a predicted liquid density. Numerical variants and restart fragments reuse
the same cell and are never counted as independent material samples.
"""

from __future__ import annotations

import argparse
from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import importlib
from itertools import product
import json
import math
import os
from pathlib import Path
import platform
import shutil
import time

from polychemmd import __version__
from polychemmd.analysis.sampling import DensitySamplingSpec, analyze_density_sampling
from polychemmd.geometry import MolecularGeometry


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE_CELL_PATH = "evidence/geometry_v03/cell.json"
RECIPE_PATH = "examples/tiled_decane.json"
RECIPE = {"repeats": [2, 6, 6], "spacing_angstrom": [13., 5.1, 5.1],
          "rotation_matrix": None, "min_distance_angstrom": 1.5,
          "max_atoms": 8192, "max_distance_checks": 5000000}
PER_STAGE_TIMEOUT_SECONDS = 60.0
TOTAL_WALL_BUDGET_SECONDS = 240.0
MAX_ENGINE_INVOCATIONS = 5
SETTINGS = {
    "baseline": {"cutoff_angstrom": 10., "ewald_alpha_per_angstrom": .30, "reciprocal_nmax": 16},
    "cutoff_12": {"cutoff_angstrom": 12., "ewald_alpha_per_angstrom": .30, "reciprocal_nmax": 16},
    "reciprocal_20": {"cutoff_angstrom": 10., "ewald_alpha_per_angstrom": .30, "reciprocal_nmax": 20},
}
# These review margins are declared before execution, not fitted to results.
# The cutoff variant changes the truncated Hamiltonian; its difference is not
# an implementation-error tolerance or an infinite-cutoff accuracy certificate.
SENSITIVITY_REVIEW_LIMITS = {
    "cutoff_12": {"total_energy_difference_per_atom_kj_mol": .02,
                   "max_force_component_difference_kj_mol_angstrom": .05},
    "reciprocal_20": {"total_energy_difference_per_atom_kj_mol": 1e-5,
                       "max_force_component_difference_kj_mol_angstrom": 1e-4},
}
UNCHANGED_COMPONENT_ABS_TOLERANCE_KJ_MOL = 1e-7
PREPARATION = {
    "temperature_k": 300., "pressure_bar": 1., "timestep_fs": .25,
    "minimize_max_iterations": 10, "minimize_max_evaluations": 100,
    "nvt_steps": 10, "npt_steps": 10, "velocity_seed": 49217,
    "timeout_seconds": PER_STAGE_TIMEOUT_SECONDS, **SETTINGS["baseline"],
}
CONTINUATION_STEPS = 10
COMPONENTS = ("bond", "angle", "dihedral", "lj", "coulomb_realspace", "coulomb_kspace", "total")
SOURCE_MODULES = (
    "polychemmd", "polychemmd.geometry", "polychemmd.builders.tiling",
    "polychemmd.forcefields.periodic", "polychemmd.engines.periodic_alkane",
    "polychemmd.engines.lammps", "polychemmd.protocols.alkane",
    "polychemmd.protocols.restart", "polychemmd.analysis.sampling",
)


class LargeCellAuditError(RuntimeError):
    """A declared geometry, identity, resource or engine invariant failed."""


def require(condition, message):
    if not condition:
        raise LargeCellAuditError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical_bytes(record):
    return json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False)+"\n", encoding="utf-8")


def finite(value, label):
    require(type(value) in (int, float), label+" must be numeric, not boolean")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise LargeCellAuditError(label+" is unrepresentable") from exc
    require(math.isfinite(number), label+" must be finite")
    return number


def source_records():
    records = {}
    for name in SOURCE_MODULES:
        actual = Path(importlib.import_module(name).__file__).resolve()
        relative = Path("src")/Path(*name.split("."))
        expected = REPOSITORY/(relative/"__init__.py" if name == "polychemmd" else relative.with_suffix(".py"))
        require(actual.is_file() and expected.is_file(), "Run this auditor from an editable source checkout")
        require(digest(actual.read_bytes()) == digest(expected.read_bytes()), "Imported source differs from checkout: "+name)
        records[expected.relative_to(REPOSITORY).as_posix()] = {"sha256": digest(actual.read_bytes())}
    records["scripts/verify_large_cell_v07.py"] = {"sha256": digest(Path(__file__).read_bytes())}
    return records


def hardware_record():
    model = platform.processor()
    cpuinfo = Path("/proc/cpuinfo")
    if not model and cpuinfo.is_file():
        for line in cpuinfo.read_text().splitlines():
            if line.startswith("model name"):
                model = line.split(":", 1)[1].strip()
                break
    return {"platform": platform.platform(), "machine": platform.machine(),
            "processor": model or None, "logical_cpu_count": os.cpu_count(),
            "engine_parallelism": "one direct engine process, one OpenMP thread; no GPU requested",
            "timing_interpretation": "Measured on this executing host, not a benchmark of the project owner's computer"}


def audit_constructed_cell(cell):
    """Check the declared large construction without inferring its phase."""
    atoms = cell["atoms"]
    require(len(atoms) == 2304 and len(cell["molecules"]) == 72, "Declared large cell must contain 72 decanes / 2304 atoms")
    require(Counter(atom["molecule_index"] for atom in atoms) == {index: 32 for index in range(72)},
            "Large cell molecule/atom accounting changed")
    for index in range(72):
        require(Counter(atom["element"] for atom in atoms if atom["molecule_index"] == index) == {"C": 10, "H": 22},
                "Large cell chemistry differs from neutral explicit decane")
    lengths = cell["box"]["lengths_angstrom"]
    for actual, expected in zip(lengths, (26., 30.6, 30.6)):
        require(abs(finite(actual, "box length")-expected) < 1e-10, "Construction box differs from fixed recipe")
    require(len(lengths) == 3 and min(lengths) > 2*max(v["cutoff_angstrom"] for v in SETTINGS.values()),
            "Construction cannot support the largest declared minimum-image cutoff")
    density = finite(cell["construction_density_g_cm3"], "construction density")
    require(.65 <= density <= .80, "Declared construction density range is not met")
    return {"atoms": len(atoms), "molecules": 72, "formula_per_molecule": "C10H22",
            "box_lengths_angstrom": lengths, "construction_density_g_cm3": density,
            "density_interpretation": "Mass divided by requested geometric volume; not an equilibrium prediction",
            "phase_observed": "not_assessed", "construction": "ordered repeated molecular motif",
            "independently_prepared_cells": 1, "molecules_are_independent_cell_replicates": False}


def motif_contact_reference(geometry, recipe):
    """Independent finite lattice check with a bound for all farther offsets.

    This exact recipe uses an identity rotation and at least two motifs per
    axis. Every offset in {-1,0,1}^3 except zero therefore refers to a distinct
    molecule or one of its periodic images. Any farther lattice offset has at
    least one component of magnitude two; its axis gap is bounded below by
    2*spacing minus the molecular coordinate span. This check is independent of
    the builder's linked contact bins and does not enumerate the full cell.
    """
    require(recipe["rotation_matrix"] is None and min(recipe["repeats"]) >= 2,
            "Motif reference supports this identity-oriented repeated recipe")
    points = geometry.coordinates_angstrom
    require(bool(points), "Motif requires coordinates")
    spacing = recipe["spacing_angstrom"]
    spans = [max(point[axis] for point in points)-min(point[axis] for point in points) for axis in range(3)]
    outside_lower_bound = min(2*step-span for step, span in zip(spacing, spans))
    distances = []
    for offset in product((-1, 0, 1), repeat=3):
        if offset == (0, 0, 0):
            continue
        shift = tuple(index*step for index, step in zip(offset, spacing))
        for first in points:
            for second in points:
                distance = math.dist(first, tuple(value+delta for value, delta in zip(second, shift)))
                require(math.isfinite(distance), "Motif contact is unrepresentable")
                distances.append(distance)
    minimum = min(distances)
    require(outside_lower_bound > minimum, "Farther lattice offsets are not excluded by this finite reference")
    require(minimum >= recipe["min_distance_angstrom"], "Independent motif reference found a forbidden contact")
    return {"minimum_periodic_intermolecular_distance_angstrom": minimum,
            "checked_nonzero_lattice_offsets": 26, "atom_pair_distance_evaluations": len(distances),
            "farther_offset_distance_lower_bound_angstrom": outside_lower_bound,
            "method": "26 neighboring motif offsets plus coordinate-span lower bound for all farther offsets",
            "construction_geometry_only": True}


def compare_singlepoints(baseline, variant, variant_name):
    """Compare the same ordered atoms under one changed numerical control."""
    require(variant_name in SENSITIVITY_REVIEW_LIMITS, "Unknown sensitivity comparison")
    ids = baseline["atom_ids"]
    require(isinstance(ids, list) and ids and len(set(ids)) == len(ids), "Singlepoint IDs must be unique")
    require(ids == variant["atom_ids"], "Sensitivity runs changed atom identity/order")
    require(baseline["assignment_canonical_sha256"] == variant["assignment_canonical_sha256"],
            "Sensitivity runs do not use the identical assigned geometry")
    energies_a, energies_b = baseline["energy_components_kj_mol"], variant["energy_components_kj_mol"]
    require(set(energies_a) == set(energies_b) == set(COMPONENTS), "Energy component fields changed")
    for energies in (energies_a, energies_b):
        total = finite(energies["total"], "total energy")
        require(abs(math.fsum(finite(energies[name], name) for name in COMPONENTS[:-1])-total)
                <= UNCHANGED_COMPONENT_ABS_TOLERANCE_KJ_MOL, "Singlepoint total energy does not match its components")
    deltas = {name: finite(energies_b[name], name)-finite(energies_a[name], name) for name in COMPONENTS}
    require(all(math.isfinite(value) for value in deltas.values()), "Energy difference is unrepresentable")
    unchanged = ("bond", "angle", "dihedral", "coulomb_kspace") if variant_name == "cutoff_12" else (
        "bond", "angle", "dihedral", "lj", "coulomb_realspace")
    require(all(abs(deltas[name]) <= UNCHANGED_COMPONENT_ABS_TOLERANCE_KJ_MOL for name in unchanged),
            "A fixed model component changed during the single-control comparison")
    first, second = baseline["forces_kj_mol_angstrom"], variant["forces_kj_mol_angstrom"]
    require(len(first) == len(second) == len(ids), "Singlepoint force count differs from atom identities")
    errors = []
    for row_a, row_b in zip(first, second):
        require(isinstance(row_a, (list, tuple)) and isinstance(row_b, (list, tuple)) and len(row_a) == len(row_b) == 3,
                "Singlepoint force vectors require three components")
        for a, b in zip(row_a, row_b):
            error = abs(finite(a, "force")-finite(b, "force"))
            require(math.isfinite(error), "Force difference is unrepresentable")
            errors.append(error)
    maximum = max(errors)
    rms = maximum*math.sqrt(math.fsum((value/maximum)**2 for value in errors)/len(errors)) if maximum else 0.
    per_atom = abs(deltas["total"])/len(ids)
    limits = SENSITIVITY_REVIEW_LIMITS[variant_name]
    review = (per_atom <= limits["total_energy_difference_per_atom_kj_mol"]
              and maximum <= limits["max_force_component_difference_kj_mol_angstrom"])
    return {"variant": variant_name, "same_assignment_and_atom_identity": True,
            "signed_component_energy_differences_kj_mol": deltas,
            "total_energy_difference_per_atom_kj_mol": per_atom,
            "max_force_component_difference_kj_mol_angstrom": maximum,
            "rms_force_component_difference_kj_mol_angstrom": rms,
            "predeclared_development_review_limits": limits,
            "review_status": "diagnostics_passed" if review else "requires_review",
            "unchanged_component_invariants": "passed",
            "interpretation": ("Two LJ/Coulomb real-space cutoffs change the truncated Hamiltonian; this is sensitivity, not error against an exact model"
                               if variant_name == "cutoff_12" else
                               "Two finite reciprocal bounds at fixed splitting/cutoff; not a bound on the infinite Ewald error"),
            "bulk_property_convergence": "not_established"}


def public_summary(raw):
    """Omit duplicate large arrays while retaining their explicit hash identity."""
    record = deepcopy(raw)
    omitted = []
    def visit(value, location="record"):
        if isinstance(value, dict):
            for key in list(value):
                if key in {"staged_inputs", "input_script", "source_cell", "forces_kj_mol_angstrom", "atom_ids"} or (key == "atoms" and isinstance(value[key], list)):
                    content = value.pop(key)
                    omitted.append({"field": location+"."+key, "canonical_sha256": digest(canonical_bytes(content)),
                                    "interpretation": "Large or duplicated field omitted; raw scientific artifacts have separate hashes"})
                elif key in {"run_directory", "checkpoint_directory"}:
                    value[key] = "private_engine_run_not_published"
                elif key == "executable" and isinstance(value[key], str):
                    value[key] = Path(value[key]).name
                elif key == "command" and isinstance(value[key], list) and value[key]:
                    value[key][0] = Path(value[key][0]).name
                elif key == "engine_help_first_line" and isinstance(value[key], str) and any(x in value[key] for x in ("UCX", "hwloc")):
                    value[key] = "Incidental host diagnostic omitted; engine version and hashes retained"
                else:
                    visit(value[key], location+"."+key)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, location+f"[{index}]")
    visit(record)
    record["public_summary_policy"] = {"omitted_fields": omitted,
        "scope": "Host paths and duplicate/large arrays omitted explicitly; scientific aggregate values, model, source, engine and raw artifact hashes retained"}
    return record


def retain_artifacts(run, source, target, common, selected):
    """Hash every raw artifact; copy only the declared public single-frame set."""
    hashes = {"input.lammps": run["input_sha256"], **run.get("staged_input_sha256", {}),
              **run.get("staged_binary_input_sha256", {}), **run.get("artifact_sha256", {})}
    metadata = {}
    for name, expected in hashes.items():
        path = source/name
        require(path.is_file(), "Recorded raw artifact is missing: "+name)
        raw = path.read_bytes()
        require(bool(raw) and digest(raw) == expected, "Raw artifact hash changed: "+name)
        item = {"sha256": expected, "bytes": len(raw), "public_path": None,
                "availability": "local_raw_artifact_hash_only; regenerate with the auditor"}
        if name == "system.data":
            destination = common/name
            if destination.exists():
                require(destination.read_bytes() == raw, "Singlepoint/preparation system data differ unexpectedly")
            else:
                destination.write_bytes(raw)
            item.update(public_path="common/system.data", availability="committed_text_input")
        elif name in selected:
            raw.decode("utf-8")
            shutil.copyfile(path, target/name)
            item.update(public_path=target.name+"/"+name, availability="committed_text_artifact")
        metadata[name] = item
    require(all(name in metadata for name in ("input.lammps", "log.lammps")), "Current input/log hashes missing")
    return metadata


def main(argv=None):
    from polychemmd.builders.tiling import tile_molecular_geometry
    from polychemmd.engines.periodic_alkane import run_periodic_alkane_singlepoint
    from polychemmd.forcefields.periodic import assign_compact_periodic_alkanes
    from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, run_periodic_alkane_protocol
    from polychemmd.protocols.restart import PeriodicAlkaneContinuationSpec, run_periodic_alkane_continuation
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/large-cell-v07")
    parser.add_argument("--executable")
    args = parser.parse_args(argv)
    modules = source_records()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    common, private = output/"common", output/"private_engine_runs"
    common.mkdir(); private.mkdir()
    os.environ["OMP_NUM_THREADS"], os.environ["UCX_TLS"] = "1", "self"
    started = time.monotonic()
    report = {
        "schema_version": "1", "evidence_kind": "executed_large_sparse_cell_engine_audit", "status": "running",
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"polychemmd": __version__, "python": platform.python_version()},
        "hardware": hardware_record(), "source_module_hashes": modules,
        "predeclared": {"construction_recipe": RECIPE, "singlepoint_settings": SETTINGS,
                         "preparation": PREPARATION, "continuation_steps": CONTINUATION_STEPS,
                         "maximum_engine_invocations": MAX_ENGINE_INVOCATIONS,
                         "per_engine_timeout_seconds": PER_STAGE_TIMEOUT_SECONDS,
                         "total_engine_and_audit_wall_budget_seconds": TOTAL_WALL_BUDGET_SECONDS,
                         "sensitivity_development_review_limits": SENSITIVITY_REVIEW_LIMITS,
                         "unchanged_component_abs_tolerance_kj_mol": UNCHANGED_COMPONENT_ABS_TOLERANCE_KJ_MOL},
        "stages": [], "engine_requests_started": 0,
        "equilibrium_status": "insufficient_sampling", "production_ready": False,
        "density_prediction_validation": "not_performed", "external_comparison": "not_performed",
        "scientific_limits": ["One ordered tiled starting cell, not a validated liquid, crystal or amorphous phase",
                              "Construction density is a requested input; the short trajectory does not predict bulk density",
                              "Repeated molecules, numerical variants and continuation fragments are not independent prepared cells",
                              "Large-cell energy/force execution has no independent analytical-reference comparison; small-cell representation tests supply separate implementation evidence",
                              "Two numerical choices do not establish infinite-cutoff, infinite-reciprocal or bulk-property convergence"],
        "public_record_policy": "Commit fixed source geometry/recipe, one shared system.data, each stage input/log, three singlepoint force dumps and the final continuation checkpoint.dump. Other large outputs and native binaries remain hash-only; regenerate with this script.",
        "resource_calibration": {
            "prior_actual_failed_attempt": "evidence/large_cell_v07/failed_attempt.json",
            "original_minimizer_request": {"max_iterations": 200, "max_evaluations": 1000},
            "observed_failure": "The original 60-second preparation timed out during CG after the log reached iteration30, before any MD",
            "new_minimizer_request": {"max_iterations": 10, "max_evaluations": 100},
            "rationale": "Reduce minimizer work to obtain a bounded technical preparation/restart probe; no force-convergence or production claim",
            "physical_model_settings_changed": False, "sensitivity_review_limits_changed": False,
            "integration_steps_changed": False,
        },
    }
    write_json(output/"report.json", report)
    def budget():
        report["engine_requests_started"] += 1
        require(report["engine_requests_started"] <= MAX_ENGINE_INVOCATIONS, "Declared engine invocation count exhausted")
        write_json(output/"report.json", report)
        remaining = TOTAL_WALL_BUDGET_SECONDS-(time.monotonic()-started)
        require(remaining > 10., "Predeclared audit wall budget exhausted before another process")
        return min(PER_STAGE_TIMEOUT_SECONDS, remaining-10.)
    def retain(name, result, elapsed, selected):
        target = output/name; target.mkdir()
        run = result["engine_run"]
        artifacts = retain_artifacts(run, private/name, target, common, selected)
        write_json(target/"record.json", public_summary(result))
        report["stages"].append({"name": name, "status": result["status"],
                                  "measured_stage_wall_seconds": elapsed,
                                  "engine_elapsed_seconds": run["elapsed_seconds"],
                                  "engine_version_banner": run["engine_version_banner"],
                                  "engine_executable_sha256": run["executable_sha256"],
                                  "engine_identity": run.get("engine_identity"), "raw_artifacts": artifacts})
        write_json(output/"report.json", report)
    stage = "construction"
    try:
        raw_source = (REPOSITORY/SOURCE_CELL_PATH).read_bytes()
        original = json.loads(raw_source); original = original.get("payload", original)
        geometry = MolecularGeometry.from_mapping(original["molecules"][0]["source_geometry"])
        raw_recipe = (REPOSITORY/RECIPE_PATH).read_bytes()
        require(json.loads(raw_recipe) == RECIPE, "Published construction recipe differs from predeclared auditor settings")
        (common/"source_geometry.json").write_bytes(json.dumps(geometry.to_dict(), indent=2).encode()+b"\n")
        (common/"recipe.json").write_bytes(raw_recipe)
        report["source_geometry"] = {"source_path": SOURCE_CELL_PATH, "source_file_sha256": digest(raw_source),
                                      "molecule_index": 0, "canonical_geometry_sha256": digest(canonical_bytes(geometry.to_dict())),
                                      "recipe_path": RECIPE_PATH, "recipe_file_sha256": digest(raw_recipe)}
        before = time.monotonic(); cell = tile_molecular_geometry(geometry, **RECIPE).to_dict()
        report["construction_wall_seconds"] = time.monotonic()-before
        report["system"] = audit_constructed_cell(cell)
        report["independent_motif_contact_reference"] = motif_contact_reference(geometry, RECIPE)
        report["construction_canonical_sha256"] = digest(canonical_bytes(cell))
        stage = "compact_assignment"; before = time.monotonic()
        assigned = assign_compact_periodic_alkanes(cell)
        assignment = assigned.to_dict()
        assignment_bytes = canonical_bytes(assignment)
        (private/"compact_assignment.json").write_bytes(assignment_bytes)
        report["assignment"] = {"canonical_sha256": digest(assignment_bytes), "canonical_json_bytes": len(assignment_bytes),
                                "wall_seconds": time.monotonic()-before, "schema_version": assignment["schema_version"],
                                "model_id": assignment["model_id"], "special_pair_records": len(assignment["special_pairs"]),
                                "implicit_all_unordered_pair_count": len(assignment["atoms"])*(len(assignment["atoms"])-1)//2,
                                "global_nonbonded_pairs_materialized": "nonbonded_pairs" in assignment,
                                "parameter_set_sha256": digest(canonical_bytes(assignment["parameter_set"])),
                                "source_snapshot": assignment["source_snapshot"], "validation": assignment["validation"]}
        require("nonbonded_pairs" not in assignment and len(assignment["atoms"]) > 512, "Large assignment did not remove the old global-pair bottleneck")
        write_json(output/"report.json", report)
        results = {}
        for name, settings in SETTINGS.items():
            stage = name; before = time.monotonic()
            result = run_periodic_alkane_singlepoint(assigned, run_directory=private/name,
                executable=args.executable, timeout_seconds=budget(), **settings)
            results[name] = result
            retain(name, result, time.monotonic()-before, {"input.lammps", "log.lammps", "forces.dump"})
        report["numerical_sensitivity"] = {name: compare_singlepoints(results["baseline"], results[name], name)
                                            for name in ("cutoff_12", "reciprocal_20")}
        report["numerical_sensitivity_review_status"] = ("diagnostics_passed" if all(
            value["review_status"] == "diagnostics_passed" for value in report["numerical_sensitivity"].values()) else "requires_review")
        stage = "preparation"; before = time.monotonic()
        spec = PeriodicAlkaneProtocolSpec.from_mapping({**PREPARATION, "timeout_seconds": budget()})
        prepared = run_periodic_alkane_protocol(assigned, spec, run_directory=private/stage, executable=args.executable)
        retain(stage, prepared, time.monotonic()-before, {"input.lammps", "log.lammps"})
        report["preparation_minimizer"] = prepared["stages"]["minimization"]
        stage = "continuation"; before = time.monotonic()
        result = run_periodic_alkane_continuation(assigned, prepared["checkpoint"],
            PeriodicAlkaneContinuationSpec(steps=CONTINUATION_STEPS, thermo_stride_steps=10, timeout_seconds=budget()),
            checkpoint_directory=private/"preparation", run_directory=private/stage, executable=args.executable)
        retain(stage, result, time.monotonic()-before, {"input.lammps", "log.lammps", "checkpoint.dump"})
        rows = [{key: row[key] for key in ("step", "temperature_k", "pressure_bar", "density_g_cm3")}
                for row in result["thermo_rows"]]
        review = analyze_density_sampling(rows, DensitySamplingSpec(timestep_fs=spec.timestep_fs))
        require(review["sampling_status"] == "insufficient_sampling", "Short large-cell continuation cannot establish sufficient sampling")
        report["density_sampling_assessment"] = {"canonical_thermo_rows": rows, "analysis": review}
        final = result["final_state"]
        report["final_instantaneous_state"] = {key: final[key] for key in (
            "step", "atom_count", "box", "temperature_k", "pressure_bar", "volume_angstrom3",
            "density_g_cm3", "engine_density_g_cm3", "potential_energy_kj_mol", "total_energy_kj_mol")}
        report["restart_checks"] = {key: result[key] for key in ("npt_fix_restart_state", "topology_parameter_unit_readback_audit", "initial_position_velocity_image_box_audit")}
        for name, metadata in modules.items():
            require(digest((REPOSITORY/name).read_bytes()) == metadata["sha256"], "Runtime source changed during audit: "+name)
        require(time.monotonic()-started <= TOTAL_WALL_BUDGET_SECONDS, "Actual audit exceeded its declared total wall budget")
        report.update(status="passed_engine_integration", actual_engine_invocations=5,
                      md_integration_steps_executed=spec.nvt_steps+spec.npt_steps+CONTINUATION_STEPS,
                      physical_path_duration_fs=(spec.nvt_steps+spec.npt_steps+CONTINUATION_STEPS)*spec.timestep_fs,
                      minimization_iterations_are_md_steps=False)
    except Exception as exc:
        report.update(status="failed", failure_stage=stage, failure_type=type(exc).__name__, failure_reason=str(exc))
        if isinstance(getattr(exc, "record", None), dict):
            report["failed_protocol_record"] = public_summary(exc.record)
    report.update(finished_utc=datetime.now(timezone.utc).isoformat(), measured_total_wall_seconds=time.monotonic()-started)
    write_json(output/"report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output/"report.json"),
                      "measured_total_wall_seconds": report["measured_total_wall_seconds"],
                      "numerical_sensitivity_review_status": report.get("numerical_sensitivity_review_status"),
                      "failure_reason": report.get("failure_reason")}, indent=2))
    return 0 if report["status"] == "passed_engine_integration" else 1


if __name__ == "__main__":
    raise SystemExit(main())
