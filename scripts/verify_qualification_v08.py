#!/usr/bin/env python3
"""Execute one predeclared numerical and preparation-readiness audit.

The construction is one ordered cell. Static refinement points and continuation
fragments are paired observations, never independent prepared-cell replicates.
Numerical review of its initial geometry cannot certify a deformed final state.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import shutil
import time

from polychemmd import __version__
from polychemmd.geometry import MolecularGeometry

REPOSITORY = Path(__file__).resolve().parents[1]
_HELPER_PATH = Path(__file__).with_name("verify_large_cell_v07.py")
_HELPER_SPEC = importlib.util.spec_from_file_location("polychemmd_v07_audit_helpers", _HELPER_PATH)
_HELPER = importlib.util.module_from_spec(_HELPER_SPEC)
_HELPER_SPEC.loader.exec_module(_HELPER)
# Independently written project helpers are reused; their source hash is retained.
digest = _HELPER.digest
canonical_bytes = _HELPER.canonical_bytes
write_json = _HELPER.write_json
hardware_record = _HELPER.hardware_record
motif_contact_reference = _HELPER.motif_contact_reference
retain_artifacts = _HELPER.retain_artifacts

SOURCE_GEOMETRY_PATH = "evidence/large_cell_v07/common/source_geometry.json"
RECIPE_PATH = "examples/qualification_cell.json"
RECIPE = {"repeats": [3, 7, 7], "spacing_angstrom": [13., 5.1, 5.1],
          "rotation_matrix": None, "min_distance_angstrom": 1.5,
          "max_atoms": 8192, "max_distance_checks": 5000000}
TOTAL_WALL_SECONDS = 240.0
PER_RUN_TIMEOUT_SECONDS = 60.0
MAX_ENGINE_INVOCATIONS = 8
NUMERICAL_SPEC = {
    "cutoffs_angstrom": [10., 12., 14., 16., 17.],
    "ewald_alpha_per_angstrom": .30, "reciprocal_nmax": 16,
    "refined_reciprocal_nmax": 20, "min_consecutive_cutoff_pairs": 2,
    "cutoff_energy_per_atom_limit_kj_mol": .02,
    "cutoff_force_limit_kj_mol_angstrom": .05, "cutoff_pressure_limit_bar": 20.,
    "reciprocal_energy_per_atom_limit_kj_mol": 1e-5,
    "reciprocal_force_limit_kj_mol_angstrom": 1e-4,
    "reciprocal_pressure_limit_bar": .01,
    "unchanged_component_abs_tolerance_kj_mol": 1e-7,
    "per_run_timeout_seconds": PER_RUN_TIMEOUT_SECONDS,
    "total_wall_seconds": 120., "max_engine_runs": 6,
}
PREPARATION = {"temperature_k": 300., "pressure_bar": 1., "timestep_fs": .25,
    "minimize_max_iterations": 10, "minimize_max_evaluations": 100,
    "nvt_steps": 10, "npt_steps": 10, "velocity_seed": 49217,
    "timeout_seconds": PER_RUN_TIMEOUT_SECONDS, "cutoff_angstrom": 17.,
    "ewald_alpha_per_angstrom": .30, "reciprocal_nmax": 20}
CONTINUATION_STEPS = 10
SOURCE_MODULES = (
    "polychemmd", "polychemmd.geometry", "polychemmd.builders.tiling",
    "polychemmd.forcefields.periodic", "polychemmd.engines.periodic_alkane",
    "polychemmd.engines.periodic_thermo", "polychemmd.engines.lammps",
    "polychemmd.protocols.numerics", "polychemmd.protocols.alkane",
    "polychemmd.protocols.restart", "polychemmd.analysis.structure",
    "polychemmd.analysis.readiness", "polychemmd.analysis.sampling",
)


class QualificationAuditError(RuntimeError):
    """An actual identity, resource or declared evidence invariant failed."""


def require(condition, message):
    if not condition:
        raise QualificationAuditError(message)


def public_summary(raw):
    """Keep static aggregate metrics while hash-identifying duplicated rows."""
    record = _HELPER.public_summary(raw)
    omitted = record["public_summary_policy"]["omitted_fields"]
    def visit(value, location="record"):
        if isinstance(value, dict):
            for key in list(value):
                if key == "per_chain" and isinstance(value[key], list):
                    rows = value.pop(key)
                    omitted.append({"field": location+"."+key,
                        "canonical_sha256": digest(canonical_bytes(rows)),
                        "row_count": len(rows),
                        "interpretation": "Repeated static chain rows omitted here; full continued structure has a separate public JSON and earlier rows can be regenerated"})
                elif key != "public_summary_policy":
                    visit(value[key], location+"."+key)
        elif isinstance(value, list):
            for index, item in enumerate(value):
                visit(item, location+f"[{index}]")
    visit(record)
    return record


def source_records():
    records = {}
    for name in SOURCE_MODULES:
        actual = Path(importlib.import_module(name).__file__).resolve()
        relative = Path("src") / Path(*name.split("."))
        expected = REPOSITORY / (relative/"__init__.py" if name == "polychemmd" else relative.with_suffix(".py"))
        require(actual.is_file() and expected.is_file(), "Use this editable source checkout")
        require(digest(actual.read_bytes()) == digest(expected.read_bytes()), "Imported source differs from checkout: "+name)
        records[expected.relative_to(REPOSITORY).as_posix()] = {"sha256": digest(actual.read_bytes())}
    for path in (Path(__file__), _HELPER_PATH):
        records[path.relative_to(REPOSITORY).as_posix()] = {"sha256": digest(path.read_bytes())}
    return records


def audit_cell(cell):
    atoms, molecules = cell["atoms"], cell["molecules"]
    require(len(atoms) == 4704 and len(molecules) == 147, "Declared cell must contain 147 decanes / 4704 atoms")
    require(Counter(atom["molecule_index"] for atom in atoms) == {index: 32 for index in range(147)},
            "Molecule/atom accounting differs from the declaration")
    for index in range(147):
        require(Counter(atom["element"] for atom in atoms if atom["molecule_index"] == index) == {"C": 10, "H": 22},
                "Constructed chemistry changed")
    box = cell["box"]["lengths_angstrom"]
    require(len(box) == 3 and all(math.isclose(a, b, rel_tol=0, abs_tol=1e-10)
            for a, b in zip(box, (39., 35.7, 35.7))), "Declared box changed")
    require(min(box) > 2*max(NUMERICAL_SPEC["cutoffs_angstrom"]), "Finest cutoff is not admissible")
    masses = [atom["mass_g_mol"] for atom in atoms]
    density = math.fsum(masses)/(0.602214076*math.prod(box))
    require(math.isclose(density, cell["construction_density_g_cm3"], rel_tol=1e-12), "Construction mass/volume density disagrees")
    return {"atoms": 4704, "molecules": 147, "formula_per_molecule": "C10H22",
        "box_lengths_angstrom": box, "construction_density_g_cm3": density,
        "density_role": "requested geometric mass/volume, not a prediction",
        "phase": "not_assigned", "construction": "one ordered repeated molecular motif cell",
        "independently_prepared_cells": 1, "molecules_and_variants_are_independent_cells": False}


def checked_preparation_review(review):
    """A short ordered probe must never be promoted to a material prediction."""
    require(review.get("production_ready") is False, "Readiness review promoted the technical pilot")
    require(review.get("equilibrium_status") == "not_assessed", "Readiness review invented an equilibrium certificate")
    require(review.get("status") in ("requires_preparation", "insufficient_sampling", "diagnostics_passed_but_unqualified"),
            "Readiness review status is outside its evidence contract")
    return review


def assert_sources_unchanged(records):
    for name, metadata in records.items():
        require(digest((REPOSITORY/name).read_bytes()) == metadata["sha256"], "Runtime source changed during audit: "+name)


def require_numerical_phase_budget(elapsed_seconds):
    """Do not start a fixed-budget child after spending its outer allowance."""
    require(type(elapsed_seconds) in (int, float) and math.isfinite(elapsed_seconds)
            and elapsed_seconds >= 0, "Elapsed audit time must be finite and nonnegative")
    remaining = TOTAL_WALL_SECONDS-elapsed_seconds
    require(remaining >= NUMERICAL_SPEC["total_wall_seconds"]+10.,
            "Remaining outer wall budget cannot cover the fixed numerical phase plus audit overhead")
    return remaining


def retain_failed_inputs(private, output):
    """Preserve bounded partial input/log text without claiming completion."""
    retained = {}
    for name in ("input.lammps", "log.lammps"):
        for source in sorted(private.rglob(name)):
            if source.is_symlink() or not source.resolve().is_relative_to(private.resolve()):
                continue
            raw = source.read_bytes()
            if not raw:
                continue
            relative = source.relative_to(private)
            try:
                raw.decode("utf-8")
            except UnicodeDecodeError:
                retained[relative.as_posix()] = {"sha256": digest(raw), "bytes": len(raw),
                    "public_path": None, "evidence": "non_UTF8_partial_file_hash_only_not_success"}
                continue
            target = output/"failed_inputs"/relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(source, target)
            retained[relative.as_posix()] = {"sha256": digest(raw), "bytes": len(raw),
                "public_path": target.relative_to(output).as_posix(),
                "evidence": "partial_current_attempt_text_not_success_or_convergence"}
    return retained


def main(argv=None):
    from polychemmd.builders.tiling import tile_molecular_geometry
    from polychemmd.forcefields.periodic import assign_compact_periodic_alkanes
    from polychemmd.protocols.numerics import NumericalQualificationSpec, plan_numerical_qualification, run_numerical_qualification
    from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, run_periodic_alkane_protocol
    from polychemmd.protocols.restart import PeriodicAlkaneContinuationSpec, run_periodic_alkane_continuation
    from polychemmd.analysis.structure import analyze_chain_structure
    from polychemmd.analysis.readiness import PreparationReadinessSpec, analyze_preparation_readiness
    from polychemmd.analysis.sampling import DensitySamplingSpec, analyze_density_sampling

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/qualification-v08")
    parser.add_argument("--executable")
    args = parser.parse_args(argv)
    modules = source_records()
    output = Path(args.output_dir).resolve(); output.mkdir(parents=True, exist_ok=False)
    common, private = output/"common", output/"private_engine_runs"
    common.mkdir(); private.mkdir()
    os.environ["OMP_NUM_THREADS"], os.environ["UCX_TLS"] = "1", "self"
    started = time.monotonic()
    readiness_spec = PreparationReadinessSpec()
    numerical_spec = NumericalQualificationSpec.from_mapping(NUMERICAL_SPEC)
    report = {"schema_version": "1", "evidence_kind": "executed_numerical_and_preparation_qualification_audit",
        "status": "running", "started_utc": datetime.now(timezone.utc).isoformat(),
        "software": {"polychemmd": __version__, "python": platform.python_version()},
        "hardware": hardware_record(), "source_module_hashes": modules,
        "predeclared": {"construction_recipe": RECIPE, "numerical_spec": numerical_spec.to_dict(),
            "preparation": PREPARATION, "continuation_steps": CONTINUATION_STEPS,
            "preparation_readiness_spec": readiness_spec.to_dict(),
            "maximum_engine_invocations": MAX_ENGINE_INVOCATIONS,
            "per_engine_timeout_seconds": PER_RUN_TIMEOUT_SECONDS, "total_wall_seconds": TOTAL_WALL_SECONDS},
        "stages": [], "engine_requests_started": 0,
        "equilibrium_status": "not_assessed", "production_ready": False,
        "density_prediction_validation": "not_performed", "external_comparison": "not_performed",
        "scientific_limits": [
            "One ordered repeated-motif cell, not an identified liquid, crystal or amorphous polymer",
            "Refinements change a truncated Hamiltonian or finite Ewald approximation, not an exact-model error bound",
            "Initial-geometry numerical diagnostics are not transferable to an evolved final state without a fresh audit",
            "Molecules, refinement settings and trajectory fragments are not independent prepared cells",
            "Short integration probes cannot certify equilibrium or bulk density"],
        "public_record_policy": "Retain source geometry/recipe, one common system.data, every stage input/log/summary, only finest refined singlepoint force frame, final continuation checkpoint frame, and one full continued-chain static dataset. Duplicated initial/prepared per-chain rows, full assignment, other raw frames/readbacks and native binaries are hash-only and regenerable with this script."}
    write_json(output/"report.json", report)

    def budget():
        report["engine_requests_started"] += 1
        require(report["engine_requests_started"] <= MAX_ENGINE_INVOCATIONS, "Declared process budget exhausted")
        remaining = TOTAL_WALL_SECONDS-(time.monotonic()-started)
        require(remaining > 10., "Aggregate audit budget exhausted before another engine request")
        write_json(output/"report.json", report)
        return min(PER_RUN_TIMEOUT_SECONDS, remaining-10.)

    def retain(name, result, source, elapsed, selected):
        source = Path(source).resolve()
        require(source.is_relative_to(private), "Engine artifacts escaped the current private audit tree")
        target = output/name; target.mkdir()
        run = result["engine_run"]
        artifacts = retain_artifacts(run, source, target, common, selected)
        write_json(target/"record.json", public_summary(result))
        report["stages"].append({"name": name, "status": result["status"],
            "measured_stage_wall_seconds": elapsed,
            "stage_wall_time_basis": ("not_measured_per_point; numerical phase total and engine elapsed are retained"
                if elapsed is None else "measured around this complete protocol call"),
            "engine_elapsed_seconds": run["elapsed_seconds"], "engine_version_banner": run["engine_version_banner"],
            "engine_executable_sha256": run["executable_sha256"], "engine_identity": run.get("engine_identity"),
            "raw_artifacts": artifacts})
        write_json(output/"report.json", report)

    stage = "construction"
    try:
        raw_geometry, raw_recipe = (REPOSITORY/SOURCE_GEOMETRY_PATH).read_bytes(), (REPOSITORY/RECIPE_PATH).read_bytes()
        require(json.loads(raw_recipe) == RECIPE, "Published construction recipe differs from predeclared settings")
        geometry = MolecularGeometry.from_mapping(json.loads(raw_geometry))
        (common/"source_geometry.json").write_bytes(raw_geometry); (common/"recipe.json").write_bytes(raw_recipe)
        report["source_geometry"] = {"path": SOURCE_GEOMETRY_PATH, "sha256": digest(raw_geometry),
            "recipe_path": RECIPE_PATH, "recipe_sha256": digest(raw_recipe)}
        before = time.monotonic(); cell = tile_molecular_geometry(geometry, **RECIPE).to_dict()
        report["construction_wall_seconds"] = time.monotonic()-before
        report["system"] = audit_cell(cell)
        report["independent_motif_contact_reference"] = motif_contact_reference(geometry, RECIPE)
        stage = "assignment"; before = time.monotonic()
        assigned = assign_compact_periodic_alkanes(cell); assignment = assigned.to_dict()
        encoded = canonical_bytes(assignment); (private/"compact_assignment.json").write_bytes(encoded)
        report["assignment"] = {"canonical_sha256": digest(encoded), "canonical_json_bytes": len(encoded),
            "wall_seconds": time.monotonic()-before, "model_id": assignment["model_id"],
            "schema_version": assignment["schema_version"], "special_pair_records": len(assignment["special_pairs"]),
            "implicit_all_unordered_pair_count": len(assignment["atoms"])*(len(assignment["atoms"])-1)//2,
            "global_nonbonded_pairs_materialized": "nonbonded_pairs" in assignment,
            "parameter_set_sha256": digest(canonical_bytes(assignment["parameter_set"]))}
        report["initial_chain_structure"] = public_summary(analyze_chain_structure(assigned))
        require(abs(report["initial_chain_structure"]["nematic_order_parameter"]-1.) < 1e-10,
                "The fixed initial motif must retain its explicitly ordered orientation")
        stage = "numerical_qualification"
        plan = plan_numerical_qualification(assigned, numerical_spec, run_directory=private/"numerical", executable=args.executable)
        report["numerical_plan"] = public_summary(plan); write_json(output/"report.json", report)
        require_numerical_phase_budget(time.monotonic()-started)
        before = time.monotonic()
        numerical = run_numerical_qualification(assigned, numerical_spec, run_directory=private/"numerical", executable=args.executable)
        report["numerical_phase_wall_seconds"] = time.monotonic()-before
        report["numerical_review"] = public_summary(numerical)
        for point in numerical["stages"]:
            report["engine_requests_started"] += 1
            require(report["engine_requests_started"] <= MAX_ENGINE_INVOCATIONS, "Numerical phase exceeded process budget")
            result = point["record"]
            selected = {"input.lammps", "log.lammps"}
            if point["stage_id"] == "reciprocal_refinement": selected.add("forces.dump")
            retain("numerical_"+point["stage_id"], result, point["run_directory"],
                   None, selected)
        require(len(numerical["stages"]) == 6, "Declared numerical sweep requires exactly six points")
        report["numerical_review_scope"] = {"configuration": "original immutable constructed geometry",
            "transferred_to_final_state": False, "bulk_convergence": "not_established"}
        stage = "preparation"; before = time.monotonic()
        spec = PeriodicAlkaneProtocolSpec.from_mapping({**PREPARATION, "timeout_seconds": budget()})
        prepared = run_periodic_alkane_protocol(assigned, spec, run_directory=private/stage, executable=args.executable)
        retain(stage, prepared, private/stage, time.monotonic()-before, {"input.lammps", "log.lammps"})
        report["preparation_minimizer"] = prepared["stages"]["minimization"]
        report["preparation_readiness"] = public_summary(checked_preparation_review(analyze_preparation_readiness(
            assigned, prepared, readiness_spec, numerical_review=numerical)))
        report["prepared_chain_structure"] = public_summary(analyze_chain_structure(assigned, prepared["final_state"]))
        stage = "continuation"; before = time.monotonic()
        continued = run_periodic_alkane_continuation(assigned, prepared["checkpoint"],
            PeriodicAlkaneContinuationSpec(steps=CONTINUATION_STEPS, thermo_stride_steps=10, timeout_seconds=budget()),
            checkpoint_directory=private/"preparation", run_directory=private/stage, executable=args.executable)
        retain(stage, continued, private/stage, time.monotonic()-before,
               {"input.lammps", "log.lammps", "checkpoint.dump"})
        continued_structure = analyze_chain_structure(assigned, continued["final_state"])
        structure_path = output/"continued_chain_structure.json"
        write_json(structure_path, continued_structure)
        report["continued_chain_structure"] = public_summary(continued_structure)
        report["continued_chain_structure_data"] = {"public_path": structure_path.name,
            "sha256": digest(structure_path.read_bytes()), "bytes": structure_path.stat().st_size,
            "canonical_sha256": digest(canonical_bytes(continued_structure)),
            "interpretation": "One static configuration; chain rows are not independent material cells or equilibrium observations"}
        rows = [{key: row[key] for key in ("step", "temperature_k", "pressure_bar", "density_g_cm3")}
                for row in continued["thermo_rows"]]
        review = analyze_density_sampling(rows, DensitySamplingSpec(timestep_fs=spec.timestep_fs))
        require(review["sampling_status"] == "insufficient_sampling", "Short technical continuation cannot satisfy sampling")
        report["density_sampling_assessment"] = {"canonical_thermo_rows": rows, "analysis": review}
        final = continued["final_state"]
        report["final_instantaneous_state"] = {key: final[key] for key in (
            "step", "atom_count", "box", "temperature_k", "pressure_bar", "volume_angstrom3",
            "density_g_cm3", "engine_density_g_cm3", "potential_energy_kj_mol", "total_energy_kj_mol")}
        report["restart_checks"] = {key: continued[key] for key in (
            "npt_fix_restart_state", "topology_parameter_unit_readback_audit", "initial_position_velocity_image_box_audit")}
        assert_sources_unchanged(modules)
        require(time.monotonic()-started <= TOTAL_WALL_SECONDS, "Audit exceeded declared aggregate wall budget")
        report.update(status="passed_engine_integration", actual_engine_invocations=len(report["stages"]),
            md_integration_steps_executed=spec.nvt_steps+spec.npt_steps+CONTINUATION_STEPS,
            physical_path_duration_fs=(spec.nvt_steps+spec.npt_steps+CONTINUATION_STEPS)*spec.timestep_fs,
            minimization_iterations_are_md_steps=False)
    except Exception as exc:
        report.update(status="failed", failure_stage=stage, failure_type=type(exc).__name__, failure_reason=str(exc))
        if isinstance(getattr(exc, "record", None), dict):
            report["failed_phase_record"] = public_summary(exc.record)
            if stage == "numerical_qualification" and type(exc.record.get("engine_requests_started")) is int:
                report["engine_requests_started"] = exc.record["engine_requests_started"]
        report["retained_failed_inputs"] = retain_failed_inputs(private, output)
    report.update(finished_utc=datetime.now(timezone.utc).isoformat(), measured_total_wall_seconds=time.monotonic()-started)
    write_json(output/"report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output/"report.json"),
        "measured_total_wall_seconds": report["measured_total_wall_seconds"],
        "numerical_status": report.get("numerical_review", {}).get("numerical_status"),
        "preparation_status": report.get("preparation_readiness", {}).get("status"),
        "failure_reason": report.get("failure_reason")}, indent=2))
    return 0 if report["status"] == "passed_engine_integration" else 1


if __name__ == "__main__":
    raise SystemExit(main())
