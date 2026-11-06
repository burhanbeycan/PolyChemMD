#!/usr/bin/env python3
"""Execute the bounded isolated OPLS-AA alkane energy/force verification.

Use an editable installation of this checkout with chemistry and LAMMPS:
``python scripts/verify_alkane_v04.py --output-dir artifacts/alkane-v04``.
The output directory must be new. This runs LAMMPS ``run 0`` single points;
it does not integrate trajectories, relax structures or validate a bulk model.
Public records are separated from private engine process records.
"""

from __future__ import annotations

from collections import Counter
from copy import deepcopy
from datetime import datetime, timezone
import argparse
import hashlib
import importlib
from importlib.resources import files
import json
import math
import os
from pathlib import Path
import platform
import re
import shutil
import time

from polychemmd import __version__
from polychemmd.builders.conformers import generate_conformer
from polychemmd.builders.linear import build_linear
from polychemmd.chemistry import RepeatUnitSpec
from polychemmd.engines.alkane import verify_alkane_singlepoint
from polychemmd.forcefields.alkane_data import MODEL_ID, load_alkane_parameters
from polychemmd.forcefields.alkanes import assign_alkane
from polychemmd.provenance import envelope


REPOSITORY = Path(__file__).resolve().parents[1]
RECIPE_PATH = REPOSITORY / "examples" / "alkane_validation.json"
ENERGY_SUM_TOLERANCE_KJ_MOL = 1e-9
NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM = 1e-7
NET_TORQUE_TOLERANCE_KJ_MOL = 1e-7
PUBLIC_ENGINE_ARTIFACTS = ("input.lammps", "system.data", "log.lammps", "forces.dump")
SOURCE_MODULES = (
    "polychemmd", "polychemmd.geometry", "polychemmd.provenance",
    "polychemmd.chemistry.contracts", "polychemmd.builders.linear",
    "polychemmd.builders.conformers", "polychemmd.forcefields.alkane_data",
    "polychemmd.forcefields.alkanes", "polychemmd.forcefields.contracts",
    "polychemmd.forcefields.reference", "polychemmd.descriptors.torsion",
    "polychemmd.engines.alkane", "polychemmd.engines.lammps",
)


class AlkaneAuditError(RuntimeError):
    """The actual assigned model or numerical evidence failed its audit."""


def require(condition, message):
    # Explicit checks remain enabled under python -O.
    if not condition:
        raise AlkaneAuditError(message)


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def canonical_bytes(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def _source_records():
    records = {}
    for name in SOURCE_MODULES:
        module = importlib.import_module(name)
        actual = Path(module.__file__).resolve()
        relative = Path("src") / Path(*name.split("."))
        expected = REPOSITORY / (relative / "__init__.py" if name == "polychemmd"
                                 else relative.with_suffix(".py"))
        require(actual.is_file() and expected.is_file(),
                "The auditor requires accessible source files from this prepared checkout")
        actual_hash = digest(actual.read_bytes())
        require(actual_hash == digest(expected.read_bytes()),
                f"Imported {name} differs from this checkout; install it in editable mode")
        records[expected.relative_to(REPOSITORY).as_posix()] = {"sha256": actual_hash}
    records["scripts/verify_alkane_v04.py"] = {"sha256": digest(Path(__file__).read_bytes())}
    return records


def _recipe(raw):
    def reject(value):
        raise ValueError("Validation recipe must contain finite JSON values")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("Validation recipe contains duplicate JSON fields")
            result[key] = value
        return result
    data = json.loads(raw, parse_constant=reject, object_pairs_hook=unique)
    require(data["schema_version"] == "1" and data["model_id"] == MODEL_ID,
            "Unsupported alkane validation recipe/model")
    require(isinstance(data["cases"], list) and len(data["cases"]) == 3,
            "This bounded auditor requires the three declared alkane cases")
    names = set()
    for case in data["cases"]:
        name = case["name"]
        require(isinstance(name, str) and re.fullmatch(r"[A-Za-z0-9_-]+", name) is not None
                and name not in names, "Case names must be unique safe directory labels")
        names.add(name)
        dp = case["degree_of_polymerization"]
        require(type(dp) is int and 1 <= dp <= 20, "Auditor DP must be an integer in [1, 20]")
        seed = case["seed"]
        require(type(seed) is int and 0 <= seed <= 2**31-1, "Invalid conformer seed")
    require(len({case["seed"] for case in data["cases"]}) == len(data["cases"]),
            "Each case needs its own explicit conformer seed")
    for name in ("cutoff_angstrom", "timeout_seconds_per_engine",
                 "energy_abs_tol_kj_mol", "force_abs_tol_kj_mol_angstrom"):
        value = data[name]
        require(type(value) in (int, float) and math.isfinite(float(value)) and value > 0,
                f"{name} must be positive and finite")
    require(data["timeout_seconds_per_engine"] <= 120,
            "Auditor engine timeout must not exceed 120 seconds per case")
    return data


def _audit_assignment(record, geometry, carbon_count):
    require(record["model_id"] == MODEL_ID, "Assigned model identity changed")
    require(record["geometry"] == geometry.to_dict(), "Assignment changed the source geometry")
    rows = record["atom_assignments"]
    expected_atoms = 3*carbon_count + 2
    require(len(rows) == len(geometry.atoms) == expected_atoms, "n-alkane atom balance failed")
    require([row["atom_id"] for row in rows] == [atom.atom_id for atom in geometry.atoms],
            "Assignment changed source atom identities or ordering")
    types = Counter(row["type_id"] for row in rows)
    expected_types = Counter({"opls_135": 2, "opls_140": 2*carbon_count+2})
    if carbon_count > 2:
        expected_types["opls_136"] = carbon_count-2
    require(types == expected_types, "CH3/CH2/HC type populations disagree with an n-alkane")
    counts = {kind: len(record["interactions"][kind])
              for kind in ("bonds", "angles", "dihedrals", "impropers")}
    require(counts == {"bonds": 3*carbon_count+1, "angles": 6*carbon_count,
                       "dihedrals": 9*(carbon_count-1), "impropers": 0},
            "Complete explicit n-alkane topology count check failed")
    require(len(record["nonbonded_pairs"]) == expected_atoms*(expected_atoms-1)//2,
            "Nonbonded enumeration does not contain every unordered atom pair")
    partial_charge = math.fsum(row["charge_e"] for row in rows)
    formal_charge = sum(atom.formal_charge for atom in geometry.atoms)
    require(abs(partial_charge) <= 1e-12 and formal_charge == 0,
            "Assigned partial/formal charges do not give a neutral alkane")
    require(record["coverage"]["assigned_type_coverage_complete"] is True,
            "Assigned parameter coverage is incomplete")
    require(record["validation"]["production_ready"] is False,
            "Assignment must not claim a validated bulk production model")
    return {"status": "passed", "carbon_count": carbon_count,
            "atom_count": expected_atoms, "interaction_counts": counts,
            "atom_type_counts": dict(sorted(types.items())),
            "total_partial_charge_e": partial_charge, "formal_charge_e": formal_charge,
            "all_unordered_pairs_enumerated": True}


def _balance_checks(energies, forces, coordinates):
    require(len(forces) == len(coordinates) and all(len(row) == 3 for row in forces),
            "Force and coordinate arrays disagree")
    require(all(math.isfinite(v) for row in forces for v in row), "Nonfinite force")
    total_from_components = math.fsum(value for name, value in energies.items() if name != "total")
    energy_error = abs(total_from_components - energies["total"])
    require(energy_error <= ENERGY_SUM_TOLERANCE_KJ_MOL,
            "Energy components do not sum to the reported total")
    net_force = [math.fsum(row[axis] for row in forces) for axis in range(3)]
    centroid = [math.fsum(point[axis] for point in coordinates)/len(coordinates)
                for axis in range(3)]
    torques = []
    for point, force in zip(coordinates, forces):
        x, y, z = (v-c for v, c in zip(point, centroid))
        fx, fy, fz = force
        torques.append((y*fz-z*fy, z*fx-x*fz, x*fy-y*fx))
    net_torque = [math.fsum(row[axis] for row in torques) for axis in range(3)]
    require(math.hypot(*net_force) <= NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM,
            "Isolated internal forces do not balance")
    require(math.hypot(*net_torque) <= NET_TORQUE_TOLERANCE_KJ_MOL,
            "Isolated internal force torque does not balance")
    return {"energy_component_sum_error_kj_mol": energy_error,
            "net_force_kj_mol_angstrom": net_force,
            "net_force_norm_kj_mol_angstrom": math.hypot(*net_force),
            "net_torque_about_centroid_kj_mol": net_torque,
            "net_torque_norm_kj_mol": math.hypot(*net_torque)}


def _public_parity(raw):
    record = deepcopy(raw)
    run, changes = record["engine_run"], []
    if "run_directory" in run:
        run["run_directory"] = "private_engine_run_not_in_public_records"
        changes.append("engine_run.run_directory replaced with a non-host label")
    if "executable" in run:
        run["executable"] = Path(run["executable"]).name
        changes.append("engine_run.executable reduced to basename; executable hash retained")
    if "command" in run and run["command"]:
        run["command"][0] = Path(run["command"][0]).name
        changes.append("engine_run.command executable reduced to basename")
    first_line = run.get("engine_help_first_line")
    if isinstance(first_line, str) and ("UCX" in first_line or "hwloc" in first_line):
        run["engine_help_first_line"] = "Incidental host diagnostic omitted; see engine_version_banner"
        changes.append("incidental engine-help host diagnostic omitted")
    record["public_record_sanitization"] = {
        "changes": changes,
        "scope": "Only host-specific paths and incidental help diagnostics are replaced. Scientific inputs, actual engine identity/version/hash and raw artifact hashes are retained.",
        "raw_record_location": "separate private_engine_runs case directory; not part of public case records",
    }
    return record


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/alkane-v04")
    parser.add_argument("--executable")
    args = parser.parse_args(argv)
    raw_recipe = RECIPE_PATH.read_bytes()
    recipe = _recipe(raw_recipe)
    modules = _source_records()
    parameter_raw = files("polychemmd.forcefields").joinpath("alkane_parameters.json").read_bytes()
    parameter_data = load_alkane_parameters()
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    private_root = output / "private_engine_runs"
    private_root.mkdir()
    # Deliberate single-rank controls for this bounded verification process.
    os.environ["OMP_NUM_THREADS"] = "1"
    os.environ["UCX_TLS"] = "self"
    started_utc, started = datetime.now(timezone.utc).isoformat(), time.monotonic()
    summaries = []
    for case in recipe["cases"]:
        case_started = time.monotonic()
        directory = output / case["name"]
        directory.mkdir()
        chain = build_linear(RepeatUnitSpec.from_mapping(recipe["repeat_unit"]),
                             case["degree_of_polymerization"],
                             head_cap=recipe["head_cap"], tail_cap=recipe["tail_cap"])
        geometry = generate_conformer(chain, seed=case["seed"])
        geometry_record = envelope(geometry.to_dict(), raw_recipe, "molecular_conformer")
        geometry_record["case_name"] = case["name"]
        geometry_record["case_recipe"] = case
        write_json(directory / "geometry.json", geometry_record)
        assigned = assign_alkane(geometry)
        assignment = assigned.to_dict()
        assignment_record = envelope(assignment, (directory / "geometry.json").read_bytes(),
                                     "assigned_opls_aa_n_alkane")
        assignment_record["case_name"] = case["name"]
        write_json(directory / "assignment.json", assignment_record)
        carbon_count = sum(atom.element == "C" for atom in geometry.atoms)
        topology_audit = _audit_assignment(assignment, geometry, carbon_count)
        run_directory = private_root / case["name"]
        parity = verify_alkane_singlepoint(
            assigned, run_directory=run_directory, executable=args.executable,
            cutoff_angstrom=recipe["cutoff_angstrom"],
            timeout_seconds=recipe["timeout_seconds_per_engine"],
            energy_abs_tol_kj_mol=recipe["energy_abs_tol_kj_mol"],
            force_abs_tol_kj_mol_angstrom=recipe["force_abs_tol_kj_mol_angstrom"])
        require(parity["status"] == "passed" and parity["md_integration_steps"] == 0,
                "Actual engine single-point parity did not pass")
        require(parity["assignment_canonical_sha256"] == digest(canonical_bytes(assignment)),
                "Parity assignment hash differs from the actual assigned model")
        backend_diagnostics = parity["backend_numerical_diagnostics"]
        require(backend_diagnostics["reference_or_potential_modified"] is False
                and backend_diagnostics["parity_tolerances_modified"] is False,
                "Backend diagnostics must not modify the reference model or tolerances")
        checks = {
            "reference": _balance_checks(parity["reference"]["energy_components_kj_mol"],
                                         parity["reference"]["forces_kj_mol_angstrom"],
                                         geometry.coordinates_angstrom),
            "engine": _balance_checks(parity["engine_energy_components_kj_mol"],
                                      parity["engine_forces_kj_mol_angstrom"],
                                      geometry.coordinates_angstrom),
        }
        public = _public_parity(parity)
        public["independent_balance_checks"] = checks
        write_json(directory / "parity.json", public)
        for name in PUBLIC_ENGINE_ARTIFACTS:
            source = run_directory / name
            require(source.is_file() and source.stat().st_size > 0,
                    f"Actual engine artifact missing: {name}")
            shutil.copyfile(source, directory / name)
        for name, expected in parity["engine_run"]["artifact_sha256"].items():
            require(digest((directory / name).read_bytes()) == expected,
                    f"Copied actual engine artifact changed: {name}")
        for name in ("input.lammps", "system.data"):
            require(digest((directory/name).read_bytes()) == parity["input_files"][name+"_sha256"],
                    f"Copied engine scientific input changed: {name}")
        artifacts = {path.name: {"sha256": digest(path.read_bytes()), "bytes": path.stat().st_size}
                     for path in sorted(directory.iterdir()) if path.is_file()}
        summaries.append({
            "case": case, "molecular_formula": chain.molecular_formula,
            "status": "passed", "elapsed_seconds": time.monotonic()-case_started,
            "topology_audit": topology_audit, "independent_balance_checks": checks,
            "source_graph_sha256": geometry.source_graph_sha256,
            "geometry_canonical_payload_sha256": digest(canonical_bytes(geometry.to_dict())),
            "assignment_canonical_payload_sha256": digest(canonical_bytes(assignment)),
            "public_parity_canonical_payload_sha256": digest(canonical_bytes(public)),
            "energy_component_abs_errors_kj_mol": parity["energy_component_abs_errors_kj_mol"],
            "maximum_energy_component_abs_error_kj_mol": max(parity["energy_component_abs_errors_kj_mol"].values()),
            "max_force_component_abs_error_kj_mol_angstrom": parity["max_force_component_abs_error_kj_mol_angstrom"],
            "rms_force_component_error_kj_mol_angstrom": parity["rms_force_component_error_kj_mol_angstrom"],
            "engine_version_banner": parity["engine_run"]["engine_version_banner"],
            "engine_executable_sha256": parity["engine_run"]["executable_sha256"],
            "near_planar_dihedral_count": len(backend_diagnostics["near_planar_dihedrals"]),
            "observed_source_sine_threshold": backend_diagnostics["observed_source_sine_threshold"],
            "backend_diagnostic_source_url": backend_diagnostics["pinned_source_url"],
            "backend_diagnostic_threshold_scope": backend_diagnostics["threshold_scope"],
            "md_integration_steps": 0, "public_case_directory": case["name"],
            "artifact_files": artifacts,
        })
    from rdkit import __version__ as rdkit_version
    report = {
        "schema_version": "1", "evidence_kind": "executed_isolated_n_alkane_energy_force_parity",
        "status": "passed", "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic()-started,
        "command": "python scripts/verify_alkane_v04.py --output-dir NEW_DIRECTORY [--executable LAMMPS]",
        "software": {"polychemmd": __version__, "python": platform.python_version(),
                     "python_implementation": platform.python_implementation(), "rdkit": rdkit_version},
        "runtime_environment": {"OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"], "UCX_TLS": os.environ["UCX_TLS"]},
        "input_recipe": {"path": "examples/alkane_validation.json", "sha256": digest(raw_recipe),
                         "hash_basis": "exact checked-in recipe-file bytes", "record": recipe},
        "source_module_hashes": modules,
        "parameter_json_sha256": digest(parameter_raw),
        "parameter_set_canonical_sha256": digest(canonical_bytes(parameter_data["parameter_set"])),
        "source_snapshot": parameter_data["source_snapshot"],
        "canonical_payload_hash_basis": "sorted-key compact finite JSON encoded as UTF-8",
        "predeclared_balance_tolerances": {
            "energy_component_sum_kj_mol": ENERGY_SUM_TOLERANCE_KJ_MOL,
            "net_force_norm_kj_mol_angstrom": NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM,
            "net_torque_norm_kj_mol": NET_TORQUE_TOLERANCE_KJ_MOL,
        },
        "cases": summaries,
        "maximum_energy_component_abs_error_kj_mol": max(case["maximum_energy_component_abs_error_kj_mol"] for case in summaries),
        "max_force_component_abs_error_kj_mol_angstrom": max(case["max_force_component_abs_error_kj_mol_angstrom"] for case in summaries),
        "public_record_policy": "Only per-case public directories and this report are intended for publication; private_engine_runs retains unsanitized process records separately. Public parity replaces host paths and incidental help diagnostics, retaining scientific inputs, versions and hashes. Raw input/data/log/force artifacts are copied byte-for-byte.",
        "md_integration_steps": 0, "trajectory_production_performed": False,
        "quantum_calculation_performed": False, "equilibrium_status": "not_assessed",
        "bulk_production_model_validated": False, "material_property_validation": "not_performed",
        "scientific_limits": [
            "The attributed frozen subset is assigned only to supported neutral explicit n-alkanes.",
            "Actual isolated nonperiodic finite-cutoff LAMMPS run-zero energies/forces are compared; no trajectory is integrated.",
            "Energy/force agreement validates implementation parity, not force-field physical accuracy or polymer transferability.",
            "LAMMPS OPLS uses a near-planar sine derivative safeguard. Actual residuals are retained and the analytical reference, geometry, seeds and declared tolerances are not changed to imitate that safeguard; case parity records preserve its diagnostic context.",
            "Periodic bulk electrostatics, dense packing, relaxation, equilibration and physical properties remain separate milestones.",
        ],
    }
    write_json(output / "report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output/"report.json"),
                      "cases": [case["case"]["name"] for case in summaries],
                      "maximum_energy_component_abs_error_kj_mol": report["maximum_energy_component_abs_error_kj_mol"],
                      "max_force_component_abs_error_kj_mol_angstrom": report["max_force_component_abs_error_kj_mol_angstrom"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
