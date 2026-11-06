#!/usr/bin/env python3
"""Execute and audit the bounded periodic n-alkane preparation example.

Use an editable checkout with a working LAMMPS executable. The output directory
must be new. Public scientific records are separate from private process and
full-assignment records. Passing this audit does not certify equilibration.
"""

from __future__ import annotations

import argparse
from collections import Counter, deque
from copy import deepcopy
from datetime import datetime, timezone
import hashlib
import importlib
from importlib.resources import files
import json
import math
import os
from pathlib import Path
import platform
import shutil
import time

from polychemmd import __version__
from polychemmd.forcefields.alkane_data import MODEL_ID, load_alkane_parameters
from polychemmd.forcefields.periodic import assign_periodic_alkanes
from polychemmd.engines.periodic_alkane import verify_periodic_alkane_singlepoint
from polychemmd.protocols.alkane import PeriodicAlkaneProtocolSpec, run_periodic_alkane_protocol
from polychemmd.provenance import envelope


REPOSITORY = Path(__file__).resolve().parents[1]
RECIPE_PATH = REPOSITORY / "examples" / "periodic_validation.json"
COMPONENTS = ("bond", "angle", "dihedral", "lj", "coulomb_realspace", "coulomb_kspace")
ENERGY_SUM_TOLERANCE_KJ_MOL = 1e-8
NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM = 1e-7
TRANSLATION_ENERGY_TOLERANCE_KJ_MOL = 1e-7
TRANSLATION_FORCE_TOLERANCE_KJ_MOL_ANGSTROM = 1e-7
AVOGADRO_MOL_INVERSE = 6.02214076e23
PARITY_ARTIFACTS = ("input.lammps", "system.data", "log.lammps", "forces.dump")
PROTOCOL_ARTIFACTS = ("input.lammps", "system.data", "log.lammps", "final.dump", "prepared.data")
SOURCE_MODULES = (
    "polychemmd", "polychemmd.geometry", "polychemmd.provenance",
    "polychemmd.forcefields.alkane_data", "polychemmd.forcefields.alkanes",
    "polychemmd.forcefields.contracts", "polychemmd.forcefields.reference",
    "polychemmd.forcefields.periodic", "polychemmd.forcefields.periodic_reference",
    "polychemmd.descriptors.torsion", "polychemmd.engines.lammps",
    "polychemmd.engines.alkane", "polychemmd.engines.periodic_alkane",
    "polychemmd.protocols.alkane",
)


class PeriodicAuditError(RuntimeError):
    """An actual scientific input, output or declared audit gate failed."""


def require(condition, message):
    # Checks must remain active with python -O.
    if not condition:
        raise PeriodicAuditError(message)


def canonical_bytes(record):
    return json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def digest(raw):
    return hashlib.sha256(raw).hexdigest()


def write_json(path, record):
    path.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def source_records():
    records = {}
    for name in SOURCE_MODULES:
        actual = Path(importlib.import_module(name).__file__).resolve()
        relative = Path("src") / Path(*name.split("."))
        expected = REPOSITORY / (relative / "__init__.py" if name == "polychemmd"
                                 else relative.with_suffix(".py"))
        require(actual.is_file() and expected.is_file(), "Source checkout required for this auditor")
        actual_hash = digest(actual.read_bytes())
        require(actual_hash == digest(expected.read_bytes()),
                f"Imported {name} differs from the checkout; install this checkout in editable mode")
        records[expected.relative_to(REPOSITORY).as_posix()] = {"sha256": actual_hash}
    records["scripts/verify_periodic_v05.py"] = {"sha256": digest(Path(__file__).read_bytes())}
    return records


def load_recipe(raw):
    def reject(value):
        raise ValueError("Recipe contains nonfinite JSON")
    def unique(pairs):
        result = {}
        for key, value in pairs:
            require(key not in result, "Recipe contains duplicate fields")
            result[key] = value
        return result
    recipe = json.loads(raw, parse_constant=reject, object_pairs_hook=unique)
    require(recipe["schema_version"] == "1" and recipe["model_id"] == MODEL_ID,
            "Unsupported periodic validation recipe/model")
    path = Path(recipe["source_cell_path"])
    require(not path.is_absolute() and ".." not in path.parts,
            "Source cell path must stay within this checkout")
    for name in ("cutoff_angstrom", "ewald_alpha_per_angstrom", "energy_abs_tol_kj_mol",
                 "force_abs_tol_kj_mol_angstrom", "timeout_seconds_per_engine"):
        value = recipe[name]
        require(type(value) in (int, float) and math.isfinite(float(value)) and value > 0,
                f"{name} must be positive and finite")
    require(recipe["timeout_seconds_per_engine"] <= 120, "Engine timeout exceeds auditor budget")
    require(type(recipe["reciprocal_nmax"]) is int and 1 <= recipe["reciprocal_nmax"] <= 20,
            "Reciprocal nmax is outside this bounded auditor")
    translation = recipe["boundary_translation_angstrom"]
    require(isinstance(translation, list) and len(translation) == 3
            and all(type(v) in (int, float) and math.isfinite(float(v)) for v in translation),
            "Boundary translation must contain three finite numbers")
    # This also validates the pilot's complete recorded resource/model request.
    spec = PeriodicAlkaneProtocolSpec(**recipe["protocol"])
    for name in ("cutoff_angstrom", "ewald_alpha_per_angstrom", "reciprocal_nmax"):
        require(getattr(spec, name) == recipe[name], "Parity and preparation model settings differ")
    return recipe, spec


def translated_cell(source, translation):
    """Translate a physical periodic cell, choosing a whole-molecule image.

    Source geometries and rotations stay fixed; centroids are shifted modulo
    the box. The result is the same physical cell, not an independent replicate.
    """
    cell = deepcopy(source)
    lengths = cell["box"]["lengths_angstrom"]
    molecule_deltas = {}
    for molecule in cell["molecules"]:
        index = molecule["molecule_index"]
        old = molecule["placement"]["translated_centroid_angstrom"]
        new = [(v+d) % length for v, d, length in zip(old, translation, lengths)]
        delta = [a-b for a, b in zip(new, old)]
        molecule_deltas[index] = delta
        molecule["placement"]["translated_centroid_angstrom"] = new
    image_by_id = {}
    for atom in cell["atoms"]:
        point = [v+d for v, d in zip(atom["coordinates_angstrom"],
                                     molecule_deltas[atom["molecule_index"]])]
        images = [math.floor(value/length) for value, length in zip(point, lengths)]
        wrapped = [value-image*length for value, image, length in zip(point, images, lengths)]
        require(all(0 <= value < length for value, length in zip(wrapped, lengths)),
                "Translated wrapped coordinate is outside the half-open cell")
        atom["coordinates_angstrom"] = point
        atom["wrapped_coordinates_angstrom"] = wrapped
        atom["image_flags"] = images
        image_by_id[atom["cell_atom_id"]] = images
    crossings = sum(image_by_id[bond["cell_atom_ids"][0]] != image_by_id[bond["cell_atom_ids"][1]]
                    for bond in cell["bonds"])
    require(crossings > 0, "Declared translated fixture did not put a bond across a face")
    transform = {
        "method": "global translation with whole-molecule modulo-centroid image selection",
        "source_cell_canonical_sha256": digest(canonical_bytes(source)),
        "global_translation_angstrom": list(translation),
        "whole_molecule_translations_angstrom": {str(k): v for k, v in molecule_deltas.items()},
        "physical_interpretation": "same periodic configuration; not an independent sampling replicate",
        "bonds_crossing_cell_faces": crossings,
        "software": {"polychemmd": __version__, "python": platform.python_version()},
    }
    cell.setdefault("preparation", {})["audit_periodic_translation"] = transform
    return cell, transform


def audit_assignment(record, source):
    require(record["model_id"] == MODEL_ID and record["source_cell"] == source,
            "Assigned model/source cell changed")
    atoms = record["atoms"]
    ids = [row["cell_atom_id"] for row in atoms]
    require(len(ids) == source["atom_count"] and len(set(ids)) == len(ids), "Global atom identities differ")
    require([row["atom_map"] for row in atoms] == list(range(1, len(atoms)+1)),
            "Global atom maps must be sequential and unique in this fixture")
    for original, current in zip(source["atoms"], atoms):
        require(all(current[name] == value for name, value in original.items()),
                "Assignment changed a source atom/coordinate/lineage field")
    by_id = {row["cell_atom_id"]: row for row in atoms}
    adjacency = {identity: set() for identity in ids}
    for bond in record["interactions"]["bonds"]:
        first, second = bond["atom_ids"]
        require(first in by_id and second in by_id and first != second, "Invalid global bond endpoints")
        require(by_id[first]["molecule_index"] == by_id[second]["molecule_index"],
                "A covalent bond incorrectly joins separate molecules")
        adjacency[first].add(second)
        adjacency[second].add(first)
    expected = Counter()
    type_counts = Counter(row["type_id"] for row in atoms)
    for molecule in record["molecules"]:
        index = molecule["molecule_index"]
        current = [row for row in atoms if row["molecule_index"] == index]
        carbons = sum(row["element"] == "C" for row in current)
        require(len(current) == 3*carbons+2, "Finite n-alkane atomic balance failed")
        expected.update({"bonds": 3*carbons+1, "angles": 6*carbons,
                         "dihedrals": 9*(carbons-1), "impropers": 0})
        expected_types = Counter({"opls_135": 2, "opls_136": carbons-2, "opls_140": 2*carbons+2})
        require(Counter(row["type_id"] for row in current) == +expected_types,
                "CH3/CH2/H type populations disagree with the molecule")
        for atom in current:
            carbon_neighbours = sum(by_id[other]["element"] == "C"
                                    for other in adjacency[atom["cell_atom_id"]])
            expected_type = ("opls_140" if atom["element"] == "H" else
                             "opls_135" if carbon_neighbours == 1 else "opls_136")
            expected_charge = {"opls_135": -.18, "opls_136": -.12, "opls_140": .06}[expected_type]
            require(atom["type_id"] == expected_type and atom["charge_e"] == expected_charge,
                    "An atom's chemical environment/type/fixed charge assignment changed")
        require(abs(math.fsum(row["charge_e"] for row in current)) <= 1e-12,
                "An assigned molecule is nonneutral")
    actual = {name: len(rows) for name, rows in record["interactions"].items()}
    require(actual == dict(expected), "Global interaction counts disagree with finite n-alkanes")
    lengths = record["box"]["lengths_angstrom"]
    reconstruction_error = 0.0
    for atom in atoms:
        for value, wrapped, image, length in zip(atom["coordinates_angstrom"],
                                               atom["wrapped_coordinates_angstrom"],
                                               atom["image_flags"], lengths):
            require(type(image) is int and 0 <= wrapped < length, "Invalid image/wrapped coordinate")
            error = abs(value-(wrapped+image*length))
            reconstruction_error = max(reconstruction_error, error)
            require(error <= 1e-8, "Unwrapped/image reconstruction failed")
    pairs = record["nonbonded_pairs"]
    require(len(pairs) == len(ids)*(len(ids)-1)//2, "Missing or extra unordered global pairs")
    pair_by_ids = {frozenset(pair["atom_ids"]): pair for pair in pairs}
    require(len(pair_by_ids) == len(pairs), "Duplicated global pair")
    pair_counts = Counter()
    for i, first in enumerate(ids):
        distances, queue = {first: 0}, deque([first])
        while queue:
            current = queue.popleft()
            for neighbour in adjacency[current]:
                if neighbour not in distances:
                    distances[neighbour] = distances[current]+1
                    queue.append(neighbour)
        for second in ids[i+1:]:
            pair = pair_by_ids[frozenset((first, second))]
            same = by_id[first]["molecule_index"] == by_id[second]["molecule_index"]
            distance = distances[second] if same else None
            scale = (0.0 if distance <= 2 else .5 if distance == 3 else 1.0) if same else 1.0
            require(pair["same_molecule"] is same and pair["topological_distance"] == distance,
                    "Pair molecule/shortest-path identity failed")
            require(pair["lj_scale"] == scale and pair["coulomb_scale"] == scale,
                    "Pair exclusion or scaled 1-4 treatment changed")
            pair_counts["intermolecular" if not same else "excluded" if scale == 0 else
                        "scaled_14" if scale == .5 else "intramolecular_unscaled"] += 1
    charge = math.fsum(row["charge_e"] for row in atoms)
    require(abs(charge) <= 1e-12, "The assigned periodic cell is nonneutral")
    return {"status": "passed", "atom_count": len(atoms), "molecule_count": len(record["molecules"]),
            "interaction_counts": actual, "atom_type_counts": dict(type_counts),
            "pair_counts": dict(pair_counts), "total_partial_charge_e": charge,
            "maximum_image_reconstruction_error_angstrom": reconstruction_error,
            "independent_pair_distance_method": "breadth-first search of global covalent graph"}


def balance_checks(energies, forces):
    error = abs(math.fsum(energies[name] for name in COMPONENTS)-energies["total"])
    require(error <= ENERGY_SUM_TOLERANCE_KJ_MOL, "Energy components do not sum to total")
    require(all(math.isfinite(value) for row in forces for value in row), "Nonfinite force output")
    net = [math.fsum(row[axis] for row in forces) for axis in range(3)]
    require(math.hypot(*net) <= NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM, "Periodic total force does not balance")
    return {"energy_component_sum_error_kj_mol": error, "net_force_kj_mol_angstrom": net,
            "net_force_norm_kj_mol_angstrom": math.hypot(*net),
            "net_torque_test": "not_applied; a fixed periodic lattice does not have continuous rotation invariance"}


def public_record(raw):
    record = deepcopy(raw)
    changes = []
    def visit(value, location="record"):
        if isinstance(value, list):
            for i, item in enumerate(value):
                visit(item, f"{location}[{i}]")
        elif isinstance(value, dict):
            if "run_directory" in value:
                value["run_directory"] = "private_engine_run_not_published"
                changes.append(location+".run_directory replaced with a non-host label")
            if isinstance(value.get("executable"), str):
                value["executable"] = Path(value["executable"]).name
                changes.append(location+".executable reduced to basename; identity hash retained")
            if isinstance(value.get("command"), list) and value["command"]:
                value["command"][0] = Path(value["command"][0]).name
            first_line = value.get("engine_help_first_line")
            if isinstance(first_line, str) and any(word in first_line for word in ("UCX", "hwloc")):
                value["engine_help_first_line"] = "Incidental host diagnostic omitted; see engine_version_banner"
                changes.append(location+".engine_help_first_line omitted incidental host diagnostic")
            for key, item in value.items():
                visit(item, location+"."+key)
    visit(record)
    record["public_record_sanitization"] = {
        "changes": changes,
        "scope": "Only host paths/incidental help diagnostics replaced; scientific values, engine identity, settings and raw artifact hashes are retained.",
    }
    return record


def copy_artifacts(run, run_directory, target, names):
    for name in names:
        source = run_directory/name
        require(source.is_file() and source.stat().st_size > 0, f"Actual scientific artifact missing: {name}")
        shutil.copyfile(source, target/name)
    require(digest((target/"input.lammps").read_bytes()) == run["input_sha256"], "Input script hash changed")
    for name, expected in run.get("artifact_sha256", {}).items():
        if name in names:
            require(digest((target/name).read_bytes()) == expected, f"Copied artifact hash changed: {name}")
    for name, expected in run.get("staged_input_sha256", {}).items():
        if name in names:
            require(digest((target/name).read_bytes()) == expected, f"Copied staged input changed: {name}")


def artifact_records(directory):
    return {path.name: {"sha256": digest(path.read_bytes()), "bytes": path.stat().st_size}
            for path in sorted(directory.iterdir()) if path.is_file()}


def assignment_summary(record, audit):
    return {"kind": "periodic_assignment_summary", "model_id": record["model_id"],
            "full_assignment_canonical_sha256": digest(canonical_bytes(record)),
            "full_assignment_location": "private_validation_inputs; regenerate from public cell and frozen source model",
            "parameter_set_canonical_sha256": digest(canonical_bytes(record["parameter_set"])),
            "source_cell_canonical_sha256": digest(canonical_bytes(record["source_cell"])),
            "source_snapshot": record["source_snapshot"], "box": record["box"],
            "independent_topology_audit": audit,
            "atom_identity_map": [{"engine_id": i+1, **{name: atom[name] for name in
                                  ("cell_atom_id", "atom_map", "source_atom_id", "source_atom_map",
                                   "molecule_index", "type_id", "charge_e")}}
                                  for i, atom in enumerate(record["atoms"])],
            "validation": record["validation"]}


def translation_invariance(original, moved):
    checks = {}
    for name, energy_key, force_key in (
            ("reference", None, None),
            ("engine", "engine_energy_components_kj_mol", "engine_forces_kj_mol_angstrom")):
        first = original["reference"] if name == "reference" else original
        second = moved["reference"] if name == "reference" else moved
        energies_a = first["energy_components_kj_mol"] if name == "reference" else first[energy_key]
        energies_b = second["energy_components_kj_mol"] if name == "reference" else second[energy_key]
        forces_a = first["forces_kj_mol_angstrom"] if name == "reference" else first[force_key]
        forces_b = second["forces_kj_mol_angstrom"] if name == "reference" else second[force_key]
        errors = {key: abs(energies_a[key]-energies_b[key]) for key in (*COMPONENTS, "total")}
        maximum_force = max(abs(a-b) for left, right in zip(forces_a, forces_b) for a, b in zip(left, right))
        require(max(errors.values()) <= TRANSLATION_ENERGY_TOLERANCE_KJ_MOL,
                f"{name} global periodic-translation energy invariance failed")
        require(maximum_force <= TRANSLATION_FORCE_TOLERANCE_KJ_MOL_ANGSTROM,
                f"{name} global periodic-translation force invariance failed")
        checks[name] = {"energy_component_abs_errors_kj_mol": errors,
                        "max_force_component_abs_error_kj_mol_angstrom": maximum_force}
    return {"status": "passed", "interpretation": "two representations of the same physical cell, not independent samples",
            "predeclared_energy_abs_tolerance_kj_mol": TRANSLATION_ENERGY_TOLERANCE_KJ_MOL,
            "predeclared_force_abs_tolerance_kj_mol_angstrom": TRANSLATION_FORCE_TOLERANCE_KJ_MOL_ANGSTROM,
            "checks": checks}


def audit_final_state(protocol, assignment, spec):
    require(protocol["equilibrium_status"] in ("insufficient_sampling", "not_assessed"),
            "Short preparation pilot cannot certify equilibrium")
    require(protocol["production_ready"] is False and protocol["density_prediction_validation"] == "not_performed",
            "Preparation pilot must not claim production/material validation")
    stages = protocol["stages"]
    require(stages["nvt"]["status"] == stages["npt"]["status"] == "executed", "Requested MD pilot stages did not execute")
    require(stages["nvt"]["steps"] == spec.nvt_steps and stages["npt"]["steps"] == spec.npt_steps,
            "Actual pilot step counts differ from declared resource request")
    state = protocol["final_state"]
    require(state["atom_count"] == len(assignment["atoms"]) == len(state["atoms"]), "Pilot lost atoms")
    by_id = {atom["cell_atom_id"]: atom for atom in assignment["atoms"]}
    require({atom["cell_atom_id"] for atom in state["atoms"]} == set(by_id), "Pilot atom identity changed")
    lengths = state["box"]["lengths_angstrom"]
    origin = state["box"]["origin_angstrom"]
    require(all(length > 2*spec.cutoff_angstrom for length in lengths), "Final box violates half-box cutoff gate")
    for atom in state["atoms"]:
        source = by_id[atom["cell_atom_id"]]
        require(atom["molecule_index"] == source["molecule_index"] and atom["charge_e"] == source["charge_e"],
                "Pilot molecule/charge identity changed")
        for value, wrapped, image, length, lower in zip(atom["coordinates_angstrom"],
                                               atom["wrapped_coordinates_angstrom"], atom["image_flags"], lengths, origin):
            require(type(image) is int and lower <= wrapped < lower+length
                    and abs(value-(wrapped+image*length)) <= 1e-8, "Pilot final periodic image identity failed")
    volume = math.prod(lengths)
    mass = math.fsum(atom["mass_g_mol"] for atom in assignment["atoms"])
    density = mass/(AVOGADRO_MOL_INVERSE*volume*1e-24)
    require(math.isclose(volume, state["volume_angstrom3"], rel_tol=1e-9, abs_tol=1e-8), "Final volume mismatch")
    require(math.isclose(density, state["density_g_cm3"], rel_tol=1e-12, abs_tol=1e-12),
            "Final canonical density disagrees with exact-SI mass/volume arithmetic")
    require(math.isfinite(state["engine_density_g_cm3"]) and state["engine_density_g_cm3"] > 0,
            "Actual native engine density must be retained separately")
    return {"status": "passed", "atom_identity_preserved": True,
            "requested_md_steps": spec.nvt_steps+spec.npt_steps,
            "final_volume_angstrom3": volume, "final_density_mass_volume_g_cm3": density,
            "raw_engine_density_g_cm3": state["engine_density_g_cm3"],
            "density_conventions": "Canonical exact-SI mass/volume and actual engine-native thermo density are distinct recorded fields",
            "density_interpretation": "actual short-pilot state arithmetic; not an equilibrated prediction",
            "equilibrium_status": protocol["equilibrium_status"],
            "minimization_status": stages["minimization"]["status"]}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/periodic-v05")
    parser.add_argument("--executable")
    parser.add_argument("--skip-protocol", action="store_true")
    args = parser.parse_args(argv)
    recipe_raw = RECIPE_PATH.read_bytes()
    recipe, spec = load_recipe(recipe_raw)
    modules = source_records()
    source_path = REPOSITORY/recipe["source_cell_path"]
    source_raw = source_path.read_bytes()
    source_record = json.loads(source_raw)
    source = source_record["payload"] if "payload" in source_record else source_record
    translated, transformation = translated_cell(source, recipe["boundary_translation_angstrom"])
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=False)
    private_runs, private_inputs = output/"private_engine_runs", output/"private_validation_inputs"
    private_runs.mkdir(); private_inputs.mkdir()
    os.environ["OMP_NUM_THREADS"], os.environ["UCX_TLS"] = "1", "self"
    started_utc, started = datetime.now(timezone.utc).isoformat(), time.monotonic()
    summaries, actual_parities, assignments = [], [], []
    for name, cell in (("original", source), ("boundary_translated", translated)):
        directory = output/name; directory.mkdir()
        write_json(directory/"cell.json", envelope(cell, source_raw, "periodic_validation_source_cell"))
        assigned = assign_periodic_alkanes(cell)
        assignment = assigned.to_dict()
        assignments.append(assigned)
        write_json(private_inputs/(name+"-assignment.json"), assignment)
        topology = audit_assignment(assignment, cell)
        write_json(directory/"assignment_summary.json", assignment_summary(assignment, topology))
        run_dir = private_runs/name
        parity = verify_periodic_alkane_singlepoint(
            assigned, run_directory=run_dir, executable=args.executable,
            timeout_seconds=recipe["timeout_seconds_per_engine"],
            cutoff_angstrom=recipe["cutoff_angstrom"],
            ewald_alpha_per_angstrom=recipe["ewald_alpha_per_angstrom"],
            reciprocal_nmax=recipe["reciprocal_nmax"],
            energy_abs_tol_kj_mol=recipe["energy_abs_tol_kj_mol"],
            force_abs_tol_kj_mol_angstrom=recipe["force_abs_tol_kj_mol_angstrom"])
        require(parity["status"] == "passed", "Actual periodic energy/force parity failed")
        checks = {"reference": balance_checks(parity["reference"]["energy_components_kj_mol"],
                                               parity["reference"]["forces_kj_mol_angstrom"]),
                  "engine": balance_checks(parity["engine_energy_components_kj_mol"],
                                            parity["engine_forces_kj_mol_angstrom"])}
        public = public_record(parity); public["independent_balance_checks"] = checks
        write_json(directory/"parity.json", public)
        copy_artifacts(parity["engine_run"], run_dir, directory, PARITY_ARTIFACTS)
        actual_parities.append(parity)
        summaries.append({"case_name": name, "status": "passed", "topology_audit": topology,
                          "cell_canonical_sha256": digest(canonical_bytes(cell)),
                          "assignment_canonical_sha256": digest(canonical_bytes(assignment)),
                          "energy_component_abs_errors_kj_mol": parity["energy_component_abs_errors_kj_mol"],
                          "maximum_energy_component_abs_error_kj_mol": max(parity["energy_component_abs_errors_kj_mol"].values()),
                          "max_force_component_abs_error_kj_mol_angstrom": parity["max_force_component_abs_error_kj_mol_angstrom"],
                          "engine_version_banner": parity["engine_run"]["engine_version_banner"],
                          "engine_executable_sha256": parity["engine_run"]["executable_sha256"],
                          "independent_balance_checks": checks, "artifact_files": artifact_records(directory)})
    invariance = translation_invariance(*actual_parities)
    protocol_summary = None
    if not args.skip_protocol:
        directory = output/"preparation_pilot"; directory.mkdir()
        run_dir = private_runs/"preparation_pilot"
        protocol = run_periodic_alkane_protocol(assignments[0], spec,
                                               run_directory=run_dir, executable=args.executable)
        final_audit = audit_final_state(protocol, assignments[0].to_dict(), spec)
        public = public_record(protocol); public["independent_final_state_audit"] = final_audit
        write_json(directory/"protocol.json", public)
        copy_artifacts(protocol["engine_run"], run_dir, directory, PROTOCOL_ARTIFACTS)
        protocol_summary = {"status": protocol["status"], "specification": spec.to_dict(),
                            "stages": protocol["stages"], "independent_final_state_audit": final_audit,
                            "equilibrium_status": protocol["equilibrium_status"],
                            "density_prediction_validation": protocol["density_prediction_validation"],
                            "artifact_files": artifact_records(directory)}
    parameters = load_alkane_parameters()
    report = {
        "schema_version": "1", "evidence_kind": "executed_periodic_ewald_parity_and_bounded_preparation_audit",
        "status": "passed", "started_utc": started_utc, "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic()-started,
        "command_template": "python scripts/verify_periodic_v05.py --output-dir NEW_DIRECTORY [--executable LAMMPS] [--skip-protocol]",
        "software": {"polychemmd": __version__, "python": platform.python_version(),
                     "python_implementation": platform.python_implementation()},
        "runtime_environment": {"OMP_NUM_THREADS": os.environ["OMP_NUM_THREADS"], "UCX_TLS": os.environ["UCX_TLS"]},
        "input_recipe": {"path": RECIPE_PATH.relative_to(REPOSITORY).as_posix(), "sha256": digest(recipe_raw), "record": recipe},
        "original_cell_input": {"path": recipe["source_cell_path"], "sha256": digest(source_raw)},
        "source_module_hashes": modules,
        "parameter_json_sha256": digest(files("polychemmd.forcefields").joinpath("alkane_parameters.json").read_bytes()),
        "source_snapshot": parameters["source_snapshot"],
        "translation_case": transformation, "translation_invariance": invariance,
        "predeclared_independent_tolerances": {"component_sum_kj_mol": ENERGY_SUM_TOLERANCE_KJ_MOL,
                                               "net_force_norm_kj_mol_angstrom": NET_FORCE_TOLERANCE_KJ_MOL_ANGSTROM,
                                               "translation_energy_abs_kj_mol": TRANSLATION_ENERGY_TOLERANCE_KJ_MOL,
                                               "translation_force_component_abs_kj_mol_angstrom": TRANSLATION_FORCE_TOLERANCE_KJ_MOL_ANGSTROM},
        "cases": summaries, "bounded_protocol": protocol_summary,
        "maximum_energy_component_abs_error_kj_mol": max(c["maximum_energy_component_abs_error_kj_mol"] for c in summaries),
        "max_force_component_abs_error_kj_mol_angstrom": max(c["max_force_component_abs_error_kj_mol_angstrom"] for c in summaries),
        "protocol_execution_skipped": args.skip_protocol,
        "md_integration_steps": 0 if args.skip_protocol else spec.nvt_steps+spec.npt_steps,
        "equilibrium_status": "not_assessed" if args.skip_protocol else "insufficient_sampling",
        "production_ready": False, "material_property_validation": "not_performed",
        "public_record_policy": "Publish only case directories, preparation_pilot if present and report.json. private_engine_runs and private_validation_inputs are excluded; full assignments can be regenerated from public cells and frozen model. Scientific engine artifacts are byte-preserved; public reports replace host paths/incidental help diagnostics only.",
        "scientific_limits": [
            "The two parity cases are lattice-equivalent representations of the same periodic cell, not independent samples.",
            "Periodic Ewald numerical parity and bounded process execution do not establish physical force-field accuracy.",
            "The selected finite real-space cutoff and reciprocal representation are explicit numerical model settings, not an independently converged bulk property claim.",
            "The short minimization/NVT/NPT pilot, if executed, is preparation evidence with insufficient sampling; its density is state arithmetic, not a validated equilibrium prediction.",
            "No new force-field fit, quantum calculation, production trajectory or experimental evaluation was performed.",
        ],
    }
    write_json(output/"report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output/"report.json"),
                      "maximum_energy_component_abs_error_kj_mol": report["maximum_energy_component_abs_error_kj_mol"],
                      "max_force_component_abs_error_kj_mol_angstrom": report["max_force_component_abs_error_kj_mol_angstrom"],
                      "md_integration_steps": report["md_integration_steps"],
                      "equilibrium_status": report["equilibrium_status"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
