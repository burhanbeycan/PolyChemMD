#!/usr/bin/env python3
"""Execute and independently audit the small v0.3 coordinate examples.

Run from a prepared PolyChemMD checkout after installing ``.[chemistry]``:
``python scripts/verify_geometry_v03.py --output-dir artifacts/geometry-v03``.
The output directory must be new. This bounded geometry check runs no MD/QM,
fits no parameters and does not validate an equilibrium polymer property.
"""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
from itertools import combinations, product
import json
import math
from pathlib import Path
import platform
import subprocess
import sys
import time

from polychemmd import __version__
from polychemmd.builders.conformers import generate_conformer
from polychemmd.builders.linear import build_linear
from polychemmd.chemistry import RepeatUnitSpec
from polychemmd.geometry import MolecularGeometry
from polychemmd.provenance import envelope


REPOSITORY = Path(__file__).resolve().parents[1]
ABS_TOLERANCE_ANGSTROM = 1e-9
REL_TOLERANCE = 1e-10
MATRIX_TOLERANCE = 1e-12
AVOGADRO_MOL_INVERSE = 6.02214076e23
NEIGHBOUR_IMAGES = tuple(product((-1, 0, 1), repeat=3))


class GeometryAuditError(RuntimeError):
    """An independently checked output did not satisfy its stated contract."""


def require(condition, message):
    # Use explicit checks so ``python -O`` cannot disable the audit.
    if not condition:
        raise GeometryAuditError(message)


def close(first, second, *, absolute=ABS_TOLERANCE_ANGSTROM, relative=REL_TOLERANCE):
    return math.isfinite(first) and math.isfinite(second) and math.isclose(
        first, second, abs_tol=absolute, rel_tol=relative)


def canonical_bytes(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"),
                      allow_nan=False).encode("utf-8")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def chain_from_mapping(data):
    return build_linear(RepeatUnitSpec.from_mapping(data["repeat_unit"]),
                        data["degree_of_polymerization"],
                        head_cap=data["head_cap"], tail_cap=data["tail_cap"])


def run_cli(command, input_path, output_path, *, timeout_seconds):
    arguments = [sys.executable, "-m", "polychemmd", command, str(input_path),
                 "--output", str(output_path)]
    started = time.monotonic()
    completed = subprocess.run(arguments, cwd=REPOSITORY, check=False,
                               capture_output=True, text=True,
                               timeout=timeout_seconds)
    require(completed.returncode == 0,
            f"{command} failed with status {completed.returncode}: {completed.stderr}")
    require(output_path.is_file() and output_path.stat().st_size > 0,
            f"{command} did not produce a fresh nonempty output")
    return {"command": command, "returncode": completed.returncode,
            "timeout_seconds": timeout_seconds,
            "elapsed_seconds": time.monotonic() - started,
            "output_file": output_path.name}


def determinant(matrix):
    a, b, c = matrix
    return (a[0] * (b[1] * c[2] - b[2] * c[1])
            - a[1] * (b[0] * c[2] - b[2] * c[0])
            + a[2] * (b[0] * c[1] - b[1] * c[0]))


def independently_audit_cell(cell):
    """Audit serialized coordinates without using the packer's distance helpers."""
    require(cell["kind"] == "orthorhombic_rigid_packed_cell", "unexpected cell kind")
    require(cell["coordinate_units"] == "angstrom", "unexpected coordinate units")
    require(cell["box"]["periodic"] == [True, True, True], "all axes must be periodic")
    require(cell["box"]["origin_angstrom"] == [0.0, 0.0, 0.0], "unexpected box origin")
    require(cell["simulation_performed"] is False, "geometric preparation cannot claim MD")
    require(cell["equilibrium_status"] == "not_assessed", "unexpected equilibrium claim")
    require(cell["phase"] == "not_assigned", "unexpected phase claim")
    require(cell["material_validation"] == "not_performed", "unexpected material validation")
    require(cell["validation"]["force_field_assigned"] is False, "unexpected force-field claim")

    box = tuple(cell["box"]["lengths_angstrom"])
    require(len(box) == 3 and all(math.isfinite(v) and v > 0 for v in box),
            "invalid box dimensions")
    cutoff = cell["packing"]["min_distance_angstrom"]
    require(math.isfinite(cutoff) and cutoff > 0, "invalid geometric cutoff")
    atoms, bonds, molecules = cell["atoms"], cell["bonds"], cell["molecules"]
    require(len(atoms) == cell["atom_count"], "atom count mismatch")
    require(len(bonds) == cell["bond_count"], "bond count mismatch")
    require(len(molecules) == cell["molecule_count"], "molecule count mismatch")
    require(len({a["cell_atom_id"] for a in atoms}) == len(atoms), "duplicate global atom IDs")
    require([a["atom_map"] for a in atoms] == list(range(1, len(atoms) + 1)),
            "global maps must be unique and sequential")

    grouped, source_graph_hashes = [], []
    maximum_pair_error = 0.0
    maximum_coordinate_reconstruction_error = 0.0
    maximum_placement_transform_error = 0.0
    minimum_self_image_distance = None
    intramolecular_pairs = self_image_pairs = 0
    maximum_matrix_error = 0.0
    for index, molecule in enumerate(molecules):
        require(molecule["molecule_index"] == index, "molecule index mismatch")
        source = MolecularGeometry.from_mapping(molecule["source_geometry"])
        source_graph_hashes.append(source.source_graph_sha256)
        require(molecule["source_graph_sha256"] == source.source_graph_sha256,
                "source graph hash mismatch")
        current = [a for a in atoms if a["molecule_index"] == index]
        require(len(current) == len(source.atoms), "source/current molecule sizes differ")
        source_by_id = {a.atom_id: a for a in source.atoms}
        source_positions = dict(zip(source_by_id, source.coordinates_angstrom))
        current_by_id = {a["source_atom_id"]: a for a in current}
        require(set(current_by_id) == set(source_by_id), "source atom IDs were lost or duplicated")
        ordered = [current_by_id[a.atom_id] for a in source.atoms]

        matrix = molecule["placement"]["rotation_matrix"]
        require(len(matrix) == 3 and all(len(row) == 3 for row in matrix),
                "invalid rotation matrix dimensions")
        require(close(determinant(matrix), 1.0, absolute=MATRIX_TOLERANCE,
                      relative=MATRIX_TOLERANCE), "placement is not a proper rotation")
        for i, row in enumerate(matrix):
            for j, other in enumerate(matrix):
                error = abs(math.fsum(a * b for a, b in zip(row, other)) - float(i == j))
                maximum_matrix_error = max(maximum_matrix_error, error)
                require(error <= MATRIX_TOLERANCE, "rotation rows are not orthonormal")

        centroid = molecule["placement"]["source_centroid_angstrom"]
        translation = molecule["placement"]["translated_centroid_angstrom"]
        for record in ordered:
            original = source_by_id[record["source_atom_id"]]
            require(record["cell_atom_id"] == f"molecule:{index}/{original.atom_id}",
                    "cell/source atom ID relationship changed")
            require(record["source_atom_map"] == original.atom_map, "source atom map changed")
            for name in ("element", "isotope", "formal_charge", "mass_g_mol"):
                require(record[name] == getattr(original, name), f"atomic {name} changed")
            require(record["lineage"] == original.lineage, "source atomic lineage changed")
            position = record["coordinates_angstrom"]
            wrapped = record["wrapped_coordinates_angstrom"]
            images = record["image_flags"]
            require(len(position) == len(wrapped) == len(images) == 3,
                    "coordinate or image vector has wrong size")
            require(all(type(image) is int for image in images), "image flags must be integers")
            for x, w, image, length in zip(position, wrapped, images, box):
                require(math.isfinite(x) and math.isfinite(w) and 0 <= w < length,
                        "wrapped coordinate is outside the half-open box")
                error = abs(x - (w + image * length))
                maximum_coordinate_reconstruction_error = max(
                    maximum_coordinate_reconstruction_error, error)
                require(close(x, w + image * length), "unwrapped/wrapped/image identity failed")
            centred = tuple(x - c for x, c in zip(source_positions[original.atom_id], centroid))
            predicted = tuple(math.fsum(a * b for a, b in zip(row, centred)) + offset
                              for row, offset in zip(matrix, translation))
            for actual, expected in zip(position, predicted):
                maximum_placement_transform_error = max(
                    maximum_placement_transform_error, abs(actual - expected))
                require(close(actual, expected), "stored rigid transform does not reproduce coordinates")

        for first, second in combinations(range(len(ordered)), 2):
            original = math.dist(source.coordinates_angstrom[first], source.coordinates_angstrom[second])
            moved = math.dist(ordered[first]["coordinates_angstrom"], ordered[second]["coordinates_angstrom"])
            maximum_pair_error = max(maximum_pair_error, abs(original - moved))
            require(close(original, moved), "an intramolecular pair distance changed")
            intramolecular_pairs += 1

        for axis, length in enumerate(box):
            span = max(a["coordinates_angstrom"][axis] for a in ordered) - min(
                a["coordinates_angstrom"][axis] for a in ordered)
            require(span <= length - cutoff + ABS_TOLERANCE_ANGSTROM,
                    "self-image span bound failed")
        for image in NEIGHBOUR_IMAGES:
            if image == (0, 0, 0):
                continue
            for first in ordered:
                for second in ordered:
                    shifted = tuple(v + n * length for v, n, length in zip(
                        second["coordinates_angstrom"], image, box))
                    distance = math.dist(first["coordinates_angstrom"], shifted)
                    require(distance >= cutoff - ABS_TOLERANCE_ANGSTROM,
                            "a molecule contacts a nonzero periodic self image")
                    minimum_self_image_distance = (distance if minimum_self_image_distance is None
                                                   else min(minimum_self_image_distance, distance))
                    self_image_pairs += 1

        current_bonds = [b for b in bonds if b["molecule_index"] == index]
        require(len(current_bonds) == len(source.bonds), "source/current bond counts differ")
        bond_by_pair = {frozenset(b["cell_atom_ids"]): b for b in current_bonds}
        require(len(bond_by_pair) == len(current_bonds), "duplicate global bond pairs")
        for original in source.bonds:
            expected_ids = [f"molecule:{index}/{atom_id}" for atom_id in original.atom_ids]
            current_bond = bond_by_pair.get(frozenset(expected_ids))
            require(current_bond is not None, "source bond was lost")
            expected_maps = {current_by_id[atom_id]["atom_map"] for atom_id in original.atom_ids}
            require(set(current_bond["atom_maps"]) == expected_maps, "bond atom maps changed")
            for name in ("order", "aromatic", "stereo"):
                require(current_bond[name] == getattr(original, name), f"bond {name} changed")
        grouped.append(ordered)

    require(sum(len(group) for group in grouped) == len(atoms), "unassigned cell atoms")
    minimum_inter = None
    intermolecular_pairs = inter_image_distances = 0
    for left, right in combinations(grouped, 2):
        for first in left:
            for second in right:
                x, y = first["wrapped_coordinates_angstrom"], second["wrapped_coordinates_angstrom"]
                distance = min(math.dist(x, tuple(v + n * length for v, n, length in zip(y, image, box)))
                               for image in NEIGHBOUR_IMAGES)
                require(distance >= cutoff - ABS_TOLERANCE_ANGSTROM,
                        "a distinct-molecule periodic contact violates the cutoff")
                minimum_inter = distance if minimum_inter is None else min(minimum_inter, distance)
                intermolecular_pairs += 1
                inter_image_distances += len(NEIGHBOUR_IMAGES)
    reported = cell["packing"]["minimum_inter_molecular_distance_angstrom"]
    require((reported is None and minimum_inter is None) or
            (reported is not None and minimum_inter is not None and close(reported, minimum_inter)),
            "reported minimum intermolecular separation disagrees with brute-force images")

    mass = math.fsum(a["mass_g_mol"] for a in atoms)
    volume = math.prod(box)
    density = mass / (AVOGADRO_MOL_INVERSE * volume * 1e-24)
    require(close(mass, cell["total_molar_mass_g_mol"]), "total mass mismatch")
    require(close(volume, cell["box"]["volume_angstrom3"]), "box volume mismatch")
    require(close(density, cell["construction_density_g_cm3"], absolute=1e-14),
            "construction mass/volume density mismatch")
    require(sum(a["formal_charge"] for a in atoms) == cell["formal_charge"], "total charge mismatch")
    return {
        "status": "passed",
        "molecule_count": len(molecules), "atom_count": len(atoms), "bond_count": len(bonds),
        "source_graph_sha256": source_graph_hashes,
        "audit_methods": ["independent 27-image intermolecular enumeration",
                          "all-atom pairs across 26 nonzero adjacent self images plus span bound",
                          "all intramolecular pair distances", "proper orthonormal rotation matrices",
                          "source/global identity, lineage and bonds", "image reconstruction and rigid-transform parity"],
        "intermolecular_atom_pairs": intermolecular_pairs,
        "intermolecular_image_distances_evaluated": inter_image_distances,
        "self_image_distances_evaluated": self_image_pairs,
        "intramolecular_pairs_checked": intramolecular_pairs,
        "minimum_inter_molecular_distance_angstrom": minimum_inter,
        "minimum_nonzero_neighbour_self_image_distance_angstrom": minimum_self_image_distance,
        "maximum_intramolecular_distance_error_angstrom": maximum_pair_error,
        "maximum_image_reconstruction_error_angstrom": maximum_coordinate_reconstruction_error,
        "maximum_rigid_transform_error_angstrom": maximum_placement_transform_error,
        "maximum_rotation_orthogonality_error": maximum_matrix_error,
        "construction_density_g_cm3": density,
        "density_interpretation": "supplied atomic mass divided by prescribed cell volume; not predicted",
        "tolerances_declared_before_audit": {"distance_absolute_angstrom": ABS_TOLERANCE_ANGSTROM,
                                             "relative": REL_TOLERANCE,
                                             "matrix_absolute_and_relative": MATRIX_TOLERANCE},
    }


def run_pilots(output_dir):
    """Small actual coordinate realizations, with no expected/fabricated outputs."""
    cases = (("pe-dp20", "[*:1]CC[*:2]", 20, 11001),
             ("pvdf-dp8", "[*:1]CC(F)(F)[*:2]", 8, 11002),
             ("ps-dp4", "[*:1]CC(c1ccccc1)[*:2]", 4, 11003))
    results = []
    for name, smiles, dp, seed in cases:
        recipe = {"chain": {"repeat_unit": {"label": name + " finite repeat graph",
                                             "smiles": smiles, "head_map": 1, "tail_map": 2},
                            "degree_of_polymerization": dp, "head_cap": "H", "tail_cap": "H"},
                  "seed": seed, "max_iterations": 1000}
        recipe_bytes = canonical_bytes(recipe)
        chain = chain_from_mapping(recipe["chain"])
        started = time.monotonic()
        geometry = generate_conformer(chain, seed=seed, max_iterations=recipe["max_iterations"])
        payload = geometry.to_dict()
        filename = f"pilot-{name}.json"
        write_json(output_dir / filename, envelope(payload, recipe_bytes, "finite_graph_coordinate_pilot"))
        results.append({"name": name, "recipe": recipe, "input_recipe_sha256": digest(recipe_bytes),
                        "input_hash_basis": "canonical JSON of the expert-defined finite repeat-graph recipe",
                        "elapsed_seconds": time.monotonic() - started,
                        "molecular_formula": chain.molecular_formula,
                        "atom_count": len(geometry.atoms), "bond_count": len(geometry.bonds),
                        "formal_charge": sum(a.formal_charge for a in geometry.atoms),
                        "source_graph_sha256": geometry.source_graph_sha256,
                        "conformer_payload_sha256": digest(canonical_bytes(payload)),
                        "coordinate_validation": geometry.provenance["validation"],
                        "artifact_file": filename,
                        "evidence_level": "generated_finite_molecular_geometry",
                        "tacticity": "not_specified_or_inferred",
                        "simulation_performed": False, "material_validation": "not_performed"})
    return results


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="artifacts/geometry-v03")
    parser.add_argument("--skip-pilots", action="store_true",
                        help="Omit the three additional finite-graph coordinate realizations")
    args = parser.parse_args(argv)
    output_dir = Path(args.output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=False)
    started_utc = datetime.now(timezone.utc).isoformat()
    clock = time.monotonic()
    recipes = {name: REPOSITORY / "examples" / name for name in ("conformer_pe.json", "packed_pe.json")}
    steps = [run_cli("build-conformer", recipes["conformer_pe.json"], output_dir / "conformer.json", timeout_seconds=45),
             run_cli("pack-cell", recipes["packed_pe.json"], output_dir / "cell.json", timeout_seconds=180),
             run_cli("export-xyz", output_dir / "cell.json", output_dir / "cell.xyz", timeout_seconds=15)]
    conformer_envelope = json.loads((output_dir / "conformer.json").read_text(encoding="utf-8"))
    cell_envelope = json.loads((output_dir / "cell.json").read_text(encoding="utf-8"))
    conformer, cell = conformer_envelope["payload"], cell_envelope["payload"]
    conformer_geometry = MolecularGeometry.from_mapping(conformer)
    require(conformer_geometry.provenance["simulation_performed"] is False,
            "conformer must not claim a simulation")
    recipe_records = {}
    graphs = {}
    for name, path in recipes.items():
        raw = path.read_bytes()
        source = json.loads(raw)
        graph = chain_from_mapping(source["chain"]).to_dict()
        graphs[name] = digest(canonical_bytes(graph))
        recipe_records[name] = {"path": "examples/" + name, "sha256": digest(raw),
                                "hash_basis": "exact checked-in recipe file bytes",
                                "actual_graph_sha256": graphs[name]}
    require(conformer_envelope["input_sha256"] == recipe_records["conformer_pe.json"]["sha256"],
            "conformer envelope input hash mismatch")
    require(cell_envelope["input_sha256"] == recipe_records["packed_pe.json"]["sha256"],
            "cell envelope input hash mismatch")
    require(conformer["source_graph_sha256"] == graphs["conformer_pe.json"],
            "conformer source graph hash disagrees with actual reconstructed graph")
    audit = independently_audit_cell(cell)
    require(all(h == graphs["packed_pe.json"] for h in audit["source_graph_sha256"]),
            "a cell molecule's source graph differs from the packing recipe")
    xyz_lines = (output_dir / "cell.xyz").read_text(encoding="utf-8").splitlines()
    require(int(xyz_lines[0]) == cell["atom_count"] and len(xyz_lines) == cell["atom_count"] + 2,
            "XYZ viewer export atom count mismatch")
    pilots = [] if args.skip_pilots else run_pilots(output_dir)
    from rdkit import __version__ as rdkit_version
    artifacts = {p.name: {"sha256": digest(p.read_bytes()), "bytes": p.stat().st_size}
                 for p in sorted(output_dir.iterdir()) if p.is_file()}
    report = {
        "schema_version": "1", "evidence_kind": "executed_geometry_construction_and_independent_audit",
        "status": "passed", "started_utc": started_utc,
        "finished_utc": datetime.now(timezone.utc).isoformat(),
        "elapsed_seconds": time.monotonic() - clock,
        "software": {"polychemmd": __version__, "python": platform.python_version(),
                     "python_implementation": platform.python_implementation(), "rdkit": rdkit_version},
        "auditor_source_sha256": digest(Path(__file__).read_bytes()),
        "commands_actually_executed": steps, "input_recipes": recipe_records,
        "conformer": {"atom_count": len(conformer_geometry.atoms),
                      "bond_count": len(conformer_geometry.bonds),
                      "source_graph_sha256": conformer_geometry.source_graph_sha256,
                      "canonical_payload_sha256": digest(canonical_bytes(conformer)),
                      "validation": conformer_geometry.provenance["validation"]},
        "cell_canonical_payload_sha256": digest(canonical_bytes(cell)),
        "canonical_payload_hash_basis": "sorted-key compact finite JSON encoded as UTF-8",
        "cell_independent_audit": audit, "additional_finite_graph_coordinate_pilots": pilots,
        "artifact_files": artifacts,
        "simulation_performed": False, "qm_calculation_performed": False,
        "production_force_field_assigned": False, "equilibrium_status": "not_assessed",
        "material_property_validation": "not_performed",
        "scientific_limits": [
            "Actual finite molecular coordinates and a geometric starting cell were constructed; no MD/QM ran.",
            "All reported density is prescribed-volume construction arithmetic, not a physical-property prediction.",
            "Conformer generation is not an energy-ranked search or a validated polymer conformational ensemble.",
            "PE, PVDF and PS pilots are H-capped finite repeat graphs; tacticity is not specified or inferred.",
            "XYZ is a viewer export; topology, images and provenance remain in the JSON companions.",
        ],
    }
    write_json(output_dir / "report.json", report)
    print(json.dumps({"status": report["status"], "report": str(output_dir / "report.json"),
                      "atoms_in_cell": cell["atom_count"], "molecules_in_cell": cell["molecule_count"],
                      "additional_coordinate_pilots": len(pilots)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
