"""Command-line entry points for the independently implemented research core."""

import argparse
import importlib.resources
import importlib.util
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import sys
import tempfile

from . import __version__
from .records import MaterialSpec
from .provenance import envelope


def _read(path):
    raw = Path(path).read_bytes()
    def reject_nonfinite(value):
        raise ValueError(f"non-finite JSON value: {value}")
    def unique_object(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise ValueError(f"duplicate JSON field: {key}")
            result[key] = value
        return result
    return raw, json.loads(raw, parse_constant=reject_nonfinite, object_pairs_hook=unique_object)


def _check_output_path(input_path, output_path):
    if output_path is None:
        return
    source, target = Path(input_path), Path(output_path)
    if source.resolve() == target.resolve() or (target.exists() and source.samefile(target)):
        raise ValueError("output must not overwrite the input file or an alias of it")


def _emit(data, output=None):
    text = json.dumps(data, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    _emit_text(text, output)


def _emit_text(text, output=None):
    if output:
        target = Path(output)
        target.parent.mkdir(parents=True, exist_ok=True)
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=target.parent,
                                             prefix=".polychemmd-", delete=False) as stream:
                temporary = Path(stream.name)
                stream.write(text)
                stream.flush()
                os.fsync(stream.fileno())
            os.replace(temporary, target)
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
    else:
        print(text, end="")


def _chain_from_mapping(data):
    from .chemistry import RepeatUnitSpec
    from .builders import build_linear
    expected = {"repeat_unit", "degree_of_polymerization", "head_cap", "tail_cap"}
    if not isinstance(data, dict) or set(data) != expected:
        raise ValueError(f"linear build requires exactly {sorted(expected)}")
    return build_linear(RepeatUnitSpec.from_mapping(data["repeat_unit"]),
                        data["degree_of_polymerization"],
                        head_cap=data["head_cap"], tail_cap=data["tail_cap"])


def _seed(value, name):
    if type(value) is not int or not 0 <= value <= 2**31-1:
        raise ValueError(f"{name} must be an integer in [0, 2147483647]")
    return value


def _positive_recipe_number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be positive and finite")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be positive and finite") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _geometry_from_record(data):
    from .geometry import MolecularGeometry
    if isinstance(data, dict) and "record_kind" in data and "payload" in data:
        data = data["payload"]
    return MolecularGeometry.from_mapping(data)


def _periodic_assignment_from_record(data, *, compact=False):
    from .forcefields import (AssignedPeriodicAlkaneCell, assign_periodic_alkanes,
                             CompactAssignedPeriodicAlkaneCell, assign_compact_periodic_alkanes)
    if isinstance(data, dict) and "record_kind" in data and "payload" in data:
        data = data["payload"]
    if isinstance(data, dict) and data.get("kind") == "compact_assigned_periodic_opls_aa_n_alkane_cell":
        return CompactAssignedPeriodicAlkaneCell(json.dumps(data, allow_nan=False))
    if isinstance(data, dict) and data.get("kind") == "assigned_periodic_opls_aa_n_alkane_cell":
        if compact:
            original = AssignedPeriodicAlkaneCell(json.dumps(data, allow_nan=False))
            return assign_compact_periodic_alkanes(original.to_dict()["source_cell"])
        return AssignedPeriodicAlkaneCell(json.dumps(data, allow_nan=False))
    return (assign_compact_periodic_alkanes if compact else assign_periodic_alkanes)(data)


def _prepare_cell(data):
    from .builders.conformers import generate_conformer
    from .builders.packing import pack_cell
    expected = {"chain", "conformer_seeds", "packing_seed", "box_lengths_angstrom",
                "min_distance_angstrom", "max_attempts_per_molecule", "max_distance_checks"}
    if not isinstance(data, dict) or set(data) != expected:
        raise ValueError(f"pack-cell input requires exactly {sorted(expected)}")
    seeds = data["conformer_seeds"]
    if not isinstance(seeds, list) or not 1 <= len(seeds) <= 128:
        raise ValueError("conformer_seeds must contain 1 to 128 explicit seeds")
    for seed in seeds:
        _seed(seed, "conformer seed")
    if len(set(seeds)) != len(seeds):
        raise ValueError("conformer_seeds must be distinct; these are starting conformers, not equilibrated cell replicates")
    _seed(data["packing_seed"], "packing_seed")
    box = data["box_lengths_angstrom"]
    if not isinstance(box, list) or len(box) != 3:
        raise ValueError("box_lengths_angstrom must have three positive lengths")
    box = [_positive_recipe_number(v, "box length") for v in box]
    cutoff = _positive_recipe_number(data["min_distance_angstrom"], "min_distance_angstrom")
    if any(length < cutoff for length in box):
        raise ValueError("box length must be at least the geometric cutoff")
    if not math.isfinite(math.prod(box)) or math.prod(box) <= 0:
        raise ValueError("box volume must be positive and representable")
    for name in ("max_attempts_per_molecule", "max_distance_checks"):
        if type(data[name]) is not int or data[name] <= 0:
            raise ValueError(f"{name} must be a positive integer")
    chain = _chain_from_mapping(data["chain"])
    if len(chain.lineage)*len(seeds) > 10000:
        raise ValueError("total atom count exceeds the pack-cell CLI limit of 10000")
    geometries = [generate_conformer(chain, seed=seed) for seed in seeds]
    payload = pack_cell(geometries, box_lengths_angstrom=box, seed=data["packing_seed"],
                        min_distance_angstrom=cutoff,
                        max_attempts_per_molecule=data["max_attempts_per_molecule"],
                        max_distance_checks=data["max_distance_checks"]).to_dict()
    payload["preparation"] = {
        "chain_input": data["chain"], "conformer_seeds": seeds,
        "interpretation": "separately seeded starting conformers, not independent equilibrated cells",
    }
    return payload


def doctor():
    packages = {}
    for name in ("numpy", "scipy", "rdkit", "psi4", "pyscf"):
        try:
            packages[name] = importlib.util.find_spec(name) is not None
        except (ImportError, ValueError):
            packages[name] = False
    return {
        "software_version": __version__,
        "optional_python_packages": packages,
        "optional_executables": {name: shutil.which(name) for name in ("lmp", "lmp_serial", "gmx", "xtb", "cp2k.psmp")},
        "scope": "Presence check only; this does not verify a working engine or an implemented adapter.",
        "md_adapter_status": "implemented_process_adapter",
        "polymer_md_protocol_status": "bounded_periodic_preparation_and_restart",
        "density_sampling_status": "diagnostics_only_no_equilibrium_certificate",
        "periodic_assignment_representations": {"legacy_exhaustive_v1": 512, "compact_sparse_pairs_v2": 8192},
        "polymer_equilibrium_validation": "not_implemented",
        "qm_adapter_status": "planned",
    }


def main(argv=None):
    parser = argparse.ArgumentParser(prog="polychemmd", description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    commands.add_parser("doctor", help="Report optional dependency presence; execute no simulations")
    commands.add_parser("catalogue", help="Show planned property families and evidence status")
    validate = commands.add_parser("validate", help="Validate a material metadata record (not molecular valence)")
    validate.add_argument("input")
    validate.add_argument("--output")
    torsion = commands.add_parser("torsion-describe", help="Compute an analytic Fourier-torsion descriptor demo")
    torsion.add_argument("input")
    torsion.add_argument("--n-grid", type=int, default=720)
    torsion.add_argument("--temperature-k", type=float, default=300.)
    torsion.add_argument("--output")
    repeat = commands.add_parser("validate-repeat", help="Validate a mapped, already transformed repeat unit")
    repeat.add_argument("input")
    repeat.add_argument("--output")
    build = commands.add_parser("build-linear", help="Assemble a supported finite molecular graph; no coordinates or MD")
    build.add_argument("input")
    build.add_argument("--output")
    conformer = commands.add_parser("build-conformer", help="Generate one seeded, checked 3D realization; no MD or energy ranking")
    conformer.add_argument("input")
    conformer.add_argument("--output")
    pack = commands.add_parser("pack-cell", help="Prepare rigid molecules in an orthorhombic periodic cell; no equilibration")
    pack.add_argument("input")
    pack.add_argument("--output")
    tiled = commands.add_parser("tile-cell", help="Construct an ordered molecular motif cell with periodic contact checks; phase unvalidated")
    tiled.add_argument("input", help="Molecular geometry JSON or result envelope")
    tiled.add_argument("--recipe", required=True, help="Explicit repeats, spacings, rotation and geometric budgets")
    tiled.add_argument("--output")
    xyz = commands.add_parser("export-xyz", help="Export coordinates for viewing; topology and provenance remain in companion JSON")
    xyz.add_argument("input")
    xyz.add_argument("--wrapped", action="store_true", help="Use wrapped positions for packed cells")
    xyz.add_argument("--output")
    assignment = commands.add_parser("parameterize-alkane", help="Assign the frozen OPLS-AA numerical subset to a neutral explicit n-alkane")
    assignment.add_argument("input", help="Molecular geometry JSON or its result envelope")
    assignment.add_argument("--output")
    evaluate = commands.add_parser("alkane-energy", help="Evaluate isolated finite-cutoff alkane energies and analytic forces; no engine")
    evaluate.add_argument("input")
    evaluate.add_argument("--cutoff-angstrom", type=float, default=12.0)
    evaluate.add_argument("--output")
    parity = commands.add_parser("alkane-lammps", help="Verify a run-zero LAMMPS single point against independent energies and forces")
    parity.add_argument("input")
    parity.add_argument("--run-dir", required=True)
    parity.add_argument("--executable")
    parity.add_argument("--timeout", type=float, default=30.)
    parity.add_argument("--cutoff-angstrom", type=float, default=12.0)
    parity.add_argument("--plan", action="store_true", help="Return generated inputs without running or writing an engine directory")
    parity.add_argument("--output")
    cell_assignment = commands.add_parser("parameterize-cell", help="Assemble supported n-alkane assignments with periodic atom/molecule/image identity")
    cell_assignment.add_argument("input")
    cell_assignment.add_argument("--compact", action="store_true", help="Use sparse special pairs and an 8192-atom cell budget")
    cell_assignment.add_argument("--output")
    periodic_energy = commands.add_parser("periodic-energy", help="Evaluate the explicit neutral periodic Ewald verification model")
    periodic_energy.add_argument("input")
    periodic_energy.add_argument("--cutoff-angstrom", type=float, default=10.0)
    periodic_energy.add_argument("--ewald-alpha", type=float, default=.30)
    periodic_energy.add_argument("--reciprocal-nmax", type=int, default=16)
    periodic_energy.add_argument("--output")
    periodic_engine = commands.add_parser("periodic-lammps", help="Compare actual periodic Ewald energies/forces or review generated inputs")
    periodic_engine.add_argument("input")
    periodic_engine.add_argument("--run-dir", required=True)
    periodic_engine.add_argument("--executable")
    periodic_engine.add_argument("--timeout", type=float, default=60.)
    periodic_engine.add_argument("--cutoff-angstrom", type=float, default=10.)
    periodic_engine.add_argument("--ewald-alpha", type=float, default=.30)
    periodic_engine.add_argument("--reciprocal-nmax", type=int, default=16)
    periodic_engine.add_argument("--plan", action="store_true")
    periodic_engine.add_argument("--output")
    singlepoint = commands.add_parser("periodic-singlepoint", help="Run an engine-only periodic single point; no independent parity claim")
    singlepoint.add_argument("input")
    singlepoint.add_argument("--run-dir", required=True)
    singlepoint.add_argument("--executable")
    singlepoint.add_argument("--timeout", type=float, default=60.)
    singlepoint.add_argument("--cutoff-angstrom", type=float, default=10.)
    singlepoint.add_argument("--ewald-alpha", type=float, default=.30)
    singlepoint.add_argument("--reciprocal-nmax", type=int, default=16)
    singlepoint.add_argument("--plan", action="store_true")
    singlepoint.add_argument("--output")
    preparation = commands.add_parser("periodic-prepare", help="Plan/run bounded minimization and NVT/NPT preparation; no equilibrium certification")
    preparation.add_argument("input")
    preparation.add_argument("--protocol", required=True, help="Explicit preparation-spec JSON")
    preparation.add_argument("--run-dir", required=True)
    preparation.add_argument("--executable")
    preparation.add_argument("--plan", action="store_true")
    preparation.add_argument("--output")
    continuation = commands.add_parser("periodic-resume", help="Plan/run a hash-bound native NPT restart; no reseeding or equilibrium claim")
    continuation.add_argument("input", help="Initial assigned cell or its source packed cell")
    continuation.add_argument("--checkpoint", required=True, help="Checkpoint manifest or preparation/continuation result JSON")
    continuation.add_argument("--checkpoint-dir", required=True, help="Directory containing the native checkpoint binary")
    continuation.add_argument("--protocol", required=True, help="Bounded continuation-spec JSON")
    continuation.add_argument("--run-dir", required=True)
    continuation.add_argument("--executable")
    continuation.add_argument("--plan", action="store_true")
    continuation.add_argument("--output")
    campaign = commands.add_parser("density-plan", help="Review a declared state, independent cells, sampling gates and campaign budget; run nothing")
    campaign.add_argument("input")
    campaign.add_argument("--output")
    sampling = commands.add_parser("density-sampling", help="Check a physical-time series against predeclared density diagnostics; no equilibrium certificate")
    sampling.add_argument("input")
    sampling.add_argument("--protocol", required=True, help="Explicit density sampling gates JSON")
    sampling.add_argument("--output")
    convert = commands.add_parser("torsion-import", help="Explicitly convert a supported LAMMPS torsion convention")
    convert.add_argument("input")
    convert.add_argument("--source", choices=("lammps-fourier", "lammps-opls"), required=True)
    convert.add_argument("--output")
    coverage = commands.add_parser("parameter-check", help="Check assigned-type parameter coverage; not physical validation")
    coverage.add_argument("input")
    coverage.add_argument("--require-complete", action="store_true")
    coverage.add_argument("--output")
    engine = commands.add_parser("lammps-run", help="Execute a trusted script with bounded, current-run artifact checks")
    engine.add_argument("input")
    engine.add_argument("--run-dir", required=True)
    engine.add_argument("--executable")
    engine.add_argument("--timeout", type=float, default=30.)
    engine.add_argument("--expect", action="append", default=[])
    engine.add_argument("--plan", action="store_true")
    engine.add_argument("--output")
    density = commands.add_parser("density-analyze", help="Analyze supplied g/cm^3 data; no equilibrium or material validation claim")
    density.add_argument("input")
    density.add_argument("--block-size", type=int, required=True)
    density.add_argument("--discard", type=int, default=0)
    density.add_argument("--output")
    args = parser.parse_args(argv)
    try:
        if args.command == "doctor":
            _emit(doctor())
        elif args.command == "catalogue":
            raw = importlib.resources.files("polychemmd").joinpath("property_catalogue.json").read_text(encoding="utf-8")
            _emit(json.loads(raw))
        elif args.command == "validate":
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            record = MaterialSpec.from_mapping(data)
            _emit(envelope({"record": record.to_dict(), "metadata_validation": "passed",
                            "molecular_graph_validation": "not_performed", "simulation_performed": False},
                           raw, "material_metadata_validation"), args.output)
        elif args.command == "torsion-describe":
            from .descriptors.torsion import TorsionPotential
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            result = TorsionPotential.from_mapping(data).describe(n_grid=args.n_grid, temperature_k=args.temperature_k)
            _emit(envelope({"source_kind": "analytic_demo", "material_validated": False,
                            "simulation_performed": False, "descriptor": result}, raw,
                           "analytic_torsion_descriptor"), args.output)
        elif args.command == "validate-repeat":
            from .chemistry import RepeatUnitSpec
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            spec = RepeatUnitSpec.from_mapping(data)
            _emit(envelope({"repeat_unit": spec.to_dict(), "graph_validation": "passed", "simulation_performed": False},
                           raw, "repeat_unit_validation"), args.output)
        elif args.command == "build-linear":
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            result = _chain_from_mapping(data)
            payload = result.to_dict()
            payload["simulation_performed"] = False
            _emit(envelope(payload, raw, "linear_molecular_graph"), args.output)
        elif args.command == "build-conformer":
            from .builders.conformers import generate_conformer
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            if not isinstance(data, dict) or set(data) != {"chain", "seed", "max_iterations"}:
                raise ValueError("build-conformer input requires exactly chain, seed and max_iterations")
            _seed(data["seed"], "seed")
            if type(data["max_iterations"]) is not int or not 0 < data["max_iterations"] <= 2**31-1:
                raise ValueError("max_iterations must be a positive signed 32-bit integer")
            geometry = generate_conformer(_chain_from_mapping(data["chain"]),
                                          seed=data["seed"], max_iterations=data["max_iterations"])
            _emit(envelope(geometry.to_dict(), raw, "molecular_conformer"), args.output)
        elif args.command == "pack-cell":
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            _emit(envelope(_prepare_cell(data), raw, "periodic_cell_preparation"), args.output)
        elif args.command == "tile-cell":
            from .builders.tiling import tile_molecular_geometry
            for source in (args.input, args.recipe):
                _check_output_path(source, args.output)
            raw, data = _read(args.input)
            recipe_raw, recipe = _read(args.recipe)
            expected = {"repeats", "spacing_angstrom", "rotation_matrix", "min_distance_angstrom",
                        "max_atoms", "max_distance_checks"}
            if not isinstance(recipe, dict) or set(recipe) != expected:
                raise ValueError(f"tiling recipe requires exactly {sorted(expected)}")
            payload = tile_molecular_geometry(_geometry_from_record(data), **recipe).to_dict()
            payload["recipe_input_sha256"] = hashlib.sha256(recipe_raw).hexdigest()
            _emit(envelope(payload, raw, "ordered_molecular_cell_construction"), args.output)
        elif args.command == "export-xyz":
            from .geometry import MolecularGeometry
            from .io.xyz import molecular_xyz, cell_xyz
            _check_output_path(args.input, args.output)
            _, data = _read(args.input)
            if isinstance(data, dict) and "payload" in data and "record_kind" in data:
                data = data["payload"]
            if not isinstance(data, dict):
                raise ValueError("export-xyz requires a molecular geometry or packed-cell JSON record")
            if data.get("kind") == "molecular_geometry":
                if args.wrapped:
                    raise ValueError("--wrapped applies only to periodic cells")
                text = molecular_xyz(MolecularGeometry.from_mapping(data))
            elif data.get("kind") == "orthorhombic_rigid_packed_cell":
                text = cell_xyz(data, wrapped=args.wrapped)
            else:
                raise ValueError("unsupported structure kind for XYZ export")
            _emit_text(text, args.output)
        elif args.command in ("parameterize-alkane", "alkane-energy", "alkane-lammps"):
            from .forcefields import assign_alkane, evaluate_alkane
            _check_output_path(args.input, args.output)
            if args.command == "alkane-lammps" and args.output is not None:
                if Path(args.output).resolve().is_relative_to(Path(args.run_dir).resolve()):
                    raise ValueError("summary output must be outside the engine run directory to preserve its artifacts")
            raw, data = _read(args.input)
            assigned = assign_alkane(_geometry_from_record(data))
            if args.command == "parameterize-alkane":
                _emit(envelope(assigned.to_dict(), raw, "alkane_parameter_assignment"), args.output)
            elif args.command == "alkane-energy":
                _emit(envelope(evaluate_alkane(assigned, cutoff_angstrom=args.cutoff_angstrom),
                               raw, "alkane_energy_force_evaluation"), args.output)
            else:
                from .engines import plan_lammps, prepare_alkane_singlepoint, verify_alkane_singlepoint
                if args.plan:
                    spec = prepare_alkane_singlepoint(assigned, run_directory=args.run_dir,
                                                     executable=args.executable, timeout_seconds=args.timeout,
                                                     cutoff_angstrom=args.cutoff_angstrom)
                    payload = {"model_id": assigned.to_dict()["model_id"],
                               "engine_plan": plan_lammps(spec), "input_script": spec.input_script,
                               "staged_inputs": dict(spec.input_files),
                               "md_integration_steps": 0, "material_property_validation": "not_performed"}
                    kind = "alkane_singlepoint_plan"
                else:
                    payload = verify_alkane_singlepoint(assigned, run_directory=args.run_dir,
                                                       executable=args.executable, timeout_seconds=args.timeout,
                                                       cutoff_angstrom=args.cutoff_angstrom)
                    kind = "alkane_engine_energy_force_parity"
                _emit(envelope(payload, raw, kind), args.output)
        elif args.command in ("parameterize-cell", "periodic-energy", "periodic-lammps", "periodic-singlepoint", "periodic-prepare"):
            _check_output_path(args.input, args.output)
            if args.command in ("periodic-lammps", "periodic-singlepoint", "periodic-prepare"):
                if args.output is not None and Path(args.output).resolve().is_relative_to(Path(args.run_dir).resolve()):
                    raise ValueError("summary output must be outside the engine run directory to preserve its artifacts")
            protocol_raw = protocol_data = None
            if args.command == "periodic-prepare":
                _check_output_path(args.protocol, args.output)
                protocol_raw, protocol_data = _read(args.protocol)
            raw, data = _read(args.input)
            assigned = _periodic_assignment_from_record(data, compact=getattr(args, "compact", False))
            if args.command == "parameterize-cell":
                payload, kind = assigned.to_dict(), "periodic_alkane_parameter_assignment"
            elif args.command == "periodic-energy":
                from .forcefields import evaluate_periodic_alkanes
                payload = evaluate_periodic_alkanes(assigned, cutoff_angstrom=args.cutoff_angstrom,
                                                   ewald_alpha_per_angstrom=args.ewald_alpha,
                                                   reciprocal_nmax=args.reciprocal_nmax)
                kind = "periodic_ewald_energy_force_evaluation"
            elif args.command in ("periodic-lammps", "periodic-singlepoint"):
                from .engines import (plan_lammps, prepare_periodic_alkane_singlepoint,
                                      verify_periodic_alkane_singlepoint, run_periodic_alkane_singlepoint)
                options = {"run_directory": args.run_dir, "executable": args.executable,
                           "timeout_seconds": args.timeout, "cutoff_angstrom": args.cutoff_angstrom,
                           "ewald_alpha_per_angstrom": args.ewald_alpha,
                           "reciprocal_nmax": args.reciprocal_nmax}
                if args.plan:
                    prepare_options = dict(options)
                    if args.command == "periodic-singlepoint":
                        prepare_options["native_readback"] = True
                    spec = prepare_periodic_alkane_singlepoint(assigned, **prepare_options)
                    payload = {"engine_plan": plan_lammps(spec), "input_script": spec.input_script,
                               "staged_inputs": dict(spec.input_files), "md_integration_steps": 0,
                               "bulk_production_model_validated": False}
                    if args.command == "periodic-singlepoint":
                        payload["independent_reference_status"] = "not_performed"
                    kind = "periodic_ewald_engine_singlepoint_plan" if args.command == "periodic-singlepoint" else "periodic_ewald_singlepoint_plan"
                else:
                    function = (run_periodic_alkane_singlepoint if args.command == "periodic-singlepoint"
                                else verify_periodic_alkane_singlepoint)
                    payload = function(assigned, **options)
                    kind = "periodic_ewald_engine_singlepoint" if args.command == "periodic-singlepoint" else "periodic_ewald_engine_parity"
            else:
                from .protocols.alkane import (PeriodicAlkaneProtocolSpec,
                                              plan_periodic_alkane_protocol, run_periodic_alkane_protocol)
                spec = PeriodicAlkaneProtocolSpec.from_mapping(protocol_data)
                function = plan_periodic_alkane_protocol if args.plan else run_periodic_alkane_protocol
                payload = function(assigned, spec, run_directory=args.run_dir, executable=args.executable)
                payload["protocol_input_sha256"] = hashlib.sha256(protocol_raw).hexdigest()
                kind = "periodic_preparation_plan" if args.plan else "bounded_periodic_preparation"
            _emit(envelope(payload, raw, kind), args.output)
        elif args.command == "periodic-resume":
            from .protocols.restart import (PeriodicAlkaneContinuationSpec,
                plan_periodic_alkane_continuation, run_periodic_alkane_continuation)
            for source in (args.input, args.checkpoint, args.protocol):
                _check_output_path(source, args.output)
            if args.output is not None and any(Path(args.output).resolve().is_relative_to(Path(directory).resolve())
                                               for directory in (args.run_dir, args.checkpoint_dir)):
                raise ValueError("summary output must be outside engine/checkpoint directories to preserve their artifacts")
            raw, data = _read(args.input)
            checkpoint_raw, checkpoint = _read(args.checkpoint)
            protocol_raw, protocol = _read(args.protocol)
            if isinstance(checkpoint, dict) and "record_kind" in checkpoint and "payload" in checkpoint:
                checkpoint = checkpoint["payload"]
            if isinstance(checkpoint, dict) and "checkpoint" in checkpoint:
                checkpoint = checkpoint["checkpoint"]
            spec = PeriodicAlkaneContinuationSpec.from_mapping(protocol)
            function = plan_periodic_alkane_continuation if args.plan else run_periodic_alkane_continuation
            payload = function(_periodic_assignment_from_record(data), checkpoint, spec,
                               checkpoint_directory=args.checkpoint_dir,
                               run_directory=args.run_dir, executable=args.executable)
            payload["checkpoint_input_sha256"] = hashlib.sha256(checkpoint_raw).hexdigest()
            payload["protocol_input_sha256"] = hashlib.sha256(protocol_raw).hexdigest()
            _emit(envelope(payload, raw, "periodic_continuation_plan" if args.plan else "bounded_periodic_continuation"), args.output)
        elif args.command == "density-plan":
            from .protocols.campaign import DensityCampaignSpec, plan_density_campaign
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            _emit(envelope(plan_density_campaign(DensityCampaignSpec.from_mapping(data)), raw,
                           "density_campaign_plan"), args.output)
        elif args.command == "density-sampling":
            from .analysis.sampling import DensitySamplingSpec, analyze_density_sampling
            for source in (args.input, args.protocol):
                _check_output_path(source, args.output)
            raw, data = _read(args.input)
            protocol_raw, protocol = _read(args.protocol)
            expected = {"units", "rows", "source_reference", "source_kind"}
            if not isinstance(data, dict) or set(data) != expected:
                raise ValueError(f"sampling input requires exactly {sorted(expected)}")
            if data["units"] != {"density": "g/cm^3", "temperature": "K", "pressure": "bar"}:
                raise ValueError("sampling rows require explicit g/cm^3, K and bar units")
            if data["source_kind"] not in ("analytic_demo", "simulation", "experiment") or not isinstance(data["source_reference"], str) or not data["source_reference"].strip():
                raise ValueError("provide a valid source_kind and nonempty source_reference")
            payload = analyze_density_sampling(data["rows"], DensitySamplingSpec.from_mapping(protocol))
            payload["source_kind_reported_by_input"] = data["source_kind"]
            payload["source_reference_reported_by_input"] = data["source_reference"]
            payload["protocol_input_sha256"] = hashlib.sha256(protocol_raw).hexdigest()
            _emit(envelope(payload, raw, "density_sampling_diagnostics"), args.output)
        elif args.command == "torsion-import":
            from .descriptors.conventions import import_lammps_fourier, import_lammps_opls
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            potential = (import_lammps_fourier if args.source == "lammps-fourier" else import_lammps_opls)(data)
            _emit(envelope({"source_convention": args.source,
                            "converted_potential": {"convention": "sum_k_1_plus_cos_nphi_minus_phase",
                                                    "terms": [term.as_dict() for term in potential.terms]},
                            "parameter_physical_validation": "not_performed", "simulation_performed": False},
                           raw, "torsion_convention_conversion"), args.output)
        elif args.command == "parameter-check":
            from .forcefields import ParameterSet, TopologyRequirements, check_coverage
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            if not isinstance(data, dict) or set(data) != {"parameter_set", "topology"}:
                raise ValueError("parameter-check input requires exactly parameter_set and topology")
            report = check_coverage(ParameterSet.from_mapping(data["parameter_set"]),
                                    TopologyRequirements.from_mapping(data["topology"]))
            if args.require_complete:
                report.require_complete()
            _emit(envelope(report.to_dict(), raw, "assigned_parameter_coverage"), args.output)
        elif args.command == "lammps-run":
            from .engines import LammpsRunSpec, plan_lammps, run_lammps
            _check_output_path(args.input, args.output)
            if args.output is not None and Path(args.output).resolve().is_relative_to(Path(args.run_dir).resolve()):
                raise ValueError("summary output must be outside the engine run directory to preserve its artifacts")
            raw = Path(args.input).read_bytes()
            spec = LammpsRunSpec(input_script=raw.decode("utf-8"), run_directory=args.run_dir,
                                 executable=args.executable, timeout_seconds=args.timeout,
                                 expected_artifacts=tuple(dict.fromkeys(["log.lammps"]+args.expect)))
            result = plan_lammps(spec) if args.plan else run_lammps(spec).to_dict()
            _emit(envelope(result, raw, "engine_execution_plan" if args.plan else "engine_execution"), args.output)
        elif args.command == "density-analyze":
            from .analysis import analyze_density
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            if not isinstance(data, dict) or set(data) != {"units", "densities", "source_reference", "source_kind"}:
                raise ValueError("density input requires units, densities, source_reference and source_kind")
            if data["units"] != "g/cm^3":
                raise ValueError("density-analyze requires explicitly converted g/cm^3; reduced LJ units are unsupported")
            if data["source_kind"] not in ("analytic_demo", "simulation", "experiment") or not isinstance(data["source_reference"], str) or not data["source_reference"].strip():
                raise ValueError("provide a valid density source_kind and nonempty source_reference")
            result = analyze_density(data["densities"], block_size=args.block_size, discard=args.discard)
            result["source_kind_reported_by_input"] = data["source_kind"]
            result["source_reference_reported_by_input"] = data["source_reference"]
            _emit(envelope(result, raw, "density_diagnostics"), args.output)
    except (OSError, ValueError, TypeError, KeyError, RuntimeError) as exc:
        print(f"polychemmd: {exc}", file=sys.stderr)
        return 2
    return 0
