"""Command-line entry points for the independently implemented research core."""

import argparse
import importlib.resources
import importlib.util
import json
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
        "polymer_md_protocol_status": "not_implemented",
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
            from .chemistry import RepeatUnitSpec
            from .builders import build_linear
            _check_output_path(args.input, args.output)
            raw, data = _read(args.input)
            expected = {"repeat_unit", "degree_of_polymerization", "head_cap", "tail_cap"}
            if not isinstance(data, dict) or set(data) != expected:
                raise ValueError(f"linear build requires exactly {sorted(expected)}")
            result = build_linear(RepeatUnitSpec.from_mapping(data["repeat_unit"]),
                                  data["degree_of_polymerization"],
                                  head_cap=data["head_cap"], tail_cap=data["tail_cap"])
            payload = result.to_dict()
            payload["simulation_performed"] = False
            _emit(envelope(payload, raw, "linear_molecular_graph"), args.output)
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
