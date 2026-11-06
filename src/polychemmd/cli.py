"""Command-line entry points for the starter's implemented capabilities."""

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
        "md_adapter_status": "planned",
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
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(f"polychemmd: {exc}", file=sys.stderr)
        return 2
    return 0
