"""Bounded, artifact-checked LAMMPS process adapter.

Input scripts are trusted scientific recipes, not untrusted sandboxed programs.
This adapter checks process execution; it does not certify equilibration or any
material property. The supported log format is ordinary single-line thermo with
Step as its first column and a completed ``Loop time of`` block.
"""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import signal
import subprocess
import time


class LammpsError(RuntimeError):
    """Explicit unavailable, failed or incomplete execution; no fake result."""

    def __init__(self, message, record=None):
        super().__init__(message)
        self.record = deepcopy(record) if record is not None else None


def _name(value):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._/-]*", value):
        raise ValueError("artifact names must be safe relative POSIX paths")
    path = PurePosixPath(value)
    if path.is_absolute() or any(p in (".", "..") for p in value.split("/")):
        raise ValueError("artifact paths must not escape the run directory")
    return value


def _finite(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be finite")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class LammpsRunSpec:
    """Immutable execution recipe with an explicitly bounded resource request.

    ``input_files`` contains (relative_name, UTF-8 text) pairs, permitting a data
    file to be staged without pointing the solver at an earlier run directory.
    """

    input_script: str
    run_directory: Path | str
    executable: Path | str | None = None
    timeout_seconds: float = 30.0
    expected_artifacts: tuple[str, ...] = ("log.lammps",)
    required_thermo_columns: tuple[str, ...] = ("Step", "Temp")
    completion_marker: str = "POLYCHEMMD_RUN_COMPLETE"
    input_files: tuple[tuple[str, str], ...] = ()
    units_style: str | None = None

    def __post_init__(self):
        if not isinstance(self.input_script, str) or not self.input_script.strip():
            raise ValueError("input_script must be nonempty UTF-8 text")
        if "\x00" in self.input_script:
            raise ValueError("input_script must not contain NUL")
        if self.units_style is not None and self.units_style not in ("lj", "real", "metal", "si", "cgs", "electron", "micro", "nano"):
            raise ValueError("unsupported units_style declaration")
        unit_analysis = _units_analysis(self.input_script)
        known_styles = unit_analysis["run_styles"] or ([unit_analysis["style"]] if unit_analysis["style"] else [])
        if self.units_style is not None and any(style != self.units_style for style in known_styles):
            raise ValueError("units_style conflicts with the inferred units active at a run/reset in input_script")
        timeout = _finite(self.timeout_seconds, "timeout_seconds")
        if not 0 < timeout <= 3600:
            raise ValueError("timeout_seconds must be in (0, 3600]")
        object.__setattr__(self, "timeout_seconds", timeout)
        object.__setattr__(self, "run_directory", Path(self.run_directory).expanduser().resolve())
        if isinstance(self.expected_artifacts, str):
            raise ValueError("expected_artifacts must be a sequence of paths, not a string")
        expected = tuple(self.expected_artifacts)
        if not expected or len(set(expected)) != len(expected):
            raise ValueError("expected_artifacts must be nonempty and unique")
        for name in expected:
            _name(name)
        object.__setattr__(self, "expected_artifacts", expected)
        columns = tuple(self.required_thermo_columns)
        if not columns or len(set(columns)) != len(columns) or "Step" not in columns:
            raise ValueError("required_thermo_columns must be unique and include Step")
        if any(not isinstance(c, str) or not c or c.split() != [c] for c in columns):
            raise ValueError("invalid thermo column name")
        object.__setattr__(self, "required_thermo_columns", columns)
        if not isinstance(self.completion_marker, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]{7,127}", self.completion_marker):
            raise ValueError("completion_marker must be a single explicit identifier")
        reserved = {"input.lammps", "log.lammps", "screen.lammps", "stdout.txt", "stderr.txt", "run_record.json", "engine_help.txt"}
        inputs = tuple(tuple(item) for item in self.input_files)
        names = []
        for item in inputs:
            if len(item) != 2:
                raise ValueError("input_files must contain (relative_name, text) pairs")
            name, content = item
            _name(name)
            if name in reserved or name in expected or not isinstance(content, str):
                raise ValueError("staged inputs must be text and distinct from outputs/reserved names")
            names.append(name)
        if len(names) != len(set(names)):
            raise ValueError("staged input names must be unique")
        object.__setattr__(self, "input_files", inputs)


@dataclass(frozen=True)
class ThermoBlock:
    columns: tuple[str, ...]
    rows: tuple[tuple[float, ...], ...]

    def to_dict(self):
        return {"columns": list(self.columns), "rows": [list(r) for r in self.rows], "loop_complete": True}


@dataclass(frozen=True)
class LammpsRunResult:
    record: dict

    def to_dict(self):
        return deepcopy(self.record)


def discover_lammps(executable=None):
    """Find one executable; never import a solver or claim it is functional."""
    names = (executable,) if executable is not None else ("lmp", "lmp_serial", "lmp_mpi")
    for name in names:
        found = shutil.which(os.fspath(name))
        if found:
            path = Path(found).resolve()
            if path.is_file() and os.access(path, os.X_OK):
                return path
    raise LammpsError("LAMMPS executable unavailable; install a working solver or supply its path")


def parse_thermo_log(text, required_columns=("Step", "Temp")):
    """Parse completed single-line thermo blocks; reject nonfinite/truncated rows.

    Explicit WARNING lines are retained by run metadata and skipped here. YAML,
    multi-line thermo and arbitrary text inside an active block are unsupported.
    ``run 0`` is valid process evidence but not equilibration evidence.
    """
    if not isinstance(text, str):
        raise ValueError("thermo log must be text")
    required = tuple(required_columns)
    blocks, columns, rows = [], None, []
    for line_number, line in enumerate(text.splitlines(), 1):
        stripped = line.strip()
        if not stripped:
            continue
        tokens = stripped.split()
        if tokens[0] == "Step":
            if columns is not None:
                raise ValueError("thermo block lacks its Loop time completion")
            if len(tokens) != len(set(tokens)) or not set(required).issubset(tokens):
                raise ValueError(f"thermo header missing required columns or duplicated at line {line_number}")
            columns, rows = tuple(tokens), []
        elif columns is not None:
            if stripped.startswith("WARNING"):
                continue
            if stripped.startswith("Loop time of"):
                if not rows:
                    raise ValueError("completed thermo block contains no numeric rows")
                blocks.append(ThermoBlock(columns, tuple(rows)))
                columns, rows = None, []
                continue
            if len(tokens) != len(columns):
                raise ValueError(f"malformed thermo row at line {line_number}")
            try:
                values = tuple(float(t) for t in tokens)
            except (ValueError, OverflowError) as exc:
                raise ValueError(f"nonnumeric thermo row at line {line_number}") from exc
            if not all(math.isfinite(v) for v in values):
                raise ValueError(f"nonfinite thermo value at line {line_number}")
            if not values[0].is_integer() or values[0] < 0:
                raise ValueError(f"invalid thermo step at line {line_number}")
            if rows and values[0] < rows[-1][0]:
                raise ValueError("thermo steps decrease inside one block")
            rows.append(values)
    if columns is not None:
        raise ValueError("thermo log ends before Loop time completion")
    if not blocks:
        raise ValueError("no completed single-line thermo blocks found")
    return tuple(blocks)


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _units_analysis(script):
    """Conservative recipe inference, explicitly not an engine interpreter.

    ``clear`` resets units to LAMMPS' default LJ units. For a straight-line
    recipe retain units at run/minimize commands, rather than reporting a later
    setting as the units of earlier thermo blocks. Multiple run-unit systems or
    any dynamic/control-flow construct leave a scalar unit label unresolved.
    """
    known = {"lj", "real", "metal", "si", "cgs", "electron", "micro", "nano"}
    dynamic = {"include", "jump", "if", "variable", "next", "label", "python", "shell", "read_restart"}
    commands = []
    for line in script.splitlines():
        tokens = line.split("#", 1)[0].strip().split()
        if not tokens:
            continue
        if (tokens[0] in dynamic or "$" in line.split("#", 1)[0] or "&" in line.split("#", 1)[0]
                or (tokens[0] == "run" and "every" in tokens)):
            return {"style": None, "run_styles": [], "reason": "dynamic_or_unsupported_recipe_control"}
        commands.append(tokens)
    current, evidence, run_styles = "lj", False, []
    for tokens in commands:
        if tokens[0] == "clear":
            if len(tokens) != 1:
                return {"style": None, "run_styles": [], "reason": "nonliteral_clear_command"}
            current, evidence = "lj", True
        elif tokens[0] == "units":
            if len(tokens) != 2 or tokens[1] not in known:
                return {"style": None, "run_styles": [], "reason": "nonliteral_units_command"}
            current, evidence = tokens[1], True
        elif tokens[0] in ("run", "minimize"):
            run_styles.append(current)
    if run_styles and len(set(run_styles)) > 1:
        return {"style": None, "run_styles": run_styles, "reason": "multiple_run_unit_systems"}
    if not evidence:
        return {"style": None, "run_styles": run_styles, "reason": "no_explicit_units_or_clear_in_recipe"}
    return {"style": run_styles[0] if run_styles else current, "run_styles": run_styles,
            "reason": "straight_line_recipe_with_clear_reset_handling"}


def _literal_units(script):
    """Compatibility helper returning only a conservatively inferred label."""
    return _units_analysis(script)["style"]


def plan_lammps(spec):
    """Return a plan only: no files, solver execution or computed observables."""
    if not isinstance(spec, LammpsRunSpec):
        raise TypeError("spec must be LammpsRunSpec")
    units = _units_analysis(spec.input_script)
    return {
        "status": "planned", "simulation_performed": False,
        "run_directory": str(spec.run_directory),
        "input_sha256": _sha(spec.input_script.encode("utf-8")),
        "timeout_seconds": spec.timeout_seconds,
        "expected_artifacts": list(dict.fromkeys(("log.lammps",) + spec.expected_artifacts)),
        "required_thermo_columns": list(spec.required_thermo_columns),
        "material_validation": "not_performed",
        "declared_units_style": spec.units_style,
        "input_units_style": units["style"],
        "input_units_style_source": units["reason"],
        "units_style_confidence": "recipe_literal_not_engine_verified" if units["style"] is not None else "unresolved",
        "recipe_run_units_styles": units["run_styles"],
        "thermo_unit_conversion": "not_performed; values are raw LAMMPS output",
        "thermo_normalization": "not_inferred; inspect the trusted recipe and engine conventions",
    }


def _write_record(directory, record):
    temporary = directory / "run_record.json.tmp"
    temporary.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    os.replace(temporary, directory / "run_record.json")


def _execute(command, directory, timeout, stdout_path, stderr_path):
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        process = subprocess.Popen(command, cwd=directory, stdin=subprocess.DEVNULL,
                                   stdout=stdout, stderr=stderr,
                                   start_new_session=(os.name == "posix"))
        try:
            return process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            if os.name == "posix":
                try:
                    os.killpg(process.pid, signal.SIGKILL)
                except ProcessLookupError:
                    pass
            else:
                process.kill()
            process.wait()
            raise


def run_lammps(spec):
    """Execute and verify one new run; raise LammpsError on any failed gate.

    LAMMPS input is trusted and may itself perform file/process operations. This
    adapter is not a security sandbox. MPI launching and restart recovery are
    separate future orchestration features.
    """
    if not isinstance(spec, LammpsRunSpec):
        raise TypeError("spec must be LammpsRunSpec")
    executable = discover_lammps(spec.executable)
    directory = spec.run_directory
    try:
        directory.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise LammpsError("run_directory already exists; reuse/stale artifacts are forbidden") from exc
    except OSError as exc:
        raise LammpsError(f"cannot create new run_directory: {exc}") from exc
    record = {
        **plan_lammps(spec), "status": "preparing", "simulation_performed": False,
        "evidence_level": "process_contract", "material_validation": "not_performed",
        "executable": str(executable),
        "started_utc": datetime.now(timezone.utc).isoformat(),
        "staged_input_sha256": {name: _sha(content.encode("utf-8")) for name, content in spec.input_files},
        "runtime_environment": {name: os.environ[name] for name in ("OMP_NUM_THREADS", "UCX_TLS") if name in os.environ},
    }
    phase = "engine_identification"
    started = time.monotonic()
    try:
        record["executable_sha256"] = _sha(executable.read_bytes())
        help_result = subprocess.run([str(executable), "-h"], cwd=directory, stdin=subprocess.DEVNULL,
                                     capture_output=True, text=True,
                                     timeout=min(10.0, spec.timeout_seconds), check=False)
        help_text = help_result.stdout + help_result.stderr
        (directory / "engine_help.txt").write_text(help_text, encoding="utf-8")
        if help_result.returncode != 0 or not help_text.strip():
            raise LammpsError("LAMMPS help/version invocation failed or returned no identity")
        record["engine_help_sha256"] = _sha(help_text.encode("utf-8"))
        record["engine_help_first_line"] = next(x.strip() for x in help_text.splitlines() if x.strip())
        record["engine_version_banner"] = next((line.strip() for line in help_text.splitlines()
                                                 if "Large-scale Atomic" in line or re.match(r"^LAMMPS\s*\(\d", line.strip())), None)
        (directory / "input.lammps").write_bytes(spec.input_script.encode("utf-8"))
        for name, content in spec.input_files:
            path = directory / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content.encode("utf-8"))
        command = [str(executable), "-in", "input.lammps", "-log", "log.lammps", "-screen", "screen.lammps"]
        record["command"] = command
        phase = "engine_run"
        start_ns = time.time_ns()
        record["simulation_performed"] = True
        returncode = _execute(command, directory, spec.timeout_seconds,
                              directory / "stdout.txt", directory / "stderr.txt")
        record["returncode"] = returncode
        if returncode != 0:
            raise LammpsError(f"LAMMPS returned nonzero exit status {returncode}")
        phase = "artifact_validation"
        artifact_hashes = {}
        for name in record["expected_artifacts"]:
            path = directory / name
            if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(directory):
                raise LammpsError(f"expected current-run artifact missing or unsafe: {name}")
            stat = path.stat()
            if stat.st_size == 0:
                raise LammpsError(f"expected artifact is empty: {name}")
            # One-second tolerance accommodates filesystems with coarse timestamps.
            if stat.st_mtime_ns < start_ns - 1_000_000_000:
                raise LammpsError(f"expected artifact is stale: {name}")
            artifact_hashes[name] = _sha(path.read_bytes())
        log = (directory / "log.lammps").read_text(encoding="utf-8")
        lines = log.splitlines()
        if any(re.match(r"^ERROR(?::| on proc \d+:)", line.strip()) for line in lines):
            raise LammpsError("LAMMPS reported an ERROR line in the current-run log")
        markers = [i for i, line in enumerate(lines) if line.strip() == spec.completion_marker]
        loops = [i for i, line in enumerate(lines) if line.strip().startswith("Loop time of")]
        if not markers or not loops or markers[-1] <= loops[-1]:
            raise LammpsError("explicit completion marker after the final thermo run is missing")
        blocks = parse_thermo_log(log, spec.required_thermo_columns)
        record.update(status="completed", evidence_level="engine_execution",
                      artifact_sha256=artifact_hashes,
                      thermo_blocks=[block.to_dict() for block in blocks],
                      warning_lines=[line for line in lines if line.strip().startswith("WARNING")],
                      completion_marker=spec.completion_marker,
                      equilibration_validation="not_performed")
    except (OSError, UnicodeError, ValueError, subprocess.SubprocessError, LammpsError) as exc:
        record.update(status="failed", failure_phase=phase, failure_reason=str(exc),
                      finished_utc=datetime.now(timezone.utc).isoformat(),
                      elapsed_seconds=time.monotonic() - started)
        _write_record(directory, record)
        raise LammpsError(f"{phase}: {exc}", record=record) from exc
    record.update(finished_utc=datetime.now(timezone.utc).isoformat(),
                  elapsed_seconds=time.monotonic() - started)
    _write_record(directory, record)
    return LammpsRunResult(record)
