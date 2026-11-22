"""Isolated, optional RadonPy reference job for the Studio. No external imports in-process."""

from __future__ import annotations

from functools import lru_cache
import json
import os
from pathlib import Path
import subprocess
import sys
import time

SUPPORTED = ("PEO", "PVA", "PVP", "PEtOx")
_TIMEOUT = 120
_WORKER = Path(__file__).with_name("_radon_runner.py").absolute()


def _interpreter() -> str:
    value = os.environ.get("POLYCHEMMD_RADON_PYTHON", sys.executable)
    executable = Path(value)
    if not executable.is_absolute() or not executable.is_file():
        raise ValueError(
            "Set POLYCHEMMD_RADON_PYTHON to an existing absolute Python executable."
        )
    # Resolving the venv symlink here would silently select the base installation.
    return str(executable)


@lru_cache(maxsize=8)
def _probe(executable: str, minute: int) -> dict:
    try:
        proc = subprocess.run(
            [executable, "-I", str(_WORKER), "version"],
            text=True,
            capture_output=True,
            timeout=15,
            check=True,
        )
        return {**json.loads(proc.stdout.splitlines()[-1]), "available": True}
    except (subprocess.SubprocessError, OSError, ValueError, IndexError) as exc:
        stderr = getattr(exc, "stderr", "") or ""
        reason = (
            str(stderr).strip().splitlines()[-1]
            if stderr.strip()
            else type(exc).__name__
        )
        return {
            "available": False,
            "reason": reason,
            "setup": "Install the optional RadonPy environment and set POLYCHEMMD_RADON_PYTHON.",
        }


def status() -> dict:
    try:
        return _probe(_interpreter(), int(time.monotonic() // 60)).copy()
    except ValueError as exc:
        return {"available": False, "reason": str(exc)}


def run_radon(config: dict, out: Path) -> dict:
    """Run the registered oligomer workflow; all paths and options are application controlled."""
    polymer = config.get("polymer", "PEO")
    dp = config.get("dp", 8)
    seed = config.get("seed", 17)
    if polymer not in SUPPORTED:
        raise ValueError("Supported RadonPy model chemistries: " + ", ".join(SUPPORTED))
    if type(dp) is not int or not 2 <= dp <= 16:
        raise ValueError("RadonPy DP must be an integer in [2,16].")
    if type(seed) is not int or not 1 <= seed <= 2147483646:
        raise ValueError("Seed must be a positive integer below 2^31−1.")
    availability = status()
    if not availability["available"]:
        raise RuntimeError(
            f"RadonPy unavailable: {availability.get('reason', '')}. {availability.get('setup', '')}"
        )
    out.mkdir(parents=True, exist_ok=True)
    request = dict(polymer=polymer, dp=dp, seed=seed)
    try:
        proc = subprocess.run(
            [_interpreter(), "-I", str(_WORKER), "oligomer"],
            cwd=out,
            input=json.dumps(request),
            text=True,
            capture_output=True,
            timeout=_TIMEOUT,
            check=False,
            env={**os.environ, "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"},
        )
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"RadonPy exceeded its {_TIMEOUT} second budget.") from exc
    (out / "radon.log").write_text(proc.stdout + "\n" + proc.stderr, encoding="utf-8")
    if proc.returncode:
        reason = (
            proc.stderr.strip().splitlines()[-1]
            if proc.stderr.strip()
            else f"exit {proc.returncode}"
        )
        raise RuntimeError(f"RadonPy oligomer job failed: {reason}")
    return json.loads((out / "radon_result.json").read_text(encoding="utf-8"))
