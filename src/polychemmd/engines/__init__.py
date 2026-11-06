"""Independent external-engine adapters; importing never launches a solver."""

from .lammps import (
    LammpsError, LammpsRunSpec, LammpsRunResult, discover_lammps,
    parse_thermo_log, plan_lammps, run_lammps,
)

__all__ = ["LammpsError", "LammpsRunSpec", "LammpsRunResult", "discover_lammps",
           "parse_thermo_log", "plan_lammps", "run_lammps"]
