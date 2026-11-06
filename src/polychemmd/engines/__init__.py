"""Independent external-engine adapters; importing never launches a solver."""

from .lammps import (
    LammpsError, LammpsRunSpec, LammpsRunResult, discover_lammps,
    parse_thermo_log, plan_lammps, run_lammps,
)
from .alkane import (AlkaneParityError, build_alkane_input,
                     prepare_alkane_singlepoint, verify_alkane_singlepoint)
from .periodic_alkane import (build_periodic_alkane_input, build_periodic_alkane_setup,
                             prepare_periodic_alkane_singlepoint, verify_periodic_alkane_singlepoint)

__all__ = ["LammpsError", "LammpsRunSpec", "LammpsRunResult", "discover_lammps",
           "parse_thermo_log", "plan_lammps", "run_lammps", "AlkaneParityError",
           "build_alkane_input", "prepare_alkane_singlepoint", "verify_alkane_singlepoint",
           "build_periodic_alkane_input", "build_periodic_alkane_setup",
           "prepare_periodic_alkane_singlepoint", "verify_periodic_alkane_singlepoint"]
