"""Finite periodic-alkane minimization/NVT/NPT preparation pilots.

Successful process execution is deliberately distinct from minimizer convergence
and thermodynamic equilibrium. This module supplies no density predictor or
automatic equilibration certificate.
"""

from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import re

from ..analysis.density import density_from_mass_volume
from ..engines.lammps import LammpsError, LammpsRunSpec, plan_lammps, run_lammps


BAR_PER_ATM = 1.01325
KCAL_TO_KJ = 4.184
_PINNED_REAL_MASS_DENSITY_FACTOR = 1.0/0.602214129
_THERMO_COLUMNS = ("Step", "Temp", "Press", "Volume", "Density", "PotEng", "TotEng",
                   "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "E_long",
                   "Lx", "Ly", "Lz")
_DUMP_COLUMNS = ("id", "mol", "q", "x", "y", "z", "ix", "iy", "iz", "fx", "fy", "fz")
_STATE_DUMP_COLUMNS = (*_DUMP_COLUMNS, "vx", "vy", "vz")
_MARKERS = ("POLYCHEMMD_MINIMIZE_BEGIN", "POLYCHEMMD_MINIMIZE_END",
            "POLYCHEMMD_NVT_BEGIN", "POLYCHEMMD_NVT_END",
            "POLYCHEMMD_NPT_BEGIN", "POLYCHEMMD_NPT_END")


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a positive finite number")
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be positive and representable") from exc
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return number


def _integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


@dataclass(frozen=True)
class PeriodicAlkaneProtocolSpec:
    temperature_k: float = 300.0
    pressure_bar: float = 1.0
    timestep_fs: float = 0.25
    minimize_max_iterations: int = 200
    minimize_max_evaluations: int = 1000
    nvt_steps: int = 100
    npt_steps: int = 100
    velocity_seed: int = 49217
    timeout_seconds: float = 60.0
    cutoff_angstrom: float = 10.0
    ewald_alpha_per_angstrom: float = 0.30
    reciprocal_nmax: int = 16

    def __post_init__(self):
        for name in ("temperature_k", "pressure_bar", "timestep_fs", "timeout_seconds",
                     "cutoff_angstrom", "ewald_alpha_per_angstrom"):
            object.__setattr__(self, name, _positive(getattr(self, name), name))
        if self.timestep_fs > 1.0:
            raise ValueError("this unconstrained explicit-H preparation pilot limits timestep_fs to 1.0")
        if self.timeout_seconds > 120:
            raise ValueError("the preparation pilot timeout must not exceed 120 seconds")
        if not math.isfinite(2*self.cutoff_angstrom):
            raise ValueError("twice the cutoff must be representable")
        _integer(self.minimize_max_iterations, "minimize_max_iterations", 1, 5000)
        _integer(self.minimize_max_evaluations, "minimize_max_evaluations", 1, 50000)
        _integer(self.nvt_steps, "nvt_steps", 1, 10000)
        _integer(self.npt_steps, "npt_steps", 1, 10000)
        if self.nvt_steps + self.npt_steps > 10000:
            raise ValueError("the preparation pilot permits at most 10000 total integration steps")
        _integer(self.velocity_seed, "velocity_seed", 1, 900000000)
        _integer(self.reciprocal_nmax, "reciprocal_nmax", 1, 20)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_mapping(cls, data):
        expected = {field.name for field in fields(cls)}
        if not isinstance(data, dict) or not set(data).issubset(expected):
            raise ValueError(f"protocol fields must be drawn from {sorted(expected)}")
        return cls(**data)


class PeriodicProtocolError(LammpsError):
    """A preparation process or its output audit failed explicitly."""


def _scientific_controls(spec):
    return {
        "boundary": "orthorhombic_three_axis_periodic",
        "potential": "sourced_OPLS_AA_subset_with_real_space_LJ_and_periodic_Ewald_Coulomb",
        "minimizer": "cg", "minimizer_maximum_atom_move_angstrom": .1,
        "minimize_energy_tolerance_dimensionless": 1e-6,
        "minimize_force_tolerance_kcal_mol_angstrom": 1e-6,
        "minimize_force_tolerance_kj_mol_angstrom": 1e-6*KCAL_TO_KJ,
        "thermostat_damping_fs": 100.0, "barostat_damping_fs": 1000.0,
        "engine_pressure_target_atm": spec.pressure_bar/BAR_PER_ATM,
        "pressure_conversion": "1 atm = 1.01325 bar",
        "engine_density_conversion_factor": _PINNED_REAL_MASS_DENSITY_FACTOR,
        "engine_density_conversion_source": "https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L175",
        "canonical_density_basis": "supplied g/mol masses divided by volume with exact SI Avogadro constant",
        "thermo_stride_steps": 10,
        "runtime_box_guard": "fix halt every minimization/integration step; explicit Boolean axis checks, hard error if any box length <= 2*cutoff",
        "runtime_bond_guard": "fix halt every minimization/integration step; hard error if bondmax >= cutoff/3, a sufficient bound keeping every 1-4 endpoint separation within cutoff",
        "minimization_budget_policy": "continue the bounded pilot, retaining the actual stop criterion; no convergence inferred",
        "ensemble_claim": "short preparation stages, not a verified thermodynamic ensemble",
    }


def _prepare(assigned, spec, run_directory, executable):
    from ..engines.periodic_alkane import build_periodic_alkane_setup
    if not isinstance(spec, PeriodicAlkaneProtocolSpec):
        raise ValueError("spec must be a PeriodicAlkaneProtocolSpec")
    setup, data = build_periodic_alkane_setup(
        assigned, cutoff_angstrom=spec.cutoff_angstrom,
        ewald_alpha_per_angstrom=spec.ewald_alpha_per_angstrom,
        reciprocal_nmax=spec.reciprocal_nmax)
    target, pressure = format(spec.temperature_k, '.17g'), format(spec.pressure_bar/BAR_PER_ATM, '.17g')
    script = setup.rstrip()+"\n"+"\n".join((
        "# Bounded preparation pilot; no automatic equilibrium/density validation.",
        "thermo 10",
        "thermo_style custom step temp press vol density pe etotal ebond eangle edihed evdwl ecoul elong lx ly lz",
        "thermo_modify norm no format float %.17g",
        f'variable polychemmd_box_invalid equal "(lx<={2*spec.cutoff_angstrom:.17g})||(ly<={2*spec.cutoff_angstrom:.17g})||(lz<={2*spec.cutoff_angstrom:.17g})"',
        "fix polychemmd_boxguard all halt 1 v_polychemmd_box_invalid != 0 error hard message yes",
        f"fix polychemmd_bondguard all halt 1 bondmax >= {spec.cutoff_angstrom/3:.17g} error hard message yes",
        "min_style cg", "min_modify dmax 0.1",
        'print "POLYCHEMMD_MINIMIZE_BEGIN"',
        f"minimize 1e-6 1e-6 {spec.minimize_max_iterations} {spec.minimize_max_evaluations}",
        'print "POLYCHEMMD_MINIMIZE_END"',
        "unfix polychemmd_boxguard",
        "unfix polychemmd_bondguard",
        "reset_timestep 0",
        "fix polychemmd_boxguard all halt 1 v_polychemmd_box_invalid != 0 error hard message yes",
        f"fix polychemmd_bondguard all halt 1 bondmax >= {spec.cutoff_angstrom/3:.17g} error hard message yes",
        f"timestep {spec.timestep_fs:.17g}",
        f"velocity all create {target} {spec.velocity_seed} mom yes rot yes dist gaussian",
        f"fix polychemmd_nvt all nvt temp {target} {target} 100.0",
        'print "POLYCHEMMD_NVT_BEGIN"', f"run {spec.nvt_steps}",
        'print "POLYCHEMMD_NVT_END"', "unfix polychemmd_nvt",
        f"fix polychemmd_npt all npt temp {target} {target} 100.0 iso {pressure} {pressure} 1000.0",
        'print "POLYCHEMMD_NPT_BEGIN"', f"run {spec.npt_steps}",
        'print "POLYCHEMMD_NPT_END"',
        "# Preserve active Nose-Hoover/barostat state before removing its fix.",
        "write_restart checkpoint.restart",
        "write_dump all custom checkpoint.dump id mol q x y z ix iy iz fx fy fz vx vy vz modify sort id format float %.17g",
        "unfix polychemmd_npt",
        "unfix polychemmd_boxguard", "unfix polychemmd_bondguard", "write_data prepared.data",
        "write_dump all custom final.dump id mol q x y z ix iy iz fx fy fz modify sort id format float %.17g",
        'print "POLYCHEMMD_RUN_COMPLETE"', ""))
    run_spec = LammpsRunSpec(script, run_directory, executable=executable,
                             timeout_seconds=spec.timeout_seconds,
                             expected_artifacts=("log.lammps", "prepared.data", "final.dump", "checkpoint.restart", "checkpoint.dump"),
                             required_thermo_columns=_THERMO_COLUMNS,
                             input_files=(("system.data", data),), units_style="real")
    record = assigned.to_dict()
    plan = {
        "schema_version": "1", "kind": "periodic_alkane_preparation_protocol",
        "status": "planned", "model_id": record["model_id"],
        "assignment_canonical_sha256": _sha(_canonical(record)),
        "protocol_spec": spec.to_dict(), "scientific_controls": _scientific_controls(spec),
        "engine_plan": plan_lammps(run_spec),
        "input_script": script, "staged_inputs": {"system.data": data},
        "input_files": {"input.lammps_sha256": _sha(script.encode()),
                        "system.data_sha256": _sha(data.encode())},
        "md_integration_steps_requested": spec.nvt_steps+spec.npt_steps,
        "simulation_performed": False, "equilibrium_status": "not_assessed",
        "production_ready": False, "density_prediction_validation": "not_performed",
        "material_property_validation": "not_performed",
    }
    return run_spec, plan, record


def plan_periodic_alkane_protocol(assigned, spec, *, run_directory, executable=None):
    """Return fully reviewable inputs; discover/run no engine and create no files."""
    return _prepare(assigned, spec, run_directory, executable)[1]


def _integer_token(token, name):
    if not re.fullmatch(r"[+-]?\d+", token):
        raise ValueError(f"{name} must be an integer token")
    return int(token)


def parse_periodic_final_dump(text, assigned_record, expected_step, *, include_velocities=False):
    """Parse one actual final frame; positions may differ from the starting cell."""
    lines = text.splitlines()
    atoms = assigned_record["atoms"]
    count = len(atoms)
    if len(lines) != count+9 or lines[0] != "ITEM: TIMESTEP" or lines[2] != "ITEM: NUMBER OF ATOMS":
        raise ValueError("final dump must contain exactly one complete frame")
    step = _integer_token(lines[1], "final step")
    if step != expected_step or _integer_token(lines[3], "atom count") != count:
        raise ValueError("final dump step or atom count differs from the requested protocol")
    if lines[4].split() != ["ITEM:", "BOX", "BOUNDS", "pp", "pp", "pp"]:
        raise ValueError("final dump must retain an orthorhombic three-axis periodic box")
    origin, lengths = [], []
    for line in lines[5:8]:
        tokens = line.split()
        if len(tokens) != 2:
            raise ValueError("final box requires two orthorhombic bounds per axis")
        low, high = map(float, tokens)
        if not all(math.isfinite(v) for v in (low,high)) or high <= low:
            raise ValueError("final box bounds must be positive-width finite values")
        origin.append(low); lengths.append(high-low)
    if not all(math.isfinite(v) and v > 0 for v in lengths):
        raise ValueError("final box lengths are unrepresentable")
    columns = _STATE_DUMP_COLUMNS if include_velocities else _DUMP_COLUMNS
    if lines[8].split() != ["ITEM:", "ATOMS", *columns]:
        raise ValueError("final dump columns differ from the required identity/coordinate/force fields")
    result = {}
    for line in lines[9:]:
        tokens = line.split()
        if len(tokens) != len(columns):
            raise ValueError("malformed final dump row")
        identity, molecule = (_integer_token(tokens[i], name) for i,name in ((0,"atom ID"),(1,"molecule ID")))
        if not 1 <= identity <= count or identity in result:
            raise ValueError("final atom IDs must be unique and preserve the original sequence")
        source = atoms[identity-1]
        charge = float(tokens[2]); position = [float(v) for v in tokens[3:6]]
        images = [_integer_token(v,"image flag") for v in tokens[6:9]]
        forces = [float(v)*KCAL_TO_KJ for v in tokens[9:12]]
        if not all(math.isfinite(v) for v in [charge,*position,*forces]):
            raise ValueError("final atom row contains nonfinite coordinates/charge/forces")
        if molecule != source["molecule_index"]+1 or not math.isclose(charge,source["charge_e"],rel_tol=0,abs_tol=1e-12):
            raise ValueError("final molecule identity or fixed partial charge changed")
        if any(not low <= value < low+length for value,low,length in zip(position,origin,lengths)):
            raise ValueError("final wrapped coordinates lie outside their half-open box")
        unwrapped = [value+image*length for value,image,length in zip(position,images,lengths)]
        if not all(math.isfinite(v) for v in unwrapped):
            raise ValueError("final unwrapped coordinates are unrepresentable")
        result[identity] = {"engine_id": identity, "cell_atom_id": source["cell_atom_id"],
                            "molecule_index": source["molecule_index"], "charge_e": charge,
                            "wrapped_coordinates_angstrom": position, "image_flags": images,
                            "coordinates_angstrom": unwrapped, "forces_kj_mol_angstrom": forces}
        if include_velocities:
            velocities = [float(v) for v in tokens[12:15]]
            if not all(math.isfinite(v) for v in velocities):
                raise ValueError("checkpoint velocities must be finite")
            result[identity]["velocities_angstrom_fs"] = velocities
    if set(result) != set(range(1,count+1)):
        raise ValueError("final dump lost source atoms")
    return {"step":step,"atom_count":count,
            "box":{"origin_angstrom":origin,"lengths_angstrom":lengths,"periodic":[True]*3},
            "atoms":[result[i] for i in range(1,count+1)]}


def _minimizer_status(segment):
    stopping = re.findall(r"^\s*Stopping criterion\s*=\s*(.+?)\s*$", segment, flags=re.MULTILINE)
    criterion = stopping[-1] if stopping else None
    lower = criterion.lower() if criterion else ""
    if lower == "force tolerance":
        status = "force_tolerance_met"
    elif lower == "energy tolerance":
        status = "energy_tolerance_met_force_unverified"
    elif "max iterations" in lower or "max force evaluations" in lower:
        status = "budget_exhausted"
    else:
        status = "inconclusive"
    counts = re.findall(r"Iterations, force evaluations\s*=\s*(\d+)\s+(\d+)",segment)
    return {"status":status,"stopping_criterion":criterion,
            "iterations":int(counts[-1][0]) if counts else None,
            "force_evaluations":int(counts[-1][1]) if counts else None,
            "force_tolerance_met":lower == "force tolerance"}


def _audit_outputs(run, spec, assignment, directory):
    log = (directory/"log.lammps").read_text(encoding="utf-8")
    lines = log.splitlines()
    positions = []
    for marker in _MARKERS:
        matches = [i for i,line in enumerate(lines) if line.strip() == marker]
        if len(matches) != 1:
            raise ValueError("protocol stage completion markers are missing or duplicated")
        positions.append(matches[0])
    if positions != sorted(positions):
        raise ValueError("protocol stage markers are out of order")
    blocks = run["thermo_blocks"]
    if len(blocks) != 3:
        raise ValueError("protocol requires exactly minimization, NVT and NPT thermo blocks")
    stages = {"minimization":_minimizer_status("\n".join(lines[positions[0]:positions[1]]))}
    expected_end = (None,spec.nvt_steps,spec.nvt_steps+spec.npt_steps)
    expected_start = (None,0,spec.nvt_steps)
    last = None
    masses = [atom["mass_g_mol"] for atom in assignment["atoms"]]
    for index,(name,block) in enumerate(zip(("minimization","nvt","npt"),blocks)):
        if tuple(block["columns"]) != _THERMO_COLUMNS or not block["rows"]:
            raise ValueError("protocol thermo columns/rows differ from the prescribed fields")
        rows = [dict(zip(block["columns"],values)) for values in block["rows"]]
        for row in rows:
            if not all(math.isfinite(v) for v in row.values()):
                raise ValueError("protocol thermo contains nonfinite output")
            if min(row[axis] for axis in ("Lx","Ly","Lz")) <= 2*spec.cutoff_angstrom:
                raise ValueError("observed periodic box violates the unique minimum-image cutoff regime")
            if row["Volume"] <= 0 or row["Density"] <= 0:
                raise ValueError("protocol volume and mass density must remain positive")
        if index:
            if rows[0]["Step"] != expected_start[index] or rows[-1]["Step"] != expected_end[index]:
                raise ValueError("an integration stage did not execute exactly its requested steps")
            steps = spec.nvt_steps if name == "nvt" else spec.npt_steps
            stages[name] = {"status":"executed","steps":steps,
                            "elapsed_time_fs":steps*spec.timestep_fs,
                            "thermo_rows":[{"step":int(row["Step"]),"temperature_k":row["Temp"],
                                            "pressure_bar":row["Press"]*BAR_PER_ATM,
                                            "density_g_cm3":density_from_mass_volume(masses,row["Volume"]),
                                            "engine_density_g_cm3":row["Density"],"volume_angstrom3":row["Volume"]}
                                           for row in rows]}
        last = rows[-1]
    header = (directory/"prepared.data").read_text(encoding="utf-8").splitlines()[0]
    if "units = real" not in header:
        raise ValueError("actual prepared.data header does not confirm LAMMPS real units")
    final = parse_periodic_final_dump((directory/"final.dump").read_text(encoding="utf-8"),
                                     assignment,spec.nvt_steps+spec.npt_steps)
    box = final["box"]["lengths_angstrom"]
    if min(box) <= 2*spec.cutoff_angstrom:
        raise ValueError("final box violates the unique minimum-image cutoff regime")
    final_by_id = {atom["cell_atom_id"]:atom for atom in final["atoms"]}
    for bond in assignment["interactions"]["bonds"]:
        length = math.dist(*(final_by_id[identity]["coordinates_angstrom"] for identity in bond["atom_ids"]))
        if not math.isfinite(length) or length <= 0 or length >= spec.cutoff_angstrom/3:
            raise ValueError("final unwrapped bond violates the supported bonded-image/cutoff bound")
    volume = math.prod(box)
    if not math.isfinite(volume) or not math.isclose(volume,last["Volume"],rel_tol=1e-10,abs_tol=1e-8):
        raise ValueError("final dump box and actual thermo volume disagree")
    density = density_from_mass_volume(masses,volume)
    expected_native_density = math.fsum(masses)/volume*_PINNED_REAL_MASS_DENSITY_FACTOR
    if not math.isclose(expected_native_density,last["Density"],rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError("actual mass/volume density disagrees with the pinned real-unit engine convention")
    final.update(temperature_k=last["Temp"],pressure_bar=last["Press"]*BAR_PER_ATM,
                 density_g_cm3=density,volume_angstrom3=volume,
                 engine_density_g_cm3=last["Density"],
                 density_conversion_audit={"canonical_basis":"exact SI Avogadro constant, 6.02214076e23 mol^-1",
                                           "engine_real_factor":_PINNED_REAL_MASS_DENSITY_FACTOR,
                                           "engine_factor_basis":"pinned LAMMPS real-unit mv2d = 1/0.602214129",
                                           "source":"https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L175",
                                           "expected_engine_density_g_cm3":expected_native_density,
                                           "engine_abs_difference_g_cm3":abs(expected_native_density-last["Density"])},
                 potential_energy_kj_mol=last["PotEng"]*KCAL_TO_KJ,
                 total_energy_kj_mol=last["TotEng"]*KCAL_TO_KJ)
    return stages,final,{"basis":"actual LAMMPS write_data header","header":header,
                         "confirmed_units_style":"real"}


def _write_report(directory, record):
    (directory/"protocol_record.json").write_text(json.dumps(record,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def run_periodic_alkane_protocol(assigned, spec, *, run_directory, executable=None):
    """Run one finite preparation pilot, preserving failed/inconclusive states."""
    run_spec,report,assignment = _prepare(assigned,spec,run_directory,executable)
    report.update(status="failed",kind="periodic_alkane_preparation_protocol_execution")
    try:
        run = run_lammps(run_spec).to_dict()
    except LammpsError as exc:
        report.update(failure_phase="engine",failure_reason=str(exc),engine_run=exc.record,
                      simulation_performed=bool(exc.record and exc.record.get("simulation_performed")))
        # An existing/stale directory or unavailable binary must not be mutated.
        if exc.record is not None and run_spec.run_directory.is_dir():
            _write_report(run_spec.run_directory,report)
        raise PeriodicProtocolError(str(exc),record=report) from exc
    report.update(engine_run=run,simulation_performed=True)
    try:
        stages,final,unit_check = _audit_outputs(run,spec,assignment,run_spec.run_directory)
        from .restart import checkpoint_manifest
        checkpoint = checkpoint_manifest(assignment,spec,run,run_spec.run_directory,
                                         step=spec.nvt_steps+spec.npt_steps)
        report.update(status="completed_preparation_pilot",stages=stages,final_state=final,
                      checkpoint=checkpoint,
                      equilibrium_status="insufficient_sampling",
                      md_integration_steps_executed=spec.nvt_steps+spec.npt_steps,
                      independent_output_unit_check=unit_check,
                      cutoff_regime_validation="runtime_hard_guard_and_observed_thermo_final_box_checks",
                      final_snapshot_interpretation="actual prepared MD state, separate from initial rigid packing record")
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as exc:
        report.update(failure_phase="protocol_output_audit",failure_reason=str(exc))
        _write_report(run_spec.run_directory,report)
        raise PeriodicProtocolError(str(exc),record=report) from exc
    _write_report(run_spec.run_directory,report)
    return report
