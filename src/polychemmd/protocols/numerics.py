"""Predeclared fixed-geometry cutoff/Ewald diagnostics with finite resources.

This protocol compares actual engine evaluations of one source geometry. Its
conditional development checks never establish thermal-state or bulk-property
convergence, and a successful review cannot enable production automatically.
"""

from copy import deepcopy
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import re
import time

from ..engines.lammps import LammpsError

_COMPONENTS = {"bond", "angle", "dihedral", "lj", "coulomb_realspace", "coulomb_kspace", "total"}
_SCOPE = "one_fixed_source_geometry_not_thermal_state_or_bulk_convergence"
_SHA = re.compile(r"[0-9a-f]{64}\Z")


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _digest(value):
    return hashlib.sha256(_canonical(value)).hexdigest()


def _number(value, name, *, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be finite numeric data")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is not representable") from exc
    if not math.isfinite(result) or positive and result <= 0:
        raise ValueError(f"{name} must be finite" + (" and positive" if positive else ""))
    return result


def _integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


@dataclass(frozen=True)
class NumericalQualificationSpec:
    cutoffs_angstrom: tuple = (10., 12., 14., 16., 17.)
    ewald_alpha_per_angstrom: float = .30
    reciprocal_nmax: int = 16
    refined_reciprocal_nmax: int = 20
    min_consecutive_cutoff_pairs: int = 2
    cutoff_energy_per_atom_limit_kj_mol: float = .02
    cutoff_force_limit_kj_mol_angstrom: float = .05
    cutoff_pressure_limit_bar: float = 20.
    reciprocal_energy_per_atom_limit_kj_mol: float = 1e-5
    reciprocal_force_limit_kj_mol_angstrom: float = 1e-4
    reciprocal_pressure_limit_bar: float = .01
    unchanged_component_abs_tolerance_kj_mol: float = 1e-7
    per_run_timeout_seconds: float = 60.
    total_wall_seconds: float = 180.
    max_engine_runs: int = 6

    def __post_init__(self):
        cutoffs = self.cutoffs_angstrom
        if not isinstance(cutoffs, (list, tuple)) or not 3 <= len(cutoffs) <= 8:
            raise ValueError("cutoffs_angstrom must contain 3 to 8 declared refinements")
        cutoffs = tuple(_number(v, "cutoff_angstrom", positive=True) for v in cutoffs)
        if any(v > 100 for v in cutoffs) or any(b <= a for a,b in zip(cutoffs, cutoffs[1:])):
            raise ValueError("cutoffs must strictly increase and cannot exceed 100 angstrom")
        object.__setattr__(self, "cutoffs_angstrom", cutoffs)
        for name in ("ewald_alpha_per_angstrom", "cutoff_energy_per_atom_limit_kj_mol",
                     "cutoff_force_limit_kj_mol_angstrom", "cutoff_pressure_limit_bar",
                     "reciprocal_energy_per_atom_limit_kj_mol", "reciprocal_force_limit_kj_mol_angstrom",
                     "reciprocal_pressure_limit_bar", "unchanged_component_abs_tolerance_kj_mol",
                     "per_run_timeout_seconds", "total_wall_seconds"):
            value = _number(getattr(self, name), name, positive=True)
            if value > 1e9:
                raise ValueError(f"{name} exceeds the supported numeric review range")
            object.__setattr__(self, name, value)
        if self.ewald_alpha_per_angstrom > 10:
            raise ValueError("ewald_alpha_per_angstrom cannot exceed 10")
        _integer(self.reciprocal_nmax, "reciprocal_nmax", 1, 19)
        _integer(self.refined_reciprocal_nmax, "refined_reciprocal_nmax", self.reciprocal_nmax+1, 20)
        _integer(self.min_consecutive_cutoff_pairs, "min_consecutive_cutoff_pairs", 2, len(cutoffs)-1)
        _integer(self.max_engine_runs, "max_engine_runs", len(cutoffs)+1, 9)
        if self.per_run_timeout_seconds > 120 or self.total_wall_seconds > 240:
            raise ValueError("review permits at most 120 seconds per run and 240 total wall seconds")
        if self.per_run_timeout_seconds > self.total_wall_seconds or self.total_wall_seconds <= 2:
            raise ValueError("total wall budget must exceed two seconds and cover one requested run timeout")

    def to_dict(self):
        result = asdict(self)
        result["cutoffs_angstrom"] = list(self.cutoffs_angstrom)
        return result

    @classmethod
    def from_mapping(cls, data):
        expected = {f.name for f in fields(cls)}
        if not isinstance(data, dict) or not set(data).issubset(expected):
            raise ValueError(f"numerical qualification fields must be drawn from {sorted(expected)}")
        return cls(**data)


class NumericalQualificationError(LammpsError):
    """An actual evaluation or declared resource/record audit failed."""


def _requests(spec):
    return [{"stage_id": f"cutoff_{i:02d}", "settings": {
                "cutoff_angstrom": cutoff, "ewald_alpha_per_angstrom": spec.ewald_alpha_per_angstrom,
                "reciprocal_nmax": spec.reciprocal_nmax}}
            for i,cutoff in enumerate(spec.cutoffs_angstrom)] + [
                {"stage_id": "reciprocal_refinement", "settings": {
                 "cutoff_angstrom": spec.cutoffs_angstrom[-1],
                 "ewald_alpha_per_angstrom": spec.ewald_alpha_per_angstrom,
                 "reciprocal_nmax": spec.refined_reciprocal_nmax}}]


def _validate_record(record, settings):
    if not isinstance(record, dict) or record.get("status") != "completed_static_thermo_singlepoint":
        raise ValueError("successful numerical records require completed static thermodynamic single points")
    if record.get("kind") != "periodic_alkane_static_thermo_singlepoint":
        raise ValueError("numerical records must be static thermodynamic engine records")
    _integer(record.get("md_integration_steps"),"md_integration_steps",0,0)
    if (record.get("simulation_performed") is not True or record.get("md_integration_steps") != 0
            or record.get("trajectory_production_performed") is not False
            or record.get("pressure_interpretation") != "static_virial_not_equilibrium_average"
            or record.get("initial_velocities") != "zeroed_and_output_verified"):
        raise ValueError("numerical records must retain the verified static zero-velocity evidence boundary")
    units = {"energy_units": "kJ/mol", "force_units": "kJ/(mol angstrom)",
             "pressure_units": "bar", "volume_units": "angstrom^3", "density_units": "g/cm^3"}
    if any(record.get(k) != v for k,v in units.items()):
        raise ValueError("numerical comparison requires canonical energy/force/pressure/volume/density units")
    for name in ("assignment_canonical_sha256", "geometry_canonical_sha256"):
        if not isinstance(record.get(name), str) or not _SHA.fullmatch(record[name]):
            raise ValueError(f"{name} must identify the evaluated fixed source")
    if not isinstance(record.get("model_id"), str) or not record["model_id"]:
        raise ValueError("record requires its explicit model identity")
    engine_run = record.get("engine_run")
    if not isinstance(engine_run,dict) or not isinstance(engine_run.get("engine_identity"),dict) or not engine_run["engine_identity"]:
        raise ValueError("successful numerical records require the recorded producer engine identity")
    _canonical(engine_run["engine_identity"])
    if engine_run.get("status") != "completed" or engine_run.get("simulation_performed") is not True:
        raise ValueError("static review requires a completed actual engine execution record")
    for name in ("temperature_k","kinetic_energy_kj_mol"):
        if _number(record.get(name),name) != 0.:
            raise ValueError("static records require zero temperature and kinetic energy")
    proof = record.get("independent_output_unit_check")
    if (not isinstance(proof,dict) or proof.get("confirmed_units_style") != "real"
            or not isinstance(proof.get("header"),str) or "units = real" not in proof["header"]
            or record.get("atom_coordinate_molecule_charge_image_identity") != "passed"
            or record.get("pressure_compute") != "compute polychemmd_virial all pressure NULL virial"):
        raise ValueError("static output identity, native unit or virial-only pressure proof is missing")
    count = _integer(record.get("atom_count"), "atom_count", 1, 8192)
    ids = record.get("atom_ids")
    if (not isinstance(ids, list) or len(ids) != count or any(not isinstance(v,str) or not v for v in ids)
            or len(set(ids)) != count):
        raise ValueError("numerical record atom IDs must be nonempty, unique and ordered")
    actual = record.get("ewald_settings")
    if not isinstance(actual, dict):
        raise ValueError("actual numerical settings are missing")
    expected = {"cutoff_angstrom": settings["cutoff_angstrom"],
                "alpha_per_angstrom": settings["ewald_alpha_per_angstrom"],
                "reciprocal_nmax": settings["reciprocal_nmax"]}
    fixed_expected = {"relative_engine_accuracy_request":1e-8,"surface_boundary":"tin_foil",
                      "realspace_coulomb_table":False,"lj_tail_correction":False}
    if set(actual) != set(expected)|set(fixed_expected):
        raise ValueError("recorded numerical controls do not match the supported static model")
    for name,value in fixed_expected.items():
        if isinstance(value,bool):
            matches = actual[name] is value
        elif isinstance(value,float):
            matches = _number(actual[name],name,positive=True) == value
        else:
            matches = actual[name] == value
        if not matches:
            raise ValueError("recorded fixed numerical control differs from the supported static model")
    for name, value in expected.items():
        if name == "reciprocal_nmax":
            _integer(actual.get(name), name, 1, 20)
        else:
            _number(actual.get(name), name, positive=True)
        if actual[name] != value:
            raise ValueError("actual evaluation settings differ from the predeclared ordered request")
    energies = record.get("energy_components_kj_mol")
    if not isinstance(energies, dict) or set(energies) != _COMPONENTS:
        raise ValueError("numerical record requires every energy component")
    energies = {k: _number(v, k) for k,v in energies.items()}
    try:
        summed = math.fsum(v for k,v in energies.items() if k != "total")
    except OverflowError as exc:
        raise ValueError("component-energy sum is unrepresentable") from exc
    if not math.isclose(summed, energies["total"], rel_tol=1e-12, abs_tol=1e-8):
        raise ValueError("recorded energy components do not sum to total")
    forces = record.get("forces_kj_mol_angstrom")
    if (not isinstance(forces, list) or len(forces) != count
            or any(not isinstance(row,(list,tuple)) or len(row)!=3 for row in forces)):
        raise ValueError("ordered numerical forces must be N by 3")
    forces = [[_number(v,"force") for v in row] for row in forces]
    pressure = _number(record.get("pressure_bar"), "pressure_bar")
    tensor = record.get("pressure_tensor_bar")
    if (not isinstance(tensor, list) or len(tensor)!=3
            or any(not isinstance(row,(list,tuple)) or len(row)!=3 for row in tensor)):
        raise ValueError("static pressure tensor must be a symmetric three by three matrix")
    tensor = [[_number(v,"pressure tensor") for v in row] for row in tensor]
    if (any(not math.isclose(tensor[i][j],tensor[j][i],rel_tol=0,abs_tol=1e-10)
            for i in range(3) for j in range(3))
            or not math.isclose(pressure,math.fsum(tensor[i][i]/3 for i in range(3)),rel_tol=1e-12,abs_tol=1e-8)):
        raise ValueError("scalar pressure and tensor symmetry/trace disagree")
    volume = _number(record.get("volume_angstrom3"),"volume_angstrom3",positive=True)
    density = _number(record.get("density_g_cm3"),"density_g_cm3",positive=True)
    return {"energies":energies,"forces":forces,"pressure":pressure,"tensor":tensor,
            "volume":volume,"density":density,"count":count}


def _difference(first, second, name):
    result = second-first
    if not math.isfinite(result):
        raise ValueError(f"{name} difference exceeds the finite numeric range")
    return result


def _comparison(first, second, first_values, second_values, spec, *, reciprocal=False):
    diffs = {k:_difference(first_values["energies"][k],second_values["energies"][k],k) for k in _COMPONENTS}
    force_diff = [_difference(a,b,"force") for row1,row2 in zip(first_values["forces"],second_values["forces"])
                  for a,b in zip(row1,row2)]
    maximum = max(abs(v) for v in force_diff)
    rms = maximum*math.sqrt(math.fsum((v/maximum)**2 for v in force_diff)/len(force_diff)) if maximum else 0.
    pressure_diff = abs(_difference(first_values["pressure"],second_values["pressure"],"pressure"))
    tensor_max = max(abs(_difference(a,b,"pressure tensor")) for row1,row2 in zip(first_values["tensor"],second_values["tensor"])
                     for a,b in zip(row1,row2))
    invariants = {"bond","angle","dihedral"}
    invariants |= {"lj","coulomb_realspace"} if reciprocal else {"coulomb_kspace"}
    checks = {"unchanged_components": all(abs(diffs[k]) <= spec.unchanged_component_abs_tolerance_kj_mol for k in invariants)}
    if reciprocal:
        energy_limit,force_limit,pressure_limit = (spec.reciprocal_energy_per_atom_limit_kj_mol,
            spec.reciprocal_force_limit_kj_mol_angstrom,spec.reciprocal_pressure_limit_bar)
    else:
        energy_limit,force_limit,pressure_limit = (spec.cutoff_energy_per_atom_limit_kj_mol,
            spec.cutoff_force_limit_kj_mol_angstrom,spec.cutoff_pressure_limit_bar)
    per_atom = abs(diffs["total"])/first_values["count"]
    checks.update(energy_per_atom=per_atom<=energy_limit,maximum_force_component=maximum<=force_limit,
                  scalar_static_pressure=pressure_diff<=pressure_limit,static_pressure_tensor=tensor_max<=pressure_limit)
    return {"from_settings":deepcopy(first["ewald_settings"]),"to_settings":deepcopy(second["ewald_settings"]),
            "signed_component_energy_differences_kj_mol":diffs,"absolute_total_energy_difference_per_atom_kj_mol":per_atom,
            "max_force_component_difference_kj_mol_angstrom":maximum,"rms_force_component_difference_kj_mol_angstrom":rms,
            "absolute_scalar_static_pressure_difference_bar":pressure_diff,"max_static_pressure_tensor_component_difference_bar":tensor_max,
            "gates":checks,"status":"diagnostics_passed" if all(checks.values()) else "requires_review",
            "interpretation":("finite reciprocal refinement at fixed cutoff/splitting; no infinite-sum bound" if reciprocal
                else "cutoff refinement changes the truncated Hamiltonian; no exact-model error claim")}


def compare_numerical_records(records, spec):
    """Review every ordered request; an omitted or failed case cannot be selected away."""
    if not isinstance(spec,NumericalQualificationSpec):
        raise ValueError("spec must be a NumericalQualificationSpec")
    requests = _requests(spec)
    if not isinstance(records,(list,tuple)) or len(records)>len(requests):
        raise ValueError("records must follow the declared requests without extra evaluations")
    result = {"schema_version":"1","kind":"periodic_numerical_record_review",
              "numerical_status":"incomplete_records","scope":_SCOPE,"spec":spec.to_dict(),
              "n_expected_records":len(requests),"n_supplied_records":len(records),
              "case_outcomes":[],"cutoff_comparisons":[],"reciprocal_comparison":None,
              "selected_settings":None,"production_ready":False,"equilibrium_status":"not_assessed",
              "bulk_property_convergence":"not_established","density_prediction_validation":"not_performed"}
    values, successful = [], True
    for request,record in zip(requests,records):
        if not isinstance(record,dict):
            raise ValueError("each numerical record must be a mapping")
        try:
            _canonical(record)
        except (TypeError,ValueError,OverflowError,RecursionError) as exc:
            raise ValueError("numerical records must contain finite JSON data") from exc
        status = record.get("status")
        if status not in ("completed_static_thermo_singlepoint","failed","unavailable","timeout","inconclusive","planned","running"):
            raise ValueError("numerical record status is missing or unsupported")
        reason = record.get("failure_reason")
        if reason is not None and not isinstance(reason,str):
            raise ValueError("numerical record failure_reason must be text or null")
        result["case_outcomes"].append({"stage_id":request["stage_id"],"status":status,
            "failure_reason":record.get("failure_reason") if status!="completed_static_thermo_singlepoint" else None})
        if status != "completed_static_thermo_singlepoint":
            successful = False; values.append(None)
        else:
            values.append(_validate_record(record,request["settings"]))
    if len(records)!=len(requests) or not successful:
        _canonical(result)
        return result
    first = records[0]
    identity = ("assignment_canonical_sha256","geometry_canonical_sha256","model_id","atom_ids","atom_count")
    engine_identity = first["engine_run"]["engine_identity"]
    fixed_settings = {k:v for k,v in first["ewald_settings"].items() if k not in ("cutoff_angstrom","alpha_per_angstrom","reciprocal_nmax")}
    for index,record in enumerate(records[1:],1):
        if record["engine_run"]["engine_identity"] != engine_identity:
            raise ValueError("producer engine identity/configuration changed during the numerical sweep")
        if any(record[k] != first[k] for k in identity):
            raise ValueError("numerical records differ in assignment, model, geometry or ordered atom identity")
        if {k:v for k,v in record["ewald_settings"].items() if k not in ("cutoff_angstrom","alpha_per_angstrom","reciprocal_nmax")} != fixed_settings:
            raise ValueError("fixed nonbonded controls changed during the numerical sweep")
        if any(not math.isclose(values[index][k],values[0][k],rel_tol=1e-12,abs_tol=1e-12) for k in ("volume","density")):
            raise ValueError("fixed-geometry volume or mass density changed")
    cutoff_count = len(spec.cutoffs_angstrom)
    comparisons = [_comparison(records[i],records[i+1],values[i],values[i+1],spec)
                   for i in range(cutoff_count-1)]
    reciprocal = _comparison(records[cutoff_count-1],records[-1],values[cutoff_count-1],values[-1],spec,reciprocal=True)
    tail = comparisons[-spec.min_consecutive_cutoff_pairs:]
    cutoff_pass = all(item["status"]=="diagnostics_passed" for item in tail)
    invariant_pass = all(item["gates"]["unchanged_components"] for item in comparisons)
    passed = cutoff_pass and invariant_pass and reciprocal["status"]=="diagnostics_passed"
    magnitudes = [item["absolute_total_energy_difference_per_atom_kj_mol"] for item in comparisons]
    result.update(numerical_status="fixed_geometry_diagnostics_passed" if passed else "requires_review",
        assignment_canonical_sha256=first["assignment_canonical_sha256"], geometry_canonical_sha256=first["geometry_canonical_sha256"],
        model_id=first["model_id"], atom_count=first["atom_count"], engine_identity=deepcopy(engine_identity), cutoff_comparisons=comparisons,
        reciprocal_comparison=reciprocal, finest_consecutive_cutoff_pairs_required=spec.min_consecutive_cutoff_pairs,
        finest_consecutive_cutoff_pairs_passed=cutoff_pass, all_unchanged_component_invariants_passed=invariant_pass,
        absolute_cutoff_energy_changes_monotone_nonincreasing=all(b<=a for a,b in zip(magnitudes,magnitudes[1:])),
        nonmonotonic_changes_are_retained=True,
        selected_settings=deepcopy(requests[-1]["settings"]) if passed else None)
    _canonical(result)
    return result


def _prepare(assigned,spec,run_directory,executable):
    from ..engines.periodic_thermo import plan_periodic_alkane_thermo_singlepoint
    if not isinstance(spec,NumericalQualificationSpec):
        raise ValueError("spec must be a NumericalQualificationSpec")
    directory = Path(run_directory).expanduser().resolve()
    stages,common,identity = [],None,None
    for request in _requests(spec):
        target = directory/request["stage_id"]
        plan = plan_periodic_alkane_thermo_singlepoint(assigned,run_directory=target,
            executable=executable,timeout_seconds=spec.per_run_timeout_seconds,**request["settings"])
        this_identity = {k:deepcopy(plan[k]) for k in ("assignment_canonical_sha256","geometry_canonical_sha256","model_id","atom_count")}
        if identity is None:
            identity = this_identity; common = plan["staged_inputs"]["system.data"]
        elif this_identity != identity or plan["staged_inputs"]["system.data"] != common:
            raise ValueError("planned numerical stages changed the fixed source or common engine data")
        stages.append({**request,"run_directory":str(target),"status":"planned",
                       "input_script":plan["input_script"],"input_files":deepcopy(plan["input_files"]),
                       "engine_plan":deepcopy(plan["engine_plan"])})
    report = {"schema_version":"1","kind":"periodic_numerical_qualification","status":"planned",
        **identity,"spec":spec.to_dict(),"spec_canonical_sha256":_digest(spec.to_dict()),
        "qualification_inputs_sha256":_digest({**identity,"spec":spec.to_dict()}),
        "scope":_SCOPE,"run_directory":str(directory),"stages":stages,
        "shared_input":{"relative_path":"common/system.data","sha256":hashlib.sha256(common.encode()).hexdigest()},
        "expected_engine_runs":len(stages),"engine_requests_started":0,"engine_execution_attempts_reported":0,"simulation_performed":False,
        "numerical_status":"not_assessed","selected_settings":None,"production_ready":False,
        "equilibrium_status":"not_assessed","bulk_property_convergence":"not_established",
        "density_prediction_validation":"not_performed","material_property_validation":"not_performed"}
    return directory,report,common


def plan_numerical_qualification(assigned,spec,*,run_directory,executable=None):
    """Create a reviewable plan only; no engine discovery, directories or writes."""
    _,report,common = _prepare(assigned,spec,run_directory,executable)
    report["staged_inputs"] = {"common/system.data":common}
    return report


def _persist(directory,report):
    text = json.dumps(report,indent=2,allow_nan=False)+"\n"
    temporary = directory/"numerical_qualification.json.tmp"
    temporary.write_text(text,encoding="utf-8")
    temporary.replace(directory/"numerical_qualification.json")


def run_numerical_qualification(assigned,spec,*,run_directory,executable=None):
    """Execute all declared static evaluations with an aggregate wall budget."""
    from ..engines.periodic_thermo import run_periodic_alkane_thermo_singlepoint
    started = time.monotonic()
    directory,report,common = _prepare(assigned,spec,run_directory,executable)
    if time.monotonic()-started >= spec.total_wall_seconds-2:
        raise NumericalQualificationError("planning exhausted the declared numerical review budget",record=report)
    try:
        directory.mkdir(parents=True,exist_ok=False)
    except OSError as exc:
        raise NumericalQualificationError(f"cannot create a fresh numerical review directory: {exc}") from exc
    (directory/"common").mkdir()
    (directory/"common/system.data").write_text(common,encoding="utf-8")
    report.update(status="running",started_monotonic_elapsed_seconds=time.monotonic()-started)
    _persist(directory,report)  # Exact settings, scripts, hashes and limits precede the first process.
    records = []
    try:
        for stage in report["stages"]:
            remaining = spec.total_wall_seconds-(time.monotonic()-started)
            if remaining <= 2:
                raise NumericalQualificationError("aggregate numerical qualification wall budget exhausted before another engine call")
            timeout = min(spec.per_run_timeout_seconds,remaining-2)
            report["engine_requests_started"] += 1
            if report["engine_requests_started"] > spec.max_engine_runs:
                raise NumericalQualificationError("predeclared engine-run budget exhausted")
            stage.update(status="running",effective_timeout_seconds=timeout)
            _persist(directory,report)
            try:
                record = run_periodic_alkane_thermo_singlepoint(assigned,
                    run_directory=stage["run_directory"],executable=executable,
                    timeout_seconds=timeout,**stage["settings"])
            except LammpsError as exc:
                stage.update(status="failed",failure_reason=str(exc),record=exc.record)
                if exc.record is not None:
                    records.append(exc.record)
                raise
            stage.update(status=record.get("status"),record=record)
            records.append(record)
            report["engine_execution_attempts_reported"] = sum(
                bool(item.get("record",{}).get("engine_run",{}).get("simulation_performed"))
                for item in report["stages"] if isinstance(item.get("record"),dict))
            report["simulation_performed"] = report["simulation_performed"] or bool(record.get("simulation_performed"))
            if any(record.get(k) != report[k] for k in ("assignment_canonical_sha256", "geometry_canonical_sha256", "model_id", "atom_count")):
                raise ValueError("actual numerical record differs from the predeclared source identity")
            if record.get("input_files") != stage["input_files"]:
                raise ValueError("actual single-point input hashes differ from the predeclared plan")
            _persist(directory,report)
            if time.monotonic()-started > spec.total_wall_seconds:
                raise NumericalQualificationError("aggregate numerical qualification wall budget exceeded")
        review = compare_numerical_records(records,spec)
        report.update(status="completed_numerical_review",review=review,
                      numerical_status=review["numerical_status"],selected_settings=review["selected_settings"],
                      engine_identity=deepcopy(review["engine_identity"]))
        report["measured_total_wall_seconds"] = time.monotonic()-started
        if report["measured_total_wall_seconds"] > spec.total_wall_seconds:
            raise NumericalQualificationError("aggregate numerical wall budget exceeded during review postprocessing")
        _persist(directory,report)
        if time.monotonic()-started > spec.total_wall_seconds:
            raise NumericalQualificationError("aggregate numerical wall budget exceeded during final record serialization")
    except (OSError,ValueError,KeyError,TypeError,ArithmeticError,LammpsError) as exc:
        try:
            partial_review = compare_numerical_records(records,spec)
        except (ValueError, TypeError, KeyError, ArithmeticError) as review_error:
            partial_review = {"numerical_status":"invalid_records", "failure_reason":str(review_error)}
        report.update(status="failed",numerical_status="execution_failed",selected_settings=None,
                      failure_reason=str(exc),partial_review=partial_review,
                      measured_total_wall_seconds=time.monotonic()-started)
        report["simulation_performed"] = any(bool(s.get("record",{}).get("simulation_performed")) for s in report["stages"] if isinstance(s.get("record"),dict))
        report["engine_execution_attempts_reported"] = sum(
            bool(item.get("record",{}).get("engine_run",{}).get("simulation_performed"))
            for item in report["stages"] if isinstance(item.get("record"),dict))
        _persist(directory,report)
        raise NumericalQualificationError(str(exc),record=report) from exc
    return report
