"""Predeclared preparation review, deliberately weaker than equilibration.

Thermal forces are not minimization residuals.  The actual minimizer stopping
criterion and any explicitly requested gross-force screen have separate roles.
Input records are audited for consistency, not cryptographically authenticated.
"""

from dataclasses import asdict, dataclass, field, fields
import hashlib
import json
import math
import re

from .density import density_from_mass_volume
from .sampling import DensitySamplingSpec, analyze_density_sampling


_MODEL = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"
_PREPARATION_KIND = "periodic_alkane_preparation_protocol_execution"
_GEOMETRY_ATOM_KEYS = ("cell_atom_id", "molecule_index", "charge_e", "coordinates_angstrom",
                       "wrapped_coordinates_angstrom", "image_flags")


def _canonical(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()


def _sha(data):
    return hashlib.sha256(_canonical(data)).hexdigest()


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is unrepresentable") from exc
    if not math.isfinite(value) or positive and value <= 0 or nonnegative and value < 0:
        raise ValueError(f"{name} is outside its finite allowed range")
    return value


def _integer(value, name, minimum=0):
    if type(value) is not int or not minimum <= value <= 10**12:
        raise ValueError(f"{name} must be an integer in [{minimum}, 1000000000000]")
    return value


def _close(a, b, name):
    if not math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10):
        raise ValueError(f"{name} is inconsistent with the actual protocol/state")


@dataclass(frozen=True)
class PreparationReadinessSpec:
    phase_target: str = "liquid_candidate"
    temperature_k: float = 300.0
    pressure_bar: float = 1.0
    timestep_fs: float = 0.25
    require_minimizer_force_convergence: bool = True
    max_force_norm_kj_mol_angstrom: float | None = None
    max_rms_force_kj_mol_angstrom: float | None = None
    max_nematic_order: float = 0.20
    sampling_spec: DensitySamplingSpec = field(default_factory=DensitySamplingSpec)

    def __post_init__(self):
        if self.phase_target != "liquid_candidate":
            raise ValueError("this review supports only the explicitly unverified liquid_candidate target")
        for name in ("temperature_k", "pressure_bar", "timestep_fs"):
            object.__setattr__(self, name, _number(getattr(self, name), name, positive=True))
        if self.timestep_fs > 1:
            raise ValueError("the supported unconstrained-H timestep must not exceed 1 fs")
        if type(self.require_minimizer_force_convergence) is not bool:
            raise ValueError("require_minimizer_force_convergence must be an explicit boolean")
        for name in ("max_force_norm_kj_mol_angstrom", "max_rms_force_kj_mol_angstrom"):
            value = getattr(self, name)
            if value is not None:
                object.__setattr__(self, name, _number(value, name, positive=True))
        maximum = _number(self.max_nematic_order, "max_nematic_order", nonnegative=True)
        if maximum > 1:
            raise ValueError("max_nematic_order must be in [0, 1]")
        object.__setattr__(self, "max_nematic_order", maximum)
        sampling = self.sampling_spec
        if isinstance(sampling, dict):
            sampling = DensitySamplingSpec.from_mapping(sampling)
        if not isinstance(sampling, DensitySamplingSpec):
            raise ValueError("sampling_spec must be a DensitySamplingSpec or its mapping")
        for actual, expected, name in ((sampling.timestep_fs, self.timestep_fs, "sampling timestep"),
                (sampling.target_temperature_k, self.temperature_k, "sampling temperature target"),
                (sampling.target_pressure_bar, self.pressure_bar, "sampling pressure target")):
            if actual != expected:
                raise ValueError(f"{name} must match the readiness declaration")
        object.__setattr__(self, "sampling_spec", sampling)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_mapping(cls, data):
        names = {item.name for item in fields(cls)}
        if not isinstance(data, dict) or not set(data).issubset(names):
            raise ValueError(f"preparation-readiness fields must be drawn from {sorted(names)}")
        return cls(**data)


def _assignment(assigned):
    from ..forcefields.periodic import AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell
    if not isinstance(assigned, (AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell)):
        raise ValueError("assigned must be a validated periodic n-alkane assignment")
    result = assigned.to_dict()
    if result.get("model_id") != _MODEL:
        raise ValueError("unsupported preparation-review model")
    return result


def _stage(stage, expected_start, expected_steps, timestep_fs, masses, name):
    if not isinstance(stage, dict) or stage.get("status") != "executed":
        raise ValueError(f"{name} requires an actual executed stage record")
    if _integer(stage.get("steps"), f"{name} steps", 1) != expected_steps:
        raise ValueError(f"{name} stage step budget differs from its protocol")
    elapsed = _number(stage.get("elapsed_time_fs"), f"{name} elapsed time", positive=True)
    _close(elapsed, expected_steps*timestep_fs, f"{name} physical duration")
    supplied = stage.get("thermo_rows")
    if not isinstance(supplied, list) or len(supplied) < 2 or len(supplied) > 100000:
        raise ValueError(f"{name} requires at least its initial/final thermo rows")
    rows = []
    for source in supplied:
        if not isinstance(source, dict):
            raise ValueError(f"{name} thermo rows must be mappings")
        row = {"step": _integer(source.get("step"), f"{name} thermo step"),
               "temperature_k": _number(source.get("temperature_k"), "temperature", nonnegative=True),
               "pressure_bar": _number(source.get("pressure_bar"), "pressure"),
               "density_g_cm3": _number(source.get("density_g_cm3"), "density", positive=True)}
        volume = _number(source.get("volume_angstrom3"), "thermo volume", positive=True)
        _close(row["density_g_cm3"], density_from_mass_volume(masses, volume), "canonical thermo density")
        rows.append(row)
    if rows[0]["step"] != expected_start or rows[-1]["step"] != expected_start+expected_steps:
        raise ValueError(f"{name} initial/final steps do not match its declared integration")
    if any(b["step"] <= a["step"] for a,b in zip(rows, rows[1:])):
        raise ValueError(f"{name} thermo steps must increase strictly")
    return rows


def _force_summary(state):
    norms = []
    for atom in state["atoms"]:
        force = atom.get("forces_kj_mol_angstrom")
        if not isinstance(force, (list, tuple)) or len(force) != 3:
            raise ValueError("full final atom records require three canonical force components")
        vector = [_number(value, "final force component") for value in force]
        norm = math.hypot(*vector)
        if not math.isfinite(norm):
            raise ValueError("final force vector norm is unrepresentable")
        norms.append(norm)
    maximum = max(norms)
    rms = maximum*math.sqrt(math.fsum((value/maximum)**2 for value in norms)/len(norms)) if maximum else 0.0
    if not math.isfinite(rms):
        raise ValueError("final force RMS is unrepresentable")
    return {"atom_count": len(norms), "max_atom_force_norm_kj_mol_angstrom": maximum,
            "rms_atom_force_norm_kj_mol_angstrom": rms,
            "rms_definition": "sqrt(mean over atoms of Fx^2+Fy^2+Fz^2)); not a component RMS",
            "interpretation": "instantaneous thermal forces; not minimization residuals or an equilibrium certificate"}


def _raw_thermo_consistency(run, nvt, npt, masses, cutoff):
    from ..protocols.alkane import BAR_PER_ATM, _THERMO_COLUMNS
    blocks = run.get("thermo_blocks")
    if not isinstance(blocks, list) or len(blocks) != 3:
        raise ValueError("full raw minimization/NVT/NPT thermo blocks are required")
    for block in blocks:
        if (not isinstance(block, dict) or tuple(block.get("columns", ())) != _THERMO_COLUMNS
                or not isinstance(block.get("rows"), list) or not block["rows"]):
            raise ValueError("raw preparation thermo blocks differ from the prescribed columns")
        for values in block["rows"]:
            if not isinstance(values, list) or len(values) != len(_THERMO_COLUMNS):
                raise ValueError("raw preparation thermo row has the wrong shape")
            for value in values:
                _number(value, "raw thermo value")
            raw = dict(zip(_THERMO_COLUMNS, values))
            lengths = [_number(raw[name], "raw box length", positive=True) for name in ("Lx", "Ly", "Lz")]
            if min(lengths) <= 2*cutoff:
                raise ValueError("raw preparation box violates the selected cutoff regime")
            _close(math.prod(lengths), _number(raw["Volume"], "raw volume", positive=True), "raw box volume")
    for block, rows in zip(blocks[1:], (nvt, npt)):
        if len(block["rows"]) != len(rows):
            raise ValueError("projected preparation stages disagree with their raw thermo row counts")
        for values, row in zip(block["rows"], rows):
            raw = dict(zip(block["columns"], values))
            if raw["Step"] != row["step"]:
                raise ValueError("projected stage step differs from the raw engine thermo")
            _close(raw["Temp"], row["temperature_k"], "projected/raw temperature")
            _close(raw["Press"]*BAR_PER_ATM, row["pressure_bar"], "projected/raw pressure")
            volume = _number(raw["Volume"], "raw thermo volume", positive=True)
            _close(density_from_mass_volume(masses, volume), row["density_g_cm3"], "projected/raw density")
            native = math.fsum(masses)/(0.602214129*volume)
            _close(raw["Density"], native, "pinned native mass density")
    for key in ("temperature_k", "pressure_bar", "density_g_cm3"):
        _close(nvt[-1][key], npt[0][key], f"NVT/NPT boundary {key}")


def _numerical_review(review, assignment, report, protocol):
    geometry = {"box": assignment["box"], "atoms": [
        {key: atom[key] for key in _GEOMETRY_ATOM_KEYS} for atom in assignment["atoms"]]}
    result = {"status": "not_supplied", "starting_geometry_canonical_sha256": _sha(geometry),
              "assignment_canonical_sha256": _sha(assignment),
              "applies_to": "original_preparation_starting_configuration_only",
              "final_configuration_numerics": "not_established; a starting-state review is not transferable through MD",
              "reported_numerical_diagnostics_passed": False, "producer_engine_identity_verified": False}
    if review is None:
        return False, result
    if (not isinstance(review, dict) or review.get("kind") != "periodic_numerical_qualification"
            or review.get("status") != "completed_numerical_review"):
        raise ValueError("numerical_review requires a completed periodic numerical qualification record")
    if review.get("assignment_canonical_sha256") != result["assignment_canonical_sha256"]:
        raise ValueError("numerical review belongs to a different assignment")
    if review.get("geometry_canonical_sha256") != result["starting_geometry_canonical_sha256"]:
        raise ValueError("numerical review belongs to a different starting geometry")
    from ..protocols.numerics import NumericalQualificationSpec, compare_numerical_records
    numerical_spec = NumericalQualificationSpec.from_mapping(review.get("spec"))
    if numerical_spec.to_dict() != review["spec"]:
        raise ValueError("numerical review must retain all predeclared controls")
    if review.get("spec_canonical_sha256") != _sha(numerical_spec.to_dict()):
        raise ValueError("numerical review specification fingerprint changed")
    stages = review.get("stages")
    if not isinstance(stages, list) or not stages:
        raise ValueError("full numerical stage records are required; a summary cannot establish readiness")
    records = []
    for stage in stages:
        if (not isinstance(stage, dict) or not isinstance(stage.get("record"), dict)
                or stage.get("status") != stage["record"].get("status")):
            raise ValueError("full consistent numerical stage records are required")
        records.append(stage["record"])
    derived = compare_numerical_records(records, numerical_spec)
    if not isinstance(review.get("review"), dict) or _canonical(derived) != _canonical(review["review"]):
        raise ValueError("numerical review summary contradicts its actual stage comparisons")
    for key in ("numerical_status", "selected_settings", "assignment_canonical_sha256",
                "geometry_canonical_sha256", "model_id", "atom_count", "engine_identity"):
        if review.get(key) != derived.get(key):
            raise ValueError(f"numerical review {key} contradicts its actual stages")
    result["full_stage_comparison_recomputed"] = True
    numerical_status = review.get("numerical_status")
    if numerical_status not in ("fixed_geometry_diagnostics_passed", "requires_review"):
        raise ValueError("numerical review has an unsupported numerical status")
    result.update(status=numerical_status,
                  reported_numerical_diagnostics_passed=numerical_status == "fixed_geometry_diagnostics_passed",
                  numerical_review_canonical_sha256=_sha(review))
    if numerical_status != "fixed_geometry_diagnostics_passed":
        if review.get("selected_settings") is not None:
            raise ValueError("a failed numerical review cannot select qualified settings")
        return False, result
    expected = {"cutoff_angstrom": protocol.cutoff_angstrom,
                "ewald_alpha_per_angstrom": protocol.ewald_alpha_per_angstrom,
                "reciprocal_nmax": protocol.reciprocal_nmax}
    settings = review.get("selected_settings")
    if not isinstance(settings, dict) or set(settings) != set(expected):
        raise ValueError("numerically reviewed settings differ from the actual preparation protocol")
    normalized = {"cutoff_angstrom": _number(settings["cutoff_angstrom"], "review cutoff", positive=True),
                  "ewald_alpha_per_angstrom": _number(settings["ewald_alpha_per_angstrom"], "review alpha", positive=True),
                  "reciprocal_nmax": _integer(settings["reciprocal_nmax"], "review reciprocal bound", 1)}
    if normalized != expected:
        raise ValueError("numerically reviewed settings differ from the actual preparation protocol")
    result["selected_settings"] = settings
    identity = review.get("engine_identity")
    actual_identity = report["engine_run"].get("engine_identity")
    if identity is not None and actual_identity is not None:
        if not isinstance(identity, dict) or not identity or _canonical(identity) != _canonical(actual_identity):
            raise ValueError("numerical review and preparation producer engine identities differ")
        result["producer_engine_identity_verified"] = True
    else:
        result["producer_engine_identity_status"] = "not_available; reported numerical pass is not adopted"
    return result["producer_engine_identity_verified"], result


def analyze_preparation_readiness(assigned, preparation_record, spec, *, numerical_review=None):
    """Review one full preparation execution record without running an engine.

    This accepts full process records, not public summaries whose atom arrays
    were intentionally omitted.  No outcome establishes phase, equilibrium,
    independent-cell agreement or material-property accuracy.
    """
    from ..protocols.alkane import PeriodicAlkaneProtocolSpec
    from .structure import analyze_chain_structure
    if not isinstance(spec, PreparationReadinessSpec):
        raise ValueError("spec must be a PreparationReadinessSpec")
    assignment = _assignment(assigned)
    report = preparation_record
    if (not isinstance(report, dict) or report.get("kind") != _PREPARATION_KIND
            or report.get("schema_version") != "1" or report.get("model_id") != _MODEL):
        raise ValueError("review requires a full supported preparation execution record")
    try:
        _canonical(report)
    except (ValueError, TypeError, OverflowError, RecursionError) as exc:
        raise ValueError("preparation execution record must contain finite JSON data") from exc
    assignment_sha = _sha(assignment)
    if report.get("assignment_canonical_sha256") != assignment_sha:
        raise ValueError("preparation record belongs to a different original assignment")
    protocol = PeriodicAlkaneProtocolSpec.from_mapping(report.get("protocol_spec"))
    if protocol.to_dict() != report["protocol_spec"]:
        raise ValueError("preparation must retain its complete original protocol controls")
    for actual, expected, name in ((protocol.temperature_k, spec.temperature_k, "temperature target"),
            (protocol.pressure_bar, spec.pressure_bar, "pressure target"),
            (protocol.timestep_fs, spec.timestep_fs, "integration timestep")):
        if actual != expected:
            raise ValueError(f"actual preparation {name} differs from the predeclared review")
    result = {"schema_version": "1", "kind": "preparation_readiness_review",
              "status": "requires_preparation", "assignment_canonical_sha256": assignment_sha,
              "preparation_record_canonical_sha256": _sha(report), "spec": spec.to_dict(),
              "production_ready": False, "equilibrium_status": "not_assessed",
              "density_prediction_validation": "not_performed", "phase_identification": "not_performed",
              "engine_executed_by_analysis": False, "independent_cell_count": 1,
              "gates": {}, "flags": [], "force_summary": None, "structure_snapshot": None,
              "density_sampling_review": None, "numerical_qualification_scope": None,
              "unverified_requirements": ["final-configuration and trajectory numerical sensitivity",
                  "independently prepared cells and preparation-history sensitivity",
                  "phase identity and unobserved slow/chain relaxation modes",
                  "finite-size sensitivity", "state-matched external density comparison"],
              "record_authentication": "consistency review of supplied records; not a cryptographic authentication of engine output"}
    if report.get("status") == "failed":
        reason = report.get("failure_reason")
        if not isinstance(reason, str) or not reason.strip():
            raise ValueError("failed preparation must retain its actual failure reason")
        result.update(failure_phase=report.get("failure_phase"), failure_reason=reason,
                      gates={"engine_execution_completed": False}, flags=["engine_execution_completed"])
        return result
    run = report.get("engine_run")
    if (report.get("status") != "completed_preparation_pilot" or report.get("simulation_performed") is not True
            or not isinstance(run, dict) or run.get("status") != "completed"
            or run.get("simulation_performed") is not True):
        raise ValueError("completed readiness analysis requires an actually completed engine execution")
    unit = report.get("independent_output_unit_check")
    if (not isinstance(unit, dict) or unit.get("confirmed_units_style") != "real"
            or not isinstance(unit.get("header"), str)
            or not re.search(r"(?:^|,)\s*units\s*=\s*real(?:,|$)", unit["header"])):
        raise ValueError("native write_data real-unit confirmation is missing")
    stages = report.get("stages")
    if not isinstance(stages, dict) or set(stages) != {"minimization", "nvt", "npt"}:
        raise ValueError("full minimization/NVT/NPT stage records are required")
    masses = [atom["mass_g_mol"] for atom in assignment["atoms"]]
    nvt = _stage(stages["nvt"], 0, protocol.nvt_steps, protocol.timestep_fs, masses, "NVT")
    npt = _stage(stages["npt"], protocol.nvt_steps, protocol.npt_steps, protocol.timestep_fs, masses, "NPT")
    _raw_thermo_consistency(run, nvt, npt, masses, protocol.cutoff_angstrom)
    if _integer(report.get("md_integration_steps_executed"), "executed MD steps", 1) != protocol.nvt_steps+protocol.npt_steps:
        raise ValueError("actual MD step accounting differs from the two stage records")
    final = report.get("final_state")
    if not isinstance(final, dict) or not isinstance(final.get("atoms"), list) or not final["atoms"]:
        raise ValueError("full final_state.atoms is required; a reduced public summary cannot establish readiness")
    if _integer(final.get("step"), "final step", 1) != protocol.nvt_steps+protocol.npt_steps:
        raise ValueError("final snapshot step differs from the preparation stages")
    structure = analyze_chain_structure(assigned, final)
    forces = _force_summary(final)
    volume = math.prod(final["box"]["lengths_angstrom"])
    _close(_number(final.get("volume_angstrom3"), "final volume", positive=True), volume, "final volume")
    density = density_from_mass_volume(masses, volume)
    _close(_number(final.get("density_g_cm3"), "final density", positive=True), density, "final canonical density")
    for key in ("temperature_k", "pressure_bar", "density_g_cm3"):
        value = _number(final.get(key), f"final {key}", positive=key == "density_g_cm3", nonnegative=key == "temperature_k")
        _close(value, npt[-1][key], f"final/NPT {key}")
    minimum = stages["minimization"]
    if not isinstance(minimum, dict) or type(minimum.get("force_tolerance_met")) is not bool:
        raise ValueError("actual minimizer stopping evidence is required")
    force_met = minimum["force_tolerance_met"]
    if (force_met != (minimum.get("status") == "force_tolerance_met")
            or force_met != (minimum.get("stopping_criterion") == "force tolerance")):
        raise ValueError("minimizer force-convergence flag contradicts its stopping criterion")
    sampling = analyze_density_sampling(npt, spec.sampling_spec)
    numerical_pass, numerical = _numerical_review(numerical_review, assignment, report, protocol)
    maximum, rms = forces["max_atom_force_norm_kj_mol_angstrom"], forces["rms_atom_force_norm_kj_mol_angstrom"]
    force_health = ((spec.max_force_norm_kj_mol_angstrom is None or maximum <= spec.max_force_norm_kj_mol_angstrom)
                    and (spec.max_rms_force_kj_mol_angstrom is None or rms <= spec.max_rms_force_kj_mol_angstrom))
    order = structure.get("nematic_order_parameter")
    orientation_pass = order is not None and _number(order, "nematic order") <= spec.max_nematic_order
    gates = {"engine_execution_completed": True, "native_units_confirmed": True,
             "protocol_steps_and_state_consistent": True,
             "minimizer_policy": force_met or not spec.require_minimizer_force_convergence,
             "finite_final_forces_and_requested_development_bounds": force_health,
             "liquid_candidate_orientation_screen": orientation_pass,
             "npt_density_sampling": sampling["sampling_status"] == "diagnostics_passed",
             "starting_configuration_numerical_review": numerical_pass}
    result.update(gates=gates, flags=[key for key, passed in gates.items() if not passed],
        force_summary={**forces, "requested_development_bounds": {
            "max_force_norm_kj_mol_angstrom": spec.max_force_norm_kj_mol_angstrom,
            "max_rms_force_kj_mol_angstrom": spec.max_rms_force_kj_mol_angstrom},
            "bound_role": "optional gross-instability development heuristics; not universal physical limits"},
        structure_snapshot=structure, density_sampling_review=sampling,
        numerical_qualification_scope=numerical, minimizer=minimum,
        minimizer_force_norms={"status": "not_recorded_by_existing_protocol", "initial": None, "final": None,
                              "note": "stopping-criterion evidence is separate from final thermal forces"},
        md_integration_steps_executed=protocol.nvt_steps+protocol.npt_steps,
        preparation_duration_fs=(protocol.nvt_steps+protocol.npt_steps)*protocol.timestep_fs,
        npt_duration_fs=protocol.npt_steps*protocol.timestep_fs,
        final_state_canonical_sha256=_sha(final))
    if not all(gates[key] for key in ("minimizer_policy", "finite_final_forces_and_requested_development_bounds",
                                      "liquid_candidate_orientation_screen", "starting_configuration_numerical_review")):
        result["status"] = "requires_preparation"
    elif not gates["npt_density_sampling"]:
        result["status"] = "insufficient_sampling"
    else:
        result["status"] = "diagnostics_passed_but_unqualified"
    _canonical(result)
    return result
