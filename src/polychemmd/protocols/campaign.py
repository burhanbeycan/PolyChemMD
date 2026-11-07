"""Declared density campaigns, resource accounting and missing-evidence gates.

Planning is dependency-free and performs no construction, engine discovery or MD.
A candidate phase and a requested physical duration are hypotheses to evaluate.
"""

from dataclasses import dataclass
import hashlib
import json
import math
import re


_MODEL = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"
_TOP = {"schema_version", "chemistry", "phase_target", "carbon_counts", "temperature_k",
        "pressure_bar", "target_density_range_g_cm3", "initialization", "replicates", "integration", "numerics",
        "sampling_gates", "resources", "calibration"}
_OPTIONAL_TOP = {"assignment_representation"}
_LEGACY_REPRESENTATION = "explicit_all_pairs_v1"
_COMPACT_REPRESENTATION = "compact_sparse_pairs_v2"
_COMPACT_MAX_ATOMS = 8192
_MAX_ATOMS_PER_CHAIN = 512
_MAX_INTRAMOLECULAR_PAIRS = 2_000_000


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _keys(value, expected, name):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(f"{name} must contain exactly {sorted(expected)}")


def _number(value, name, minimum=0.0, maximum=1e12):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be representable") from exc
    if not math.isfinite(result) or not minimum < result <= maximum:
        raise ValueError(f"{name} must be in ({minimum}, {maximum}]")
    return result


def _integer(value, name, minimum, maximum):
    if type(value) is not int or not minimum <= value <= maximum:
        raise ValueError(f"{name} must be an integer in [{minimum}, {maximum}]")
    return value


def _label(value, name):
    if not isinstance(value, str) or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]{0,95}", value):
        raise ValueError(f"{name} must be a short ASCII identifier")
    return value



def _derived_positive(operation, name):
    try:
        result = operation()
    except (ArithmeticError, ValueError) as exc:
        raise ValueError(f"{name} is outside the positive finite numeric range") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} is outside the positive finite numeric range")
    return result


def _geometry_values(data):
    counts = data["carbon_counts"]
    mass = sum(count * 12.011 + (2 * count + 2) * 1.008 for count in counts)
    volume = _derived_positive(lambda: mass / (0.602214076 * data["initialization"]["density_g_cm3"]), "initial volume")
    length = _derived_positive(lambda: volume ** (1 / 3), "initial box length")
    cutoff = max(settings["cutoff_angstrom"] for settings in data["numerics"].values())
    density_upper = _derived_positive(lambda: mass / (0.602214076 * (2 * cutoff) ** 3), "density upper bound")
    return mass, volume, length, cutoff, density_upper


def _assignment_limits(data):
    representation = data.get("assignment_representation", _LEGACY_REPRESENTATION)
    if representation not in (_LEGACY_REPRESENTATION, _COMPACT_REPRESENTATION):
        raise ValueError("unsupported assignment_representation")
    compact = representation == _COMPACT_REPRESENTATION
    return representation, (_COMPACT_MAX_ATOMS if compact else 512), (1024 if compact else 64)


def _validate(data):
    if not isinstance(data, dict) or not _TOP.issubset(data) or not set(data).issubset(_TOP | _OPTIONAL_TOP):
        raise ValueError(f"campaign requires {sorted(_TOP)}; optional fields are {sorted(_OPTIONAL_TOP)}")
    representation, max_atoms, max_chains = _assignment_limits(data)
    if data["schema_version"] != "1" or data["chemistry"] != "neutral_explicit_n_alkane":
        raise ValueError("only schema 1 neutral explicit n-alkane planning is supported")
    if data["phase_target"] != "liquid_candidate":
        raise ValueError("this bounded contract only plans an unverified liquid_candidate")
    counts = data["carbon_counts"]
    if not isinstance(counts, list) or not counts or len(counts) > max_chains:
        raise ValueError(f"carbon_counts must contain 1 to {max_chains} finite-chain carbon counts")
    for value in counts:
        _integer(value, "chain carbon count", 2, 170)
    atom_count = sum(3 * value + 2 for value in counts)
    if atom_count > max_atoms:
        raise ValueError(f"campaign exceeds the {max_atoms}-atom selected periodic assignment domain")
    if representation == _COMPACT_REPRESENTATION:
        atoms_per_chain = [3 * value + 2 for value in counts]
        if any(count > _MAX_ATOMS_PER_CHAIN for count in atoms_per_chain):
            raise ValueError("compact assignment permits at most 512 explicit atoms per chain")
        intramolecular_pairs = sum(count * (count - 1) // 2 for count in atoms_per_chain)
        if intramolecular_pairs > _MAX_INTRAMOLECULAR_PAIRS:
            raise ValueError("compact assignment exceeds the 2000000 intramolecular-pair audit budget")
    _number(data["temperature_k"], "temperature_k", maximum=5000)
    _number(data["pressure_bar"], "pressure_bar", maximum=100000)
    density_range = data["target_density_range_g_cm3"]
    if density_range is not None:
        if not isinstance(density_range, list) or len(density_range) != 2:
            raise ValueError("target_density_range_g_cm3 must be null or two planning bounds")
        low, high = (_number(value, "target density bound", maximum=5) for value in density_range)
        if high <= low:
            raise ValueError("target density bounds must be strictly increasing")
    initial = data["initialization"]
    _keys(initial, {"strategy", "density_g_cm3", "min_distance_angstrom", "packing_attempts_per_molecule"}, "initialization")
    allowed_strategies = {"dilute_rigid_then_npt"}
    if representation == _COMPACT_REPRESENTATION:
        allowed_strategies.add("seeded_rigid_then_npt")
    if not isinstance(initial["strategy"], str) or initial["strategy"] not in allowed_strategies:
        raise ValueError(f"initialization strategy must be one of {sorted(allowed_strategies)}")
    _number(initial["density_g_cm3"], "initial density", maximum=5)
    _number(initial["min_distance_angstrom"], "minimum packing distance", maximum=10)
    _integer(initial["packing_attempts_per_molecule"], "packing attempts", 1, 100000)
    reps = data["replicates"]
    if not isinstance(reps, list) or not 3 <= len(reps) <= 16:
        raise ValueError("declare 3 to 16 independently prepared cell recipes")
    identities, seeds = set(), set()
    for rep in reps:
        _keys(rep, {"cell_id", "conformer_seeds", "packing_seed", "velocity_seed"}, "replicate")
        label = _label(rep["cell_id"], "cell_id")
        if label in identities:
            raise ValueError("independent cell IDs must be unique")
        identities.add(label)
        conformers = rep["conformer_seeds"]
        if not isinstance(conformers, list) or len(conformers) != len(counts):
            raise ValueError("one independent conformer seed is required per chain and cell")
        for seed in [*conformers, rep["packing_seed"], rep["velocity_seed"]]:
            _integer(seed, "stochastic seed", 1, 900000000)
            if seed in seeds:
                raise ValueError("construction and velocity seeds must be distinct across primary cell recipes")
            seeds.add(seed)
    integ = data["integration"]
    _keys(integ, {"timestep_fs", "minimize_max_iterations", "minimize_max_evaluations", "nvt_steps",
                   "npt_steps", "production_steps", "segment_steps", "thermo_stride_steps"}, "integration")
    _number(integ["timestep_fs"], "timestep_fs", maximum=1)
    _integer(integ["minimize_max_iterations"], "minimize_max_iterations", 1, 5000)
    _integer(integ["minimize_max_evaluations"], "minimize_max_evaluations", 1, 50000)
    for key in ("nvt_steps", "npt_steps"):
        _integer(integ[key], key, 1, 9999)
    if integ["nvt_steps"] + integ["npt_steps"] > 10000:
        raise ValueError("initial preparation must fit the bounded 10000-step runner")
    _integer(integ["production_steps"], "production_steps", 1, 10000000)
    _integer(integ["segment_steps"], "segment_steps", 1, 10000)
    _integer(integ["thermo_stride_steps"], "thermo_stride_steps", 1, integ["segment_steps"])
    if integ["segment_steps"] % integ["thermo_stride_steps"] or integ["production_steps"] % integ["segment_steps"]:
        raise ValueError("production must contain full segments and segments full thermo strides")
    numerics = data["numerics"]
    _keys(numerics, {"baseline", "cutoff_variant", "ewald_variant"}, "numerics")
    for name, settings in numerics.items():
        _keys(settings, {"cutoff_angstrom", "ewald_alpha_per_angstrom", "reciprocal_nmax"}, name)
        _number(settings["cutoff_angstrom"], "cutoff_angstrom", maximum=100)
        _number(settings["ewald_alpha_per_angstrom"], "ewald_alpha_per_angstrom", maximum=10)
        _integer(settings["reciprocal_nmax"], "reciprocal_nmax", 1, 20)
    baseline, cutoff, ewald = (numerics[name] for name in ("baseline", "cutoff_variant", "ewald_variant"))
    if not (cutoff["cutoff_angstrom"] > baseline["cutoff_angstrom"]
            and cutoff["ewald_alpha_per_angstrom"] == baseline["ewald_alpha_per_angstrom"]
            and cutoff["reciprocal_nmax"] == baseline["reciprocal_nmax"]):
        raise ValueError("cutoff sensitivity must increase only the cutoff")
    if not (ewald["cutoff_angstrom"] == baseline["cutoff_angstrom"]
            and ewald["ewald_alpha_per_angstrom"] == baseline["ewald_alpha_per_angstrom"]
            and ewald["reciprocal_nmax"] > baseline["reciprocal_nmax"]):
        raise ValueError("Ewald sensitivity must increase only reciprocal_nmax")
    _, _, box_length, _, _ = _geometry_values(data)
    if box_length <= 2 * max(settings["cutoff_angstrom"] for settings in numerics.values()):
        raise ValueError("requested initial density gives a cubic box incompatible with sensitivity cutoffs")
    gates = data["sampling_gates"]
    _keys(gates, {"min_independent_cells", "min_production_duration_ps", "density_relative_drift_limit",
                  "density_relative_interval_halfwidth_limit", "temperature_absolute_mean_deviation_k",
                  "pressure_absolute_mean_deviation_bar", "cutoff_relative_density_change_limit",
                  "ewald_relative_density_change_limit", "finite_size_assessment_required"}, "sampling_gates")
    _integer(gates["min_independent_cells"], "min_independent_cells", 3, len(reps))
    for key in set(gates) - {"min_independent_cells", "finite_size_assessment_required"}:
        _number(gates[key], key, maximum=1e9)
    if type(gates["finite_size_assessment_required"]) is not bool or not gates["finite_size_assessment_required"]:
        raise ValueError("finite-size assessment must remain a required unresolved qualification gate")
    for stage in ("nvt_steps", "npt_steps", "production_steps", "segment_steps"):
        _derived_positive(lambda: integ[stage] * integ["timestep_fs"] / 1000, f"{stage} duration")
    if integ["production_steps"] * integ["timestep_fs"] / 1000 < gates["min_production_duration_ps"]:
        raise ValueError("planned production is shorter than its declared minimum duration gate")
    resources = data["resources"]
    _keys(resources, {"hardware_label", "threads", "max_wall_seconds", "max_engine_runs", "per_run_timeout_seconds"}, "resources")
    _label(resources["hardware_label"], "hardware_label")
    _integer(resources["threads"], "threads", 1, 256)
    _number(resources["max_wall_seconds"], "max_wall_seconds", maximum=86400)
    _number(resources["per_run_timeout_seconds"], "per_run_timeout_seconds", maximum=120)
    _integer(resources["max_engine_runs"], "max_engine_runs", 1, 10000)
    if resources["per_run_timeout_seconds"] > resources["max_wall_seconds"]:
        raise ValueError("one run timeout exceeds total wall budget")
    calibration = data["calibration"]
    if calibration is not None:
        _keys(calibration, {"hardware_label", "engine_executable_sha256", "atom_count", "integration_steps",
                           "wall_seconds", "cutoff_angstrom", "reciprocal_nmax"}, "calibration")
        if calibration["hardware_label"] != resources["hardware_label"]:
            raise ValueError("calibration hardware must match the declared hardware")
        if not isinstance(calibration["engine_executable_sha256"], str) or not re.fullmatch(r"[a-f0-9]{64}", calibration["engine_executable_sha256"]):
            raise ValueError("calibration must identify its actual engine executable hash")
        if _integer(calibration["atom_count"], "calibration atom_count", 1, max_atoms) != atom_count:
            raise ValueError("calibration atom count must match the planned primary cell")
        _integer(calibration["integration_steps"], "calibration integration_steps", 1, 10000000)
        _number(calibration["wall_seconds"], "calibration wall_seconds", maximum=86400)
        _derived_positive(lambda: calibration["wall_seconds"] / calibration["integration_steps"], "calibration seconds per step")
        _number(calibration["cutoff_angstrom"], "calibration cutoff_angstrom", maximum=100)
        _integer(calibration["reciprocal_nmax"], "calibration reciprocal_nmax", 1, 20)
        if any(calibration[key] != baseline[key] for key in ("cutoff_angstrom", "reciprocal_nmax")):
            raise ValueError("calibration settings must match the baseline; variants require separate measurements")
    return data


@dataclass(frozen=True)
class DensityCampaignSpec:
    """Immutable, strict JSON contract; constructing a spec launches no work."""

    payload_json: str

    def __post_init__(self):
        if not isinstance(self.payload_json, str):
            raise ValueError("payload_json must be JSON text")
        try:
            payload = json.loads(self.payload_json)
        except (ValueError, TypeError) as exc:
            raise ValueError("invalid campaign JSON") from exc
        _validate(payload)
        object.__setattr__(self, "payload_json", _canonical(payload))

    @classmethod
    def from_mapping(cls, data):
        _validate(data)
        return cls(_canonical(data))

    def to_dict(self):
        return json.loads(self.payload_json)


def plan_density_campaign(spec):
    """Return requested work and unresolved gates; never manufacture observations."""
    if not isinstance(spec, DensityCampaignSpec):
        raise ValueError("spec must be a DensityCampaignSpec")
    data = spec.to_dict()
    representation, max_atoms, _ = _assignment_limits(data)
    counts, integ, resources = data["carbon_counts"], data["integration"], data["resources"]
    atom_count = sum(3 * count + 2 for count in counts)
    mass, volume, box_length, cutoff, density_upper = _geometry_values(data)
    segments = integ["production_steps"] // integ["segment_steps"]
    per_setting_steps = integ["nvt_steps"] + integ["npt_steps"] + integ["production_steps"]
    requests = []
    for rep in data["replicates"]:
        for setting in data["numerics"]:
            requests.append({"cell_id": rep["cell_id"], "numerical_setting": setting,
                             "recipe": rep, "preparation_runs": 1,
                             "continuation_runs": segments,
                             "integration_steps_requested": per_setting_steps,
                             "production_steps_requested": integ["production_steps"],
                             "status": "planned_not_executed"})
    run_count = len(requests) * (1 + segments)
    run_budget_status = "within_declared_run_count" if run_count <= resources["max_engine_runs"] else "exceeds_declared_run_count"
    target_range = data["target_density_range_g_cm3"]
    if target_range is None:
        capacity = "blocked_target_density_hypothesis_not_declared"
    elif target_range[0] >= density_upper:
        capacity = "incompatible_target_range_with_current_cell_and_cutoff_bound"
    elif target_range[1] >= density_upper:
        capacity = "target_range_partly_incompatible_with_cell_and_cutoff_bound"
    else:
        capacity = "box_cutoff_bound_only; geometry_phase_and_sampling_unverified"
    calibration = data["calibration"]
    projection = None
    if calibration is not None:
        # This is a scheduling estimate, never a claimed measurement of future work.
        seconds_per_step = calibration["wall_seconds"] / calibration["integration_steps"]
        baseline_steps = len(data["replicates"]) * per_setting_steps
        projection = {"seconds_per_integration_step": seconds_per_step,
                      "baseline_integration_seconds_estimated": _derived_positive(lambda: baseline_steps * seconds_per_step, "baseline scheduling projection"),
                      "kind": "linear_projection_from_supplied_measurement_not_executed_work",
                      "covers": "baseline_integration_only",
                      "excludes": ["minimization", "per-process startup and checkpoint IO changes", "cutoff/Ewald variants", "finite-size comparison"],
                      "measurement_authenticated_by_planner": False,
                      "whole_campaign_budget_verified": False}
    result = {
        "schema_version": "1", "kind": "density_campaign_plan", "status": "planned",
        "campaign_spec": data, "campaign_canonical_sha256": hashlib.sha256(spec.payload_json.encode()).hexdigest(),
        "model_id": _MODEL, "simulation_performed": False, "execution_supported": False,
        "execution_boundary": "planning only; no campaign dispatcher, finite-size campaign or phase validator implemented",
        "target_state": {"chemistry": data["chemistry"], "phase_target": data["phase_target"],
                         "phase_observed": "not_assessed", "temperature_k": data["temperature_k"], "pressure_bar": data["pressure_bar"],
                         "chain_carbon_counts": counts, "molecule_count_per_cell": len(counts), "atom_count_per_cell": atom_count},
        "initial_state_request": {"density_g_cm3": data["initialization"]["density_g_cm3"],
                                  "density_role": "construction_input_not_prediction", "cubic_box_length_angstrom": box_length,
                                  "volume_angstrom3": volume, "total_molar_mass_g_mol": mass,
                                  "geometry_feasibility": "not_tested; actual conformer extent and packing can fail",
                                  "compression_to_target_phase": "not_established"},
        "state_capacity_gate": {"status": capacity,
                                "target_density_range_g_cm3": target_range,
                                "range_role": "declared_planning_hypothesis_not_measured_or_external_reference",
                                "strict_density_upper_bound_g_cm3": density_upper,
                                "bound_basis": "cell molar mass / (exact-SI Avogadro factor * (2*largest requested cutoff)^3); each box axis must be strictly longer than twice cutoff",
                                "largest_cutoff_angstrom": cutoff,
                                "current_default_atom_capacity": 512,
                                "geometric_bound_is_sufficient_for_liquid_preparation": False,
                                "required_action_if_incompatible": "increase supported cell capacity and verify numerical scaling, or independently validate a different numerical model; do not silently shrink cutoff"},
        "sampling_accounting": {"independent_cell_recipes_requested": len(data["replicates"]),
                                "independent_cells_executed": 0, "correlated_frames_executed": 0,
                                "numerical_variants_are_paired_with_each_primary_cell": True,
                                "seed_uniqueness_proves_statistical_independence": False,
                                "continued_segments_are_new_independent_cells": False,
                                "production_duration_ps_per_setting_per_cell": integ["production_steps"] * integ["timestep_fs"] / 1000,
                                "preparation_duration_ps_per_setting_per_cell": (integ["nvt_steps"] + integ["npt_steps"]) * integ["timestep_fs"] / 1000},
        "work_requests": requests,
        "resource_accounting": {"engine_runs_requested": run_count,
                                "integration_steps_requested": len(requests) * per_setting_steps,
                                "minimization_max_iterations_per_preparation": integ["minimize_max_iterations"],
                                "run_count_budget_status": run_budget_status,
                                "timeout_sum_seconds": run_count * resources["per_run_timeout_seconds"],
                                "timeout_sum_is_elapsed_time_estimate": False,
                                "wall_budget_status": "unqualified_no_complete_campaign_measurement",
                                "projection": projection},
        "qualification_gates": {"observable": "mass_density", "thresholds": data["sampling_gates"],
                                "status": "not_evaluated", "equilibrium_status": "not_assessed",
                                "required_evidence": ["fresh actual cell/topology and parameter hashes", "independent prepared-cell identities, not frame count",
                                                      "uninterrupted same-state NPT sampling with authenticated restart lineage",
                                                      "temperature/pressure state sampling and observed phase assessment",
                                                      "density stationarity, correlation-aware uncertainty and block-size sensitivity",
                                                      "agreement across independent cells", "paired cutoff and reciprocal-bound sensitivity",
                                                      "finite-size comparison", "documented finite resource consumption and retained failures"],
                                "duration_threshold_is_sufficient_for_equilibrium": False},
        "external_comparison": {"status": "deferred_until_sampling_and_state_gates_pass",
                                "required_reference_match": ["chemical_identity", "chain_length_distribution", "temperature_k", "pressure_bar", "phase", "density_definition"],
                                "comparison_performed": False},
        "production_ready": False, "density_prediction_validation": "not_performed",
        "material_property_validation": "not_performed",
    }
    if "assignment_representation" in data:
        atoms_per_chain = [3 * count + 2 for count in counts]
        result["assignment_capacity"] = {
            "assignment_representation": representation,
            "assignment_schema_version": "2" if representation == _COMPACT_REPRESENTATION else "1",
            "engine_nonbonded_rules_representation": "sparse_topological_exceptions" if representation == _COMPACT_REPRESENTATION else "enumerated_pairs",
            "atom_capacity": max_atoms,
            "max_atoms_per_chain": _MAX_ATOMS_PER_CHAIN,
            "intramolecular_pairs_requested": sum(count * (count - 1) // 2 for count in atoms_per_chain),
            "intramolecular_pair_audit_budget": _MAX_INTRAMOLECULAR_PAIRS if representation == _COMPACT_REPRESENTATION else None,
            "physical_state_validation": "not_performed",
        }
        result["state_capacity_gate"]["selected_assignment_atom_capacity"] = max_atoms
        result["sampling_accounting"]["planned_construction_recipe"] = "independently_seeded_rigid_packing_not_executed"
        result["sampling_accounting"]["tiling_one_cell_creates_independent_preparations"] = False
    return result
