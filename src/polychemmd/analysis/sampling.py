"""Predeclared density sampling review gates; never an equilibrium certificate.

The implementation uses explicit block means and an initial-positive-pair
estimate of statistical inefficiency. Such finite-series diagnostics can miss
slow modes. Passing them is deliberately weaker than numerical validation.
"""

from dataclasses import asdict, dataclass, fields
import math
import statistics

from .density import NORMAL_975_QUANTILE

_MAX_ROWS = 100000
_MAX_CORRELATION_WORK = 5000000
_ROW_FIELDS = {"step", "temperature_k", "pressure_bar", "density_g_cm3"}
_GATE_NAMES = {"physical_duration", "complete_blocks", "coarsened_blocks",
               "nonconstant_density_and_blocks", "autocorrelation_window_closed",
               "effective_samples", "block_correlation_time_separation", "block_lag_one",
               "block_length_sensitivity", "density_relative_precision",
               "density_split_half_drift", "density_four_window_drift",
               "temperature_sample_mean", "pressure_sample_mean"}


def _number(value, name, *, positive=False, nonnegative=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite real number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is not representable") from exc
    if not math.isfinite(value) or positive and value <= 0 or nonnegative and value < 0:
        raise ValueError(f"{name} is outside its finite allowed range")
    return value


def _integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name} must be an integer in [{low}, {high}]")
    return value


def _mean(values):
    result = float(statistics.mean(values))
    if not math.isfinite(result):
        raise ValueError("sampling statistic is outside its finite numeric range")
    return result


def _sd(values):
    try:
        result = statistics.stdev(values)
    except (ValueError, OverflowError) as exc:
        raise ValueError("sampling variability is outside its finite numeric range") from exc
    if not math.isfinite(result):
        raise ValueError("sampling variability is outside its finite numeric range")
    return result


@dataclass(frozen=True)
class DensitySamplingSpec:
    timestep_fs: float = 0.25
    discard_steps: int = 0
    min_duration_ps: float = 100.0
    block_size: int = 50
    min_blocks: int = 16
    max_autocorrelation_lag: int = 128
    min_effective_samples: float = 50.0
    min_block_inefficiency_ratio: float = 5.0
    max_abs_block_lag_one: float = 0.20
    max_coarsened_se_ratio: float = 1.25
    max_density_relative_half_width: float = 0.02
    max_density_relative_drift: float = 0.02
    target_temperature_k: float = 300.0
    temperature_mean_tolerance_k: float = 15.0
    target_pressure_bar: float = 1.0
    pressure_mean_tolerance_bar: float = 100.0

    def __post_init__(self):
        for name in ("timestep_fs", "min_duration_ps", "min_effective_samples",
                     "min_block_inefficiency_ratio", "max_coarsened_se_ratio",
                     "max_density_relative_half_width", "max_density_relative_drift",
                     "target_temperature_k", "temperature_mean_tolerance_k",
                     "pressure_mean_tolerance_bar"):
            object.__setattr__(self, name, _number(getattr(self, name), name, positive=True))
        object.__setattr__(self, "target_pressure_bar", _number(self.target_pressure_bar, "target_pressure_bar"))
        object.__setattr__(self, "max_abs_block_lag_one", _number(self.max_abs_block_lag_one, "max_abs_block_lag_one", positive=True))
        if self.max_abs_block_lag_one >= 1:
            raise ValueError("max_abs_block_lag_one must be below one")
        if self.max_coarsened_se_ratio <= 1:
            raise ValueError("max_coarsened_se_ratio must exceed one")
        _integer(self.discard_steps, "discard_steps", 0, 10**12)
        _integer(self.block_size, "block_size", 1, _MAX_ROWS)
        _integer(self.min_blocks, "min_blocks", 16, _MAX_ROWS)
        _integer(self.max_autocorrelation_lag, "max_autocorrelation_lag", 4, 4096)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_mapping(cls, data):
        names = {field.name for field in fields(cls)}
        if not isinstance(data, dict) or not set(data).issubset(names):
            raise ValueError(f"sampling fields must be drawn from {sorted(names)}")
        return cls(**data)


def _rows(rows):
    if not isinstance(rows, (list, tuple)) or not rows or len(rows) > _MAX_ROWS:
        raise ValueError(f"rows must be a nonempty list/tuple with at most {_MAX_ROWS} entries")
    result = []
    for row in rows:
        if not isinstance(row, dict) or set(row) != _ROW_FIELDS:
            raise ValueError(f"each row must contain exactly {sorted(_ROW_FIELDS)}")
        step = _integer(row["step"], "step", 0, 10**12)
        result.append({"step": step,
                       "temperature_k": _number(row["temperature_k"], "temperature_k", nonnegative=True),
                       "pressure_bar": _number(row["pressure_bar"], "pressure_bar"),
                       "density_g_cm3": _number(row["density_g_cm3"], "density_g_cm3", positive=True)})
    stride = result[1]["step"]-result[0]["step"] if len(result) > 1 else None
    if stride is not None and (stride <= 0 or any(b["step"]-a["step"] != stride for a,b in zip(result,result[1:]))):
        raise ValueError("sampling rows must have strictly increasing uniformly spaced steps")
    return result, stride


def _correlation(values, lag):
    # A common denominator is required for the finite-series autocovariance.
    mean = _mean(values)
    scale = max(abs(value-mean) for value in values)
    if scale == 0:
        return None
    centered = [(value-mean)/scale for value in values]
    denominator = math.fsum(value*value for value in centered)
    return math.fsum(a*b for a,b in zip(centered[:-lag],centered[lag:]))/denominator


def _inefficiency(values, requested_lag):
    n = len(values)
    limit = min(requested_lag, n-1)
    if n*limit > _MAX_CORRELATION_WORK:
        raise ValueError("autocorrelation work budget exceeded; reduce rows or requested lag")
    if max(values) == min(values):
        return {"status":"constant_series", "statistical_inefficiency":None,
                "effective_samples":None, "window_closed":False, "last_lag":0,
                "lag_correlations":[]}
    # Geyer adjacent autocovariance pairs start at lags (0,1), then
    # (2,3), ... . Stop at the first nonpositive pair; enforce decreasing
    # positive pairs and clamp g >= 1. The finite-series diagnostic has no
    # general equilibration guarantee for an MD observable.
    correlations, sums, previous = [1.0], [], math.inf
    closed = False
    for lag in range(0,limit,2):
        left = 1.0 if lag == 0 else _correlation(values,lag)
        right = _correlation(values,lag+1)
        if lag:
            correlations.append(left)
        correlations.append(right)
        pair_sum = left+right
        if pair_sum <= 0:
            closed = True
            break
        previous = min(previous,pair_sum)
        sums.append(previous)
    g = max(1.0,-1.0+2.0*math.fsum(sums))
    return {"status":"estimated" if closed else "positive_window_not_closed",
            "statistical_inefficiency":g, "effective_samples":n/g,
            "window_closed":closed, "last_lag":len(correlations)-1,
            "lag_correlations":correlations,
            "estimator":"Geyer initial monotone positive adjacent lag pairs starting at 0,1; g>=1",
            "finite_series_estimate_not_independence_proof":True}


def analyze_density_sampling(rows, spec):
    """Review one uninterrupted state/ensemble window with declared thresholds.

    ``discard_steps`` is relative to the first supplied step. All returned
    uncertainty is conditional on the stated sampling assumptions. No caller
    can request a thermodynamic-equilibrium or material-validation certificate.
    """
    if not isinstance(spec,DensitySamplingSpec):
        raise ValueError("spec must be a DensitySamplingSpec")
    supplied, stride = _rows(rows)
    first_step = supplied[0]["step"]
    kept = [row for row in supplied if row["step"]-first_step >= spec.discard_steps]
    duration = ((kept[-1]["step"]-kept[0]["step"])*spec.timestep_fs/1000) if len(kept)>1 else 0.0
    if not math.isfinite(duration):
        raise ValueError("physical sampling duration is not representable")
    result = {
        "schema_version":"1", "kind":"density_sampling_review",
        "sampling_status":"insufficient_sampling", "status":"insufficient_sampling",
        "equilibrium_status":"not_assessed", "production_ready":False,
        "density_prediction_validation":"not_performed",
        "simulation_executed_by_analysis":False, "spec":spec.to_dict(),
        "n_input_frames":len(supplied), "n_discarded_frames":len(supplied)-len(kept),
        "n_retained_frames":len(kept), "independent_cell_count":1,
        "sampling_stride_steps":stride,
        "sampling_interval_fs":stride*spec.timestep_fs if stride is not None else None,
        "retained_duration_ps":duration, "flags":[], "gates":{},
        "uncertainty_interpretation":"conditional finite-series estimates; neither a calibrated confidence interval nor independent-cell uncertainty",
        "unverified_requirements":["discard removes preparation transients", "unobserved slow modes and chain relaxation",
                                   "phase and force-field applicability", "independently prepared cell agreement",
                                   "cell size and numerical sensitivity", "state-matched external comparison"],
    }
    gates = result["gates"]
    gates["physical_duration"] = duration >= spec.min_duration_ps
    count = len(kept)//spec.block_size
    gates["complete_blocks"] = count >= spec.min_blocks
    gates["coarsened_blocks"] = count//2 >= 8
    used_n = count*spec.block_size
    used_duration = (kept[used_n-1]["step"]-kept[0]["step"])*spec.timestep_fs/1000 if used_n>1 else 0.0
    gates["physical_duration"] = used_duration >= spec.min_duration_ps
    result.update(n_complete_blocks=count, n_used_frames=used_n, used_duration_ps=used_duration,
                  n_omitted_trailing_frames=len(kept)-used_n)
    if used_n < 4 or count < 2:
        result["flags"] = [name for name,passed in gates.items() if not passed]+["too_few_frames_for_variability_diagnostics"]
        return result
    used = kept[:used_n]
    density = [row["density_g_cm3"] for row in used]
    mean = _mean(density)
    means = [_mean(density[i:i+spec.block_size]) for i in range(0,used_n,spec.block_size)]
    se = _sd(means)/math.sqrt(count)
    coarse = [_mean(means[i:i+2]) for i in range(0,count-1,2)]
    coarse_se = _sd(coarse)/math.sqrt(len(coarse)) if len(coarse)>1 else None
    autocorrelation = _inefficiency(density,spec.max_autocorrelation_lag)
    g = autocorrelation["statistical_inefficiency"]
    # Preserve failures at either direction of the block-size sensitivity test.
    ratio = coarse_se/se if coarse_se is not None and se > 0 else None
    block_lag = _correlation(means,1)
    autocorrelation_se = _sd(density)*math.sqrt(g/used_n) if g is not None else None
    half_width = NORMAL_975_QUANTILE*max(se,coarse_se or 0.0,autocorrelation_se or 0.0) if se > 0 else None
    split = count//2
    split_drift = abs(_mean(means[split:])-_mean(means[:split]))/mean
    # Four consecutive windows can reveal an excursion hidden by split halves.
    quarter_means = [_mean(density[i*used_n//4:(i+1)*used_n//4]) for i in range(4)]
    window_drift = (max(quarter_means)-min(quarter_means))/mean
    temperature = _mean([row["temperature_k"] for row in used])
    pressure = _mean([row["pressure_bar"] for row in used])
    gates.update({
        "nonconstant_density_and_blocks":max(density)>min(density) and se>0,
        "autocorrelation_window_closed":autocorrelation["window_closed"],
        "effective_samples":g is not None and autocorrelation["effective_samples"]>=spec.min_effective_samples,
        "block_correlation_time_separation":g is not None and spec.block_size>=spec.min_block_inefficiency_ratio*g,
        "block_lag_one":block_lag is not None and abs(block_lag)<=spec.max_abs_block_lag_one,
        "block_length_sensitivity":ratio is not None and 1/spec.max_coarsened_se_ratio<=ratio<=spec.max_coarsened_se_ratio,
        "density_relative_precision":half_width is not None and half_width/mean<=spec.max_density_relative_half_width,
        "density_split_half_drift":split_drift<=spec.max_density_relative_drift,
        "density_four_window_drift":window_drift<=spec.max_density_relative_drift,
        "temperature_sample_mean":abs(temperature-spec.target_temperature_k)<=spec.temperature_mean_tolerance_k,
        "pressure_sample_mean":abs(pressure-spec.target_pressure_bar)<=spec.pressure_mean_tolerance_bar,
    })
    result.update(mean_density_g_cm3=mean, block_means_g_cm3=means,
                  conditional_block_standard_error_g_cm3=se if se>0 else None,
                  coarsened_standard_error_g_cm3=coarse_se if coarse_se and coarse_se>0 else None,
                  coarsened_to_original_se_ratio=ratio, lag_one_block_correlation=block_lag,
                  conditional_autocorrelation_standard_error_g_cm3=autocorrelation_se,
                  conditional_normal_half_width_g_cm3=half_width,
                  conditional_normal_nominal_level=0.95,
                  density_relative_split_half_drift=split_drift,
                  density_relative_four_window_drift=window_drift,
                  mean_temperature_k=temperature, mean_pressure_bar=pressure,
                  autocorrelation=autocorrelation)
    result["flags"] = [name for name,passed in gates.items() if not passed]
    if not gates["physical_duration"] or not gates["complete_blocks"] or not gates["coarsened_blocks"]:
        status = "insufficient_sampling"
    elif not gates["nonconstant_density_and_blocks"]:
        status = "inconclusive"
    elif result["flags"]:
        status = "requires_review"
    else:
        status = "diagnostics_passed"
    result["status"] = result["sampling_status"] = status
    _finite_output(result)
    return result


def _finite_output(value):
    if isinstance(value,dict):
        for item in value.values():
            _finite_output(item)
    elif isinstance(value,(list,tuple)):
        for item in value:
            _finite_output(item)
    elif isinstance(value,float) and not math.isfinite(value):
        raise ValueError("sampling statistic is outside its finite numeric range")


def summarize_density_replicates(records, *, expected_cell_ids):
    """Retain every declared cell; never count frames/restarts as replicates.

    Cell identities and independence are declared provenance, not authenticated
    by this numerical function. All records must concern the same state and
    protocol. Failed/missing sampling cannot be silently dropped to report an
    apparently successful aggregate.
    """
    if (not isinstance(expected_cell_ids,(list,tuple)) or not expected_cell_ids
            or any(not isinstance(v,str) or not v.strip() for v in expected_cell_ids)
            or len(set(expected_cell_ids))!=len(expected_cell_ids)):
        raise ValueError("expected_cell_ids must be distinct nonempty strings")
    if not isinstance(records,(list,tuple)) or len(records)!=len(expected_cell_ids):
        raise ValueError("every declared independent cell, including failures, must be retained")
    names = {"independent_cell_id","state_id","protocol_id","status","sampling_analysis"}
    ids, states, protocols, specs, outcomes, values = set(), set(), set(), [], [], []
    for record in records:
        if not isinstance(record,dict) or set(record)!=names:
            raise ValueError(f"replicate records require exactly {sorted(names)}")
        for name in ("independent_cell_id","state_id","protocol_id"):
            if not isinstance(record[name],str) or not record[name].strip():
                raise ValueError(f"{name} must be a nonempty string")
        identity = record["independent_cell_id"]
        if identity in ids:
            raise ValueError("repeated frames or continuations cannot reuse an independent cell identity")
        ids.add(identity); states.add(record["state_id"]); protocols.add(record["protocol_id"])
        if record["status"] not in ("completed","failed","inconclusive"):
            raise ValueError("replicate status must be completed, failed or inconclusive")
        analysis = record["sampling_analysis"]
        if record["status"] == "completed":
            if (not isinstance(analysis,dict) or analysis.get("kind")!="density_sampling_review"
                    or analysis.get("schema_version")!="1" or not isinstance(analysis.get("spec"),dict)):
                raise ValueError("completed cells require a density_sampling_review")
            spec = DensitySamplingSpec.from_mapping(analysis["spec"]).to_dict()
            specs.append(spec)
            sampling_status = analysis.get("sampling_status")
            if (sampling_status not in ("diagnostics_passed", "insufficient_sampling", "requires_review", "inconclusive")
                    or analysis.get("status") != sampling_status
                    or analysis.get("equilibrium_status") != "not_assessed"
                    or analysis.get("production_ready") is not False
                    or analysis.get("density_prediction_validation") != "not_performed"):
                raise ValueError("sampling review status or evidence boundary is inconsistent")
            if sampling_status == "diagnostics_passed":
                gates = analysis.get("gates")
                if (not isinstance(gates,dict) or set(gates)!=_GATE_NAMES
                        or any(value is not True for value in gates.values())
                        or analysis.get("flags") != []):
                    raise ValueError("diagnostics_passed requires every predeclared review gate to pass")
                n_used = _integer(analysis.get("n_used_frames"), "n_used_frames", 1, _MAX_ROWS)
                if n_used < spec["block_size"]*spec["min_blocks"]:
                    raise ValueError("diagnostics_passed contradicts its complete-block budget")
                if _number(analysis.get("used_duration_ps"), "used_duration_ps", positive=True) < spec["min_duration_ps"]:
                    raise ValueError("diagnostics_passed contradicts its physical duration budget")
                if (abs(_number(analysis.get("mean_temperature_k"),"mean temperature")-spec["target_temperature_k"]) > spec["temperature_mean_tolerance_k"]
                        or abs(_number(analysis.get("mean_pressure_bar"),"mean pressure")-spec["target_pressure_bar"]) > spec["pressure_mean_tolerance_bar"]):
                    raise ValueError("diagnostics_passed contradicts its sample-mean state checks")
                values.append(_number(analysis.get("mean_density_g_cm3"),"cell mean density",positive=True))
        elif analysis is not None:
            raise ValueError("failed/inconclusive engine records must retain no successful sampling analysis")
        outcomes.append({"independent_cell_id":identity,"execution_status":record["status"],
                         "sampling_status":analysis.get("sampling_status") if isinstance(analysis,dict) else "not_available"})
    if ids != set(expected_cell_ids):
        raise ValueError("replicate identities do not match the predeclared cells")
    if len(states)!=1 or len(protocols)!=1 or any(spec!=specs[0] for spec in specs[1:]):
        raise ValueError("replicates must use identical state, protocol and sampling specification")
    result = {"schema_version":"1","kind":"density_independent_cell_review",
              "state_id":next(iter(states)),"protocol_id":next(iter(protocols)),
              "status":"incomplete_or_inconclusive", "n_declared_independent_cells":len(records),
              "n_diagnostics_passed_cells":len(values), "cell_outcomes":outcomes,
              "independence_status":"declared_not_verified", "frames_are_independent_cells":False,
              "equilibrium_status":"not_assessed","density_prediction_validation":"not_performed",
              "mean_density_g_cm3":None,"between_cell_standard_deviation_g_cm3":None,
              "conditional_between_cell_standard_error_g_cm3":None}
    if len(values)==len(records) and len(values)>=3 and max(values)>min(values):
        sd = _sd(values)
        result.update(status="replicate_diagnostics_passed",mean_density_g_cm3=_mean(values),
                      between_cell_standard_deviation_g_cm3=sd,
                      conditional_between_cell_standard_error_g_cm3=sd/math.sqrt(len(values)))
    _finite_output(result)
    return result
