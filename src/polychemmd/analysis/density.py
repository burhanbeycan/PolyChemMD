"""Density conversion and transparent block diagnostics, not equilibration proof.

This module runs no engine. A supplied density series must use g/cm^3, uniform
sampling, one state/ensemble and one uninterrupted production window. The caller
is responsible for recording temperature, pressure, preparation and run identity.
Independent block means are an assumption, not a consequence of choosing a block
size. Chain relaxation and agreement across independently prepared cells require
separate evidence.
"""

from __future__ import annotations

import math
from numbers import Real
import statistics
from typing import Any, Iterable

AVOGADRO_CONSTANT_MOL_INV = 6.02214076e23
ANGSTROM3_TO_CM3 = 1.0e-24
NORMAL_975_QUANTILE = 1.959963984540054


def _positive_finite(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a positive finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be representable as a finite number") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return result


def _values(values: Iterable[float], name: str) -> tuple[float, ...]:
    if isinstance(values, (str, bytes, dict)):
        raise ValueError(f"{name} must be a non-empty iterable of numbers")
    try:
        result = tuple(_positive_finite(value, name) for value in values)
    except TypeError as exc:
        raise ValueError(f"{name} must be a non-empty iterable of numbers") from exc
    if not result:
        raise ValueError(f"{name} must not be empty")
    return result


def _integer(value: Any, name: str, minimum: int) -> int:
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _mean(values: Iterable[float]) -> float:
    # statistics.mean avoids an overflowing intermediate sum for a finite mean.
    result = float(statistics.mean(values))
    if not math.isfinite(result):
        raise ValueError("density statistic is outside the finite numeric range")
    return result


def density_from_mass_volume(masses_g_mol: Iterable[float], volume_angstrom3: float) -> float:
    """Return geometric mass density in g/cm^3 from atomic molar masses.

    ``rho = sum(m_i) / (N_A * V * 1e-24)``. This conversion says nothing about
    equilibration, material identity or force-field validity.
    """
    masses = _values(masses_g_mol, "atomic molar mass")
    volume = _positive_finite(volume_angstrom3, "volume_angstrom3")
    try:
        total_mass = math.fsum(masses)
        denominator = (AVOGADRO_CONSTANT_MOL_INV * ANGSTROM3_TO_CM3) * volume
        if denominator == 0:
            raise ValueError("volume is too small for a resolvable density conversion")
        result = total_mass / denominator
    except OverflowError as exc:
        raise ValueError("density conversion exceeds the finite numeric range") from exc
    if not math.isfinite(result) or result <= 0:
        raise ValueError("density conversion is outside the positive finite numeric range")
    return result


def _lag_one_correlation(values: tuple[float, ...]) -> float | None:
    """Adjacent-block Pearson correlation; a diagnostic, not effective sample size."""
    left, right = values[:-1], values[1:]
    left_mean, right_mean = _mean(left), _mean(right)
    dx = tuple(x - left_mean for x in left)
    dy = tuple(y - right_mean for y in right)
    scale = max(*(abs(x) for x in dx), *(abs(y) for y in dy))
    if scale == 0:
        return None
    x, y = tuple(v / scale for v in dx), tuple(v / scale for v in dy)
    denominator = math.sqrt(math.fsum(v*v for v in x) * math.fsum(v*v for v in y))
    if denominator == 0:
        return None
    return max(-1.0, min(1.0, math.fsum(a*b for a, b in zip(x, y)) / denominator))


def analyze_density(
    densities_g_cm3: Iterable[float],
    block_size: int,
    discard: int = 0,
    min_blocks: int = 8,
    drift_relative_limit: float = 0.02,
) -> dict[str, Any]:
    """Compute density block diagnostics under explicitly unverified assumptions.

    Discard the requested prefix, then omit an incomplete trailing block. At least
    ``min_blocks`` complete blocks are required; the minimum configurable count
    is five. The default relative split-half drift limit is a transparent review
    heuristic, not a validated physical convergence threshold.

    The normal-approximation interval assumes independent block means and is
    reported as conditional, never as a calibrated material-property interval.
    A constant series has no uncertainty estimate and is marked inconclusive.
    No result is labeled equilibrated or a validated polymer density.
    """
    densities = _values(densities_g_cm3, "density_g_cm3")
    size = _integer(block_size, "block_size", 1)
    prefix = _integer(discard, "discard", 0)
    minimum_blocks = _integer(min_blocks, "min_blocks", 5)
    drift_limit = _positive_finite(drift_relative_limit, "drift_relative_limit")
    remaining = len(densities) - prefix
    if remaining <= 0:
        raise ValueError("discard must leave a production window")
    count = remaining // size
    if count < minimum_blocks:
        raise ValueError(f"insufficient complete blocks: {count}; require {minimum_blocks}")
    used_count = count * size
    used = densities[prefix:prefix + used_count]
    block_means = tuple(_mean(used[index:index + size]) for index in range(0, used_count, size))
    mean = _mean(block_means)
    try:
        sd = statistics.stdev(block_means)
    except (OverflowError, ValueError) as exc:
        raise ValueError("block variability is outside the finite numeric range") from exc
    if not math.isfinite(sd):
        raise ValueError("block variability is outside the finite numeric range")
    observed_variance = max(used) > min(used)
    variable_blocks = sd > 0
    standard_error = sd / math.sqrt(count) if variable_blocks else None
    if standard_error is not None and standard_error <= 0:
        raise ValueError("block uncertainty is below the resolvable numeric range")
    split = count // 2
    first_mean, second_mean = _mean(block_means[:split]), _mean(block_means[split:])
    signed_drift = second_mean - first_mean
    relative_drift = abs(signed_drift) / mean
    interval = None
    if standard_error is not None:
        half_width = NORMAL_975_QUANTILE * standard_error
        if not math.isfinite(half_width):
            raise ValueError("conditional interval is outside the finite numeric range")
        interval = {
            "lower_g_cm3": mean - half_width,
            "upper_g_cm3": mean + half_width,
            "nominal_level": 0.95,
            "method": "normal_approximation_assuming_independent_block_means",
            "calibrated": False,
        }
    coarsened = tuple(_mean(block_means[index:index + 2]) for index in range(0, count - 1, 2))
    try:
        coarse_se = statistics.stdev(coarsened) / math.sqrt(len(coarsened)) if len(coarsened) >= 3 else None
    except (OverflowError, ValueError) as exc:
        raise ValueError("coarsened-block statistic is outside the finite numeric range") from exc
    if coarse_se is not None and not math.isfinite(coarse_se):
        raise ValueError("coarsened-block statistic is outside the finite numeric range")
    flags = []
    if not observed_variance:
        flags.append("constant_density_series_does_not_establish_sampling_or_chain_relaxation")
    elif not variable_blocks:
        flags.append("identical_block_means_do_not_establish_independent_sampling")
    if relative_drift > drift_limit:
        flags.append("split_half_drift_exceeds_configured_review_limit")
    status = "inconclusive" if not variable_blocks else ("requires_review" if flags else "diagnostics_only")
    result = {
        "analysis_schema_version": "0.1",
        "observable": "mass_density",
        "units": "g/cm^3",
        "result_kind": "density_time_series_diagnostics",
        "sampling_status": status,
        "equilibrium_status": "not_assessed",
        "polymer_density_validation": "not_performed",
        "simulation_executed_by_analysis": False,
        "n_input": len(densities),
        "n_discarded_prefix": prefix,
        "n_used": used_count,
        "n_omitted_trailing": remaining - used_count,
        "block_size": size,
        "n_blocks": count,
        "mean_g_cm3": mean,
        "observed_block_standard_deviation_g_cm3": sd,
        "conditional_standard_error_g_cm3": standard_error,
        "conditional_normal_interval": interval,
        "block_means_g_cm3": list(block_means),
        "split_half": {
            "first_mean_g_cm3": first_mean,
            "second_mean_g_cm3": second_mean,
            "signed_difference_g_cm3": signed_drift,
            "relative_difference": relative_drift,
            "configured_review_limit": drift_limit,
            "limit_role": "review_heuristic_not_equilibration_proof",
        },
        "lag_one_block_correlation": _lag_one_correlation(block_means),
        "block_length_sensitivity": {
            "coarsened_block_size": 2 * size,
            "n_coarsened_blocks": len(coarsened),
            "coarsened_standard_error_g_cm3": coarse_se,
            "coarsened_to_original_standard_error_ratio": coarse_se / standard_error if coarse_se is not None and standard_error else None,
            "interpretation": "sensitivity diagnostic only; effective independence is not established",
        },
        "flags": flags,
        "assumptions_not_verified": [
            "uniform sampling interval and one state/ensemble",
            "user-selected discard removes preparation transients",
            "block duration exceeds relevant density correlation time",
            "independently prepared cells agree",
            "chain relaxation and force-field applicability",
        ],
    }
    pending = [result]
    while pending:
        value = pending.pop()
        if isinstance(value, dict):
            pending.extend(value.values())
        elif isinstance(value, list):
            pending.extend(value)
        elif isinstance(value, float) and not math.isfinite(value):
            raise ValueError("density diagnostic exceeds the finite numeric range")
    return result
