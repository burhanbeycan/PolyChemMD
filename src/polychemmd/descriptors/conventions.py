"""Explicit LAMMPS-to-canonical torsion convention conversions.

Only the documented ``fourier`` and ``opls`` functional forms in LAMMPS
``real`` units are supported. These routines implement published equations;
they do not import a force-field database or assign parameters to molecules.

Sources:
https://docs.lammps.org/dihedral_fourier.html
https://docs.lammps.org/dihedral_opls.html
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
import math
from numbers import Real
from typing import Any

from .torsion import FourierTerm, TorsionPotential

KCAL_TO_KJ = 4.184


def _fields(data: Any, expected: set[str], name: str) -> Mapping[str, Any]:
    if not isinstance(data, Mapping):
        raise ValueError(f"{name} must be an object")
    if any(not isinstance(key, str) for key in data):
        raise ValueError(f"{name} field names must be strings")
    unknown, missing = set(data) - expected, expected - set(data)
    if unknown or missing:
        raise ValueError(f"{name}: missing={sorted(missing)}, unknown={sorted(unknown)}")
    return data


def _real(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


def _convert(value: float, factor: float, name: str) -> float:
    result = value * factor
    if not math.isfinite(result) or (value != 0.0 and result == 0.0):
        raise ValueError(f"{name} conversion exceeds the usable floating-point range")
    return result


def _header(data: Mapping[str, Any], style: str) -> None:
    if data["style"] != style:
        raise ValueError(f"style must be {style!r}")
    if data["units"] != "real":
        raise ValueError("only explicitly specified LAMMPS real units are supported")


def import_lammps_fourier(data: Mapping[str, Any]) -> TorsionPotential:
    """Convert an explicit LAMMPS Fourier record to kJ/mol and radians.

    Schema: ``style='fourier', units='real', m=N, terms=[{k_kcal_mol,
    periodicity, phase_deg}, ...]``. The list length must equal ``m``. Signed
    coefficients, zero coefficients, term order and unwrapped phases are retained.

    LAMMPS permits periodicity zero. The current canonical core requires positive
    periodicity, so zero-periodicity constant-offset terms are rejected explicitly
    even if their coefficient is zero; they are never silently discarded.
    """
    data = _fields(data, {"style", "units", "m", "terms"}, "LAMMPS Fourier record")
    _header(data, "fourier")
    count = data["m"]
    if type(count) is not int or count < 1:
        raise ValueError("m must be a positive integer")
    raw_terms = data["terms"]
    if isinstance(raw_terms, (str, bytes)) or not isinstance(raw_terms, Sequence):
        raise ValueError("terms must be an array")
    if len(raw_terms) != count:
        raise ValueError("m must equal the number of terms")
    terms = []
    for raw in raw_terms:
        raw = _fields(raw, {"k_kcal_mol", "periodicity", "phase_deg"}, "LAMMPS Fourier term")
        periodicity = raw["periodicity"]
        if type(periodicity) is not int or periodicity < 0:
            raise ValueError("periodicity must be an integer; booleans are not accepted")
        if periodicity == 0:
            raise ValueError("periodicity=0 constant-offset terms are unsupported by the canonical core; do not discard them")
        coefficient = _convert(_real(raw["k_kcal_mol"], "k_kcal_mol"), KCAL_TO_KJ, "energy")
        phase = _convert(_real(raw["phase_deg"], "phase_deg"), math.pi / 180.0, "phase")
        terms.append(FourierTerm(coefficient, periodicity, phase))
    return TorsionPotential(tuple(terms))


def export_lammps_fourier(potential: TorsionPotential) -> dict[str, Any]:
    """Export every canonical term to the explicit LAMMPS real-unit schema.

    No phase normalization or coefficient combining is performed. Unrepresentable
    converted values cause an error rather than silently changing the profile.
    """
    if not isinstance(potential, TorsionPotential):
        raise ValueError("potential must be a TorsionPotential")
    return {
        "style": "fourier",
        "units": "real",
        "m": len(potential.terms),
        "terms": [
            {"k_kcal_mol": _convert(term.k_kj_mol, 1.0 / KCAL_TO_KJ, "energy"),
             "periodicity": term.periodicity,
             "phase_deg": _convert(term.phase_rad, 180.0 / math.pi, "phase")}
            for term in potential.terms
        ],
    }


def import_lammps_opls(data: Mapping[str, Any]) -> TorsionPotential:
    """Convert the four LAMMPS OPLS K values, preserving the explicit 1/2.

    Schema: ``style='opls', units='real', k_kcal_mol=[K1,K2,K3,K4]``.
    The canonical phases are 0, pi, 0, pi for periodicities 1, 2, 3, 4.
    All four terms are retained, including zero coefficients. This function
    converts this single torsion expression, not an entire OPLS force field.
    """
    data = _fields(data, {"style", "units", "k_kcal_mol"}, "LAMMPS OPLS record")
    _header(data, "opls")
    coefficients = data["k_kcal_mol"]
    if (isinstance(coefficients, (str, bytes)) or not isinstance(coefficients, Sequence)
            or len(coefficients) != 4):
        raise ValueError("k_kcal_mol must contain exactly four coefficients")
    terms = []
    for periodicity, raw in enumerate(coefficients, start=1):
        coefficient = _convert(_real(raw, f"K{periodicity}"), KCAL_TO_KJ / 2.0, "energy/prefactor")
        phase = 0.0 if periodicity % 2 else math.pi
        terms.append(FourierTerm(coefficient, periodicity, phase))
    return TorsionPotential(tuple(terms))
