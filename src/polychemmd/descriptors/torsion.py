"""A complete, explicitly defined Fourier torsion representation.

The potential is exactly ``sum(k * (1 + cos(n * phi - phase)))`` in kJ/mol,
with all angles in radians. There is no implicit one-half prefactor. Imported
force-field terms must first be converted to this convention. Signed finite
coefficients are permitted and retained; no coefficient or phase is discarded.

Grid Boltzmann probabilities describe an isolated one-dimensional model on a
uniform angular grid, not a polymer's condensed-phase torsion populations. The
discrete Shannon entropy depends on grid size and is not thermodynamic entropy.
This module uses only the Python standard library.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
import math
from numbers import Real
from typing import Any

GAS_CONSTANT_KJ_MOL_K = 0.00831446261815324
CONVENTION = "sum_k_1_plus_cos_nphi_minus_phase"


def _finite_real(value: Any, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number, not a boolean")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result


@dataclass(frozen=True)
class FourierTerm:
    """One term of the module's explicit Fourier convention.

    ``periodicity`` must be a positive Python integer. ``k_kj_mol`` is a
    coefficient rather than a non-negative barrier height: signed values are
    valid under the stated mathematical convention.
    """

    k_kj_mol: float
    periodicity: int
    phase_rad: float = 0.0

    def __post_init__(self) -> None:
        object.__setattr__(self, "k_kj_mol", _finite_real(self.k_kj_mol, "k_kj_mol"))
        if isinstance(self.periodicity, bool) or not isinstance(self.periodicity, int):
            raise ValueError("periodicity must be a positive integer")
        if self.periodicity <= 0:
            raise ValueError("periodicity must be a positive integer")
        object.__setattr__(self, "phase_rad", _finite_real(self.phase_rad, "phase_rad"))

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> FourierTerm:
        if not isinstance(data, Mapping):
            raise ValueError("each Fourier term must be a mapping")
        unknown = set(data) - {"k_kj_mol", "periodicity", "phase_rad"}
        missing = {"k_kj_mol", "periodicity"} - set(data)
        if unknown or missing:
            raise ValueError(f"invalid Fourier term fields; missing={sorted(missing)}, unknown={sorted(unknown)}")
        return cls(data["k_kj_mol"], data["periodicity"], data.get("phase_rad", 0.0))

    def as_dict(self) -> dict[str, float | int]:
        return {"k_kj_mol": self.k_kj_mol, "periodicity": self.periodicity, "phase_rad": self.phase_rad}


def _minimum_indices(values: list[float], tolerance: float) -> list[int]:
    """Return one representative per sampled cyclic minimum plateau.

    Adjacent samples within the reported tolerance form a plateau. A plateau
    is a minimum only when both exterior neighbours are above its minimum.
    A flat cyclic profile has no isolated minima and returns an empty list.
    The representatives are grid indices, not optimized continuous extrema.
    """
    size = len(values)
    if max(values) - min(values) <= tolerance:
        return []
    starts = [i for i in range(size) if abs(values[i] - values[(i - 1) % size]) > tolerance]
    if not starts:
        return []
    groups: list[list[int]] = []
    current: list[int] = []
    for offset in range(size):
        index = (starts[0] + offset) % size
        if current and abs(values[index] - values[current[-1]]) > tolerance:
            groups.append(current)
            current = []
        current.append(index)
    groups.append(current)
    minima: list[int] = []
    for group in groups:
        representative = min(group, key=lambda i: (values[i], i))
        minimum = values[representative]
        left = values[(group[0] - 1) % size]
        right = values[(group[-1] + 1) % size]
        if left > minimum + tolerance and right > minimum + tolerance:
            minima.append(representative)
    return sorted(minima)


@dataclass(frozen=True)
class TorsionPotential:
    """An immutable nonempty collection of complete Fourier terms."""

    terms: tuple[FourierTerm, ...]

    def __post_init__(self) -> None:
        try:
            terms = tuple(self.terms)
        except TypeError as exc:
            raise ValueError("terms must be a nonempty collection of FourierTerm objects") from exc
        if not terms or any(not isinstance(term, FourierTerm) for term in terms):
            raise ValueError("terms must be a nonempty collection of FourierTerm objects")
        object.__setattr__(self, "terms", terms)

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> TorsionPotential:
        """Read ``terms`` and an optional explicit ``convention``.

        Unknown fields are rejected to catch unit and phase-name mistakes.
        ``phase_rad`` can be omitted within a term and then means zero radians.
        """
        if not isinstance(data, Mapping):
            raise ValueError("torsion potential must be a mapping")
        unknown = set(data) - {"terms", "convention"}
        if unknown or "terms" not in data:
            raise ValueError(f"potential requires terms and permits only terms/convention; unknown={sorted(unknown)}")
        if data.get("convention", CONVENTION) != CONVENTION:
            raise ValueError(f"convention must be {CONVENTION!r}; convert source force-field terms explicitly")
        terms = data["terms"]
        if isinstance(terms, (str, bytes)) or not isinstance(terms, Sequence):
            raise ValueError("terms must be a nonempty sequence of mappings")
        return cls(tuple(FourierTerm.from_mapping(term) for term in terms))

    def energy(self, phi_rad: float) -> float:
        """Evaluate the periodic potential in kJ/mol at a finite angle."""
        phi = math.remainder(_finite_real(phi_rad, "phi_rad"), math.tau)
        try:
            value = math.fsum(
                term.k_kj_mol
                * (1.0 + math.cos(term.periodicity * phi - math.remainder(term.phase_rad, math.tau)))
                for term in self.terms
            )
        except (OverflowError, ValueError) as exc:
            raise ValueError("torsion energy exceeds the usable floating-point range") from exc
        if not math.isfinite(value):
            raise ValueError("torsion energy exceeds the usable floating-point range")
        return value

    def describe(self, n_grid: int = 720, temperature_k: float = 300.0) -> dict[str, Any]:
        """Sample a complete cyclic grid and compute transparent descriptors.

        The grid has no duplicated 2*pi endpoint and requires at least 36 bins
        and four bins per largest Fourier periodicity. Extrema and barriers
        are sampled-grid estimates; they are not exact optimized extrema.
        Boltzmann values are normalized probabilities of equally sized bins.
        """
        if isinstance(n_grid, bool) or not isinstance(n_grid, int) or n_grid < 36:
            raise ValueError("n_grid must be an integer of at least 36")
        if n_grid < 4 * max(term.periodicity for term in self.terms):
            raise ValueError("n_grid must provide at least four bins per largest Fourier periodicity")
        temperature = _finite_real(temperature_k, "temperature_k")
        thermal_energy = GAS_CONSTANT_KJ_MOL_K * temperature
        if temperature <= 0 or not math.isfinite(thermal_energy) or thermal_energy <= 0:
            raise ValueError("temperature_k must be positive and numerically resolvable")
        angles = [math.tau * i / n_grid for i in range(n_grid)]
        energies = [self.energy(angle) for angle in angles]
        minimum, maximum = min(energies), max(energies)
        barrier = maximum - minimum
        if not math.isfinite(barrier):
            raise ValueError("relative torsion energies exceed the usable floating-point range")
        relative = [energy - minimum for energy in energies]
        weights = [math.exp(-energy / thermal_energy) for energy in relative]
        total = math.fsum(weights)
        probabilities = [weight / total for weight in weights]
        entropy = -math.fsum(p * math.log(p) for p in probabilities if p > 0)
        scale = max(abs(minimum), abs(maximum))
        tolerance = max(64.0 * math.ulp(scale), 1e-12 * barrier)
        minima = _minimum_indices(energies, tolerance)
        return {
            "schema_version": "0.1",
            "result_kind": "analytic_torsion_grid",
            "convention": CONVENTION,
            "terms": [term.as_dict() for term in self.terms],
            "n_grid": n_grid,
            "temperature_k": temperature,
            "gas_constant_kj_mol_k": GAS_CONSTANT_KJ_MOL_K,
            "energy_min_grid_kj_mol": minimum,
            "energy_max_grid_kj_mol": maximum,
            "barrier_grid_kj_mol": barrier,
            "n_minima_grid": len(minima),
            "minima_angles_grid_rad": [angles[index] for index in minima],
            "minimum_plateau_tolerance_kj_mol": tolerance,
            "mean_relative_energy_grid_kj_mol": math.fsum(p * energy for p, energy in zip(probabilities, relative)),
            "shannon_grid_entropy_nats": entropy,
            "shannon_grid_entropy_fraction_of_uniform": entropy / math.log(n_grid),
            "angles_rad": angles,
            "relative_energies_kj_mol": relative,
            "boltzmann_bin_probabilities": probabilities,
            "interpretation": {
                "extrema": "sampled grid estimates, not optimized continuous extrema",
                "probabilities": "isolated one-dimensional model on equally sized angular bins; not condensed-phase polymer populations",
                "entropy": "discrete grid Shannon entropy, dependent on n_grid; not thermodynamic entropy",
                "provenance": "analytic model evaluation; no QM, MD, or experimental measurement was performed",
            },
        }
