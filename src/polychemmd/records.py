"""Strict metadata validation; molecular graph/valence validation is a later task.

SMILES are retained as supplied and MUST NOT be considered chemically validated
by these classes. This distinction is exposed in the public output.
"""

from dataclasses import dataclass, asdict
import math
from collections.abc import Mapping


def _mapping(value, name):
    if not isinstance(value, Mapping):
        raise ValueError(f"{name} must be an object")
    return value


def _keys(data, allowed, required, name):
    _mapping(data, name)
    unknown, missing = set(data) - set(allowed), set(required) - set(data)
    if unknown or missing:
        raise ValueError(f"{name}: unknown keys={sorted(unknown)}; missing keys={sorted(missing)}")


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _real(value, name, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        converted = float(value)
    except (OverflowError, ValueError) as exc:
        raise ValueError(f"{name} must be representable as a finite number") from exc
    if not math.isfinite(converted):
        raise ValueError(f"{name} must be a finite number")
    if positive and converted <= 0:
        raise ValueError(f"{name} must be positive")
    return converted


@dataclass(frozen=True)
class RepeatUnit:
    label: str
    smiles: str
    mole_fraction: float

    @classmethod
    def from_mapping(cls, data):
        _keys(data, {"label", "smiles", "mole_fraction"}, {"label", "smiles", "mole_fraction"}, "repeat unit")
        fraction = _real(data["mole_fraction"], "mole_fraction", positive=True)
        if fraction > 1:
            raise ValueError("mole_fraction must not exceed 1")
        return cls(_text(data["label"], "label"), _text(data["smiles"], "smiles"), fraction)


@dataclass(frozen=True)
class Microstructure:
    architecture: str
    sequence_type: str
    tacticity: str = "unspecified"
    degree_of_polymerization: int | None = None

    @classmethod
    def from_mapping(cls, data):
        _keys(data, {"architecture", "sequence_type", "tacticity", "degree_of_polymerization"}, {"architecture", "sequence_type"}, "microstructure")
        architecture = data["architecture"]
        sequence = data["sequence_type"]
        tacticity = data.get("tacticity", "unspecified")
        if architecture not in ("linear", "star", "graft", "comb", "network", "cyclic"):
            raise ValueError("unsupported architecture label")
        if sequence not in ("homopolymer", "alternating", "random", "block", "gradient", "explicit"):
            raise ValueError("unsupported sequence_type")
        if tacticity not in ("unspecified", "atactic", "isotactic", "syndiotactic", "explicit"):
            raise ValueError("unsupported tacticity")
        dp = data.get("degree_of_polymerization")
        if dp is not None and (type(dp) is not int or dp < 1):
            raise ValueError("degree_of_polymerization must be a positive integer or null")
        return cls(architecture, sequence, tacticity, dp)


@dataclass(frozen=True)
class StatePoint:
    phase: str
    temperature_k: float
    pressure_bar: float
    preparation_history: str

    @classmethod
    def from_mapping(cls, data):
        keys = {"phase", "temperature_k", "pressure_bar", "preparation_history"}
        _keys(data, keys, keys, "state")
        if data["phase"] not in ("amorphous", "melt", "solution", "mixture", "crystalline", "oriented", "network", "interface"):
            raise ValueError("unsupported phase label")
        return cls(data["phase"], _real(data["temperature_k"], "temperature_k", positive=True),
                   _real(data["pressure_bar"], "pressure_bar"),
                   _text(data["preparation_history"], "preparation_history"))


@dataclass(frozen=True)
class MaterialSpec:
    schema_version: str
    material_id: str
    chemical_family: str
    repeat_units: tuple[RepeatUnit, ...]
    microstructure: Microstructure
    state: StatePoint
    source_reference: str

    @classmethod
    def from_mapping(cls, data):
        keys = {"schema_version", "material_id", "chemical_family", "repeat_units", "microstructure", "state", "source_reference"}
        _keys(data, keys, keys, "material")
        if data["schema_version"] != "0.1":
            raise ValueError("unsupported schema_version; expected 0.1")
        if not isinstance(data["repeat_units"], list) or not data["repeat_units"]:
            raise ValueError("repeat_units must be a non-empty array")
        units = tuple(RepeatUnit.from_mapping(x) for x in data["repeat_units"])
        if len({u.label for u in units}) != len(units):
            raise ValueError("repeat-unit labels must be unique")
        if not math.isclose(sum(u.mole_fraction for u in units), 1.0, rel_tol=0, abs_tol=1e-8):
            raise ValueError("repeat-unit mole fractions must sum to 1")
        micro = Microstructure.from_mapping(data["microstructure"])
        if micro.sequence_type == "homopolymer" and len(units) != 1:
            raise ValueError("homopolymer records require exactly one repeat unit")
        if micro.sequence_type == "alternating" and (len(units) != 2 or any(not math.isclose(u.mole_fraction, .5, abs_tol=1e-8) for u in units)):
            raise ValueError("this schema defines alternating as binary 1:1; use explicit for other patterns")
        return cls("0.1", _text(data["material_id"], "material_id"),
                   _text(data["chemical_family"], "chemical_family"), units, micro,
                   StatePoint.from_mapping(data["state"]),
                   _text(data["source_reference"], "source_reference"))

    def to_dict(self):
        data = asdict(self)
        data["repeat_units"] = list(data["repeat_units"])
        return data
