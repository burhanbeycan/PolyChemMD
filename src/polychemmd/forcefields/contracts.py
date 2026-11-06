"""Explicit, independently implemented force-field parameter contracts.

These contracts validate the supplied definitions and assigned-type coverage.
They do not infer chemical atom types, generate omitted topology interactions,
assign partial charges, or establish a force field's physical accuracy.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
import json
import math
from numbers import Real
from typing import Any

from ..descriptors.torsion import TorsionPotential

KINDS = ("atoms", "bonds", "angles", "dihedrals", "impropers")
UNITS = {"energy": "kJ/mol", "length": "angstrom", "angle": "radian",
         "mass": "g/mol", "charge": "e"}
FORMS = {
    "atoms": "lj12_6",
    "bonds": "harmonic_k_delta_r_squared",
    "angles": "harmonic_k_delta_theta_squared",
    "dihedrals": "fourier_sum_k_1_plus_cos",
    "impropers": "harmonic_k_delta_chi_squared",
}


def _fields(data: Any, expected: set[str], name: str) -> Mapping[str, Any]:
    if not isinstance(data, Mapping) or any(not isinstance(k, str) for k in data):
        raise ValueError(f"{name} must be an object with string keys")
    missing, unknown = expected - set(data), set(data) - expected
    if missing or unknown:
        raise ValueError(f"{name}: missing={sorted(missing)}, unknown={sorted(unknown)}")
    return data


def _text(value: Any, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a nonempty string")
    return value


def _number(value: Any, name: str, lower: float | None = None,
            upper: float | None = None, strictly_positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite real number") from exc
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    if strictly_positive and result <= 0:
        raise ValueError(f"{name} must be positive")
    if lower is not None and result < lower or upper is not None and result > upper:
        raise ValueError(f"{name} is outside the supported range")
    return result


@dataclass(frozen=True)
class ParameterSource:
    kind: str
    reference: str

    def __post_init__(self):
        if self.kind not in ("analytic_demo", "literature", "qm_fit", "user_supplied"):
            raise ValueError("unsupported parameter source kind")
        _text(self.reference, "source reference")

    @classmethod
    def from_mapping(cls, data):
        _fields(data, {"kind", "reference"}, "parameter source")
        if data["kind"] not in ("analytic_demo", "literature", "qm_fit", "user_supplied"):
            raise ValueError("unsupported parameter source kind")
        return cls(data["kind"], _text(data["reference"], "source reference"))

    def to_dict(self):
        return {"kind": self.kind, "reference": self.reference}


@dataclass(frozen=True)
class NonbondedConvention:
    mixing_rule: str
    lj_14_scale: float
    coulomb_14_scale: float

    def __post_init__(self):
        if self.mixing_rule not in ("lorentz_berthelot", "geometric"):
            raise ValueError("unsupported mixing rule")
        object.__setattr__(self, "lj_14_scale", _number(self.lj_14_scale, "1-4 LJ scale", 0, 1))
        object.__setattr__(self, "coulomb_14_scale", _number(self.coulomb_14_scale, "1-4 Coulomb scale", 0, 1))

    @classmethod
    def from_mapping(cls, data):
        _fields(data, {"mixing_rule", "excluded_bonded_distances", "one_four"}, "nonbonded convention")
        if data["mixing_rule"] not in ("lorentz_berthelot", "geometric"):
            raise ValueError("unsupported mixing rule; explicit conversion is required")
        distances = data["excluded_bonded_distances"]
        if not isinstance(distances, list) or any(type(v) is not int for v in distances) or sorted(distances) != [1, 2]:
            raise ValueError("exclude exactly bonded distances [1, 2]; 1-4 interactions must use one_four once")
        scaling = _fields(data["one_four"], {"lj_scale", "coulomb_scale", "handling"}, "one_four")
        if scaling["handling"] != "scaled_pair_once":
            raise ValueError("1-4 handling must be scaled_pair_once; overlapping explicit and scaled terms are unsupported")
        return cls(data["mixing_rule"], _number(scaling["lj_scale"], "1-4 LJ scale", 0, 1),
                   _number(scaling["coulomb_scale"], "1-4 Coulomb scale", 0, 1))

    def to_dict(self):
        return {"mixing_rule": self.mixing_rule, "excluded_bonded_distances": [1, 2],
                "one_four": {"lj_scale": self.lj_14_scale,
                             "coulomb_scale": self.coulomb_14_scale,
                             "handling": "scaled_pair_once"}}


def _values(kind, values):
    if kind == "atoms":
        _fields(values, {"mass_g_mol", "epsilon_kj_mol", "sigma_angstrom"}, "LJ atom values")
        return {"mass_g_mol": _number(values["mass_g_mol"], "mass", strictly_positive=True),
                "epsilon_kj_mol": _number(values["epsilon_kj_mol"], "epsilon", lower=0),
                "sigma_angstrom": _number(values["sigma_angstrom"], "sigma", strictly_positive=True)}
    if kind == "bonds":
        _fields(values, {"k_kj_mol_angstrom2", "r0_angstrom"}, "harmonic bond values")
        return {"k_kj_mol_angstrom2": _number(values["k_kj_mol_angstrom2"], "bond k", strictly_positive=True),
                "r0_angstrom": _number(values["r0_angstrom"], "r0", strictly_positive=True)}
    if kind == "angles":
        _fields(values, {"k_kj_mol_rad2", "theta0_rad"}, "harmonic angle values")
        theta = _number(values["theta0_rad"], "theta0", lower=0, upper=math.pi)
        if theta == 0:
            raise ValueError("theta0 must be greater than zero")
        return {"k_kj_mol_rad2": _number(values["k_kj_mol_rad2"], "angle k", strictly_positive=True),
                "theta0_rad": theta}
    if kind == "dihedrals":
        _fields(values, {"terms"}, "Fourier values")
        potential = TorsionPotential.from_mapping(values)
        return {"terms": [term.as_dict() for term in potential.terms]}
    _fields(values, {"k_kj_mol_rad2", "chi0_rad"}, "harmonic improper values")
    return {"k_kj_mol_rad2": _number(values["k_kj_mol_rad2"], "improper k", strictly_positive=True),
            "chi0_rad": _number(values["chi0_rad"], "chi0", lower=-math.pi, upper=math.pi)}


@dataclass(frozen=True)
class ParameterDefinition:
    kind: str
    type_id: str
    form: str
    values_json: str
    source: ParameterSource

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError("unsupported parameter interaction kind")
        _text(self.type_id, "type_id")
        if self.form != FORMS[self.kind]:
            raise ValueError("unsupported explicit parameter form")
        if not isinstance(self.source, ParameterSource):
            raise ValueError("source must be a validated ParameterSource")
        values = _values(self.kind, json.loads(_text(self.values_json, "values_json")))
        object.__setattr__(self, "values_json", json.dumps(values, sort_keys=True, allow_nan=False))

    @classmethod
    def from_mapping(cls, kind, data):
        _fields(data, {"type_id", "form", "values", "source"}, f"{kind} definition")
        if data["form"] != FORMS[kind]:
            raise ValueError(f"{kind} requires explicit supported form {FORMS[kind]!r}")
        values = _values(kind, data["values"])
        return cls(kind, _text(data["type_id"], "type_id"), data["form"],
                   json.dumps(values, sort_keys=True, allow_nan=False),
                   ParameterSource.from_mapping(data["source"]))

    @property
    def values(self):
        """Return a fresh copy so public callers cannot mutate a definition."""
        return json.loads(self.values_json)

    def to_dict(self):
        return {"type_id": self.type_id, "form": self.form, "values": self.values,
                "source": self.source.to_dict()}


@dataclass(frozen=True)
class ParameterSet:
    name: str
    source: ParameterSource
    nonbonded: NonbondedConvention
    definitions: tuple[ParameterDefinition, ...]

    def __post_init__(self):
        _text(self.name, "parameter-set name")
        if not isinstance(self.source, ParameterSource) or not isinstance(self.nonbonded, NonbondedConvention):
            raise ValueError("validated source and nonbonded conventions are required")
        try:
            definitions = tuple(self.definitions)
        except TypeError as exc:
            raise ValueError("definitions must contain ParameterDefinition objects") from exc
        if any(not isinstance(entry, ParameterDefinition) for entry in definitions):
            raise ValueError("definitions must contain ParameterDefinition objects")
        keys = [(entry.kind, entry.type_id) for entry in definitions]
        if len(set(keys)) != len(keys):
            raise ValueError("duplicate parameter type_id within an interaction kind")
        object.__setattr__(self, "definitions", definitions)

    @classmethod
    def from_mapping(cls, data):
        _fields(data, {"schema_version", "name", "units", "source", "nonbonded", "parameters"}, "parameter set")
        if data["schema_version"] != "0.1":
            raise ValueError("unsupported parameter-set schema_version")
        if _fields(data["units"], set(UNITS), "units") != UNITS:
            raise ValueError(f"mixed/unsupported units; explicitly convert to {UNITS}")
        registries = _fields(data["parameters"], set(KINDS), "parameter registries")
        definitions = []
        for kind in KINDS:
            if not isinstance(registries[kind], list):
                raise ValueError(f"{kind} parameters must be an array")
            entries = [ParameterDefinition.from_mapping(kind, row) for row in registries[kind]]
            if len({row.type_id for row in entries}) != len(entries):
                raise ValueError(f"duplicate {kind} parameter type_id")
            definitions.extend(entries)
        return cls(_text(data["name"], "parameter-set name"), ParameterSource.from_mapping(data["source"]),
                   NonbondedConvention.from_mapping(data["nonbonded"]), tuple(definitions))

    def to_dict(self):
        return {"schema_version": "0.1", "name": self.name, "units": dict(UNITS),
                "source": self.source.to_dict(), "nonbonded": self.nonbonded.to_dict(),
                "parameters": {kind: [entry.to_dict() for entry in self.definitions if entry.kind == kind] for kind in KINDS}}


@dataclass(frozen=True)
class TopologyRequirement:
    kind: str
    interaction_id: str
    type_id: str

    def __post_init__(self):
        if self.kind not in KINDS:
            raise ValueError("unsupported topology interaction kind")
        _text(self.interaction_id, "interaction id")
        _text(self.type_id, "assigned type_id")


@dataclass(frozen=True)
class TopologyRequirements:
    entries: tuple[TopologyRequirement, ...]

    def __post_init__(self):
        try:
            entries = tuple(self.entries)
        except TypeError as exc:
            raise ValueError("entries must contain TopologyRequirement objects") from exc
        if any(not isinstance(row, TopologyRequirement) for row in entries):
            raise ValueError("entries must contain TopologyRequirement objects")
        if not any(row.kind == "atoms" for row in entries):
            raise ValueError("assigned topology must contain atoms")
        identities = [(row.kind, row.interaction_id) for row in entries]
        if len(set(identities)) != len(identities):
            raise ValueError("duplicate topology interaction id")
        object.__setattr__(self, "entries", entries)

    @classmethod
    def from_mapping(cls, data):
        _fields(data, set(KINDS), "assigned topology")
        entries = []
        for kind in KINDS:
            if not isinstance(data[kind], list):
                raise ValueError(f"topology {kind} must be an array")
            seen = set()
            for row in data[kind]:
                _fields(row, {"id", "type_id"}, f"topology {kind} entry")
                identity = row["id"]
                if type(identity) is int and identity > 0:
                    identity = str(identity)
                else:
                    identity = _text(identity, "interaction id")
                if identity in seen:
                    raise ValueError(f"duplicate topology {kind} id {identity!r}")
                seen.add(identity)
                entries.append(TopologyRequirement(kind, identity, _text(row["type_id"], "assigned type_id")))
        if not data["atoms"]:
            raise ValueError("assigned topology must contain atoms")
        return cls(tuple(entries))

    def to_dict(self):
        return {kind: [{"id": row.interaction_id, "type_id": row.type_id}
                       for row in self.entries if row.kind == kind] for kind in KINDS}


class MissingParametersError(ValueError):
    """The supplied assigned topology has parameter types with no definition."""

    def __init__(self, report):
        self.report = report
        details = "; ".join(f"{kind}: {', '.join(types)}" for kind, types in report.missing_types if types)
        super().__init__(f"parameter coverage is incomplete ({details})")


@dataclass(frozen=True)
class CoverageReport:
    parameter_set_name: str
    topology_counts: tuple[tuple[str, int], ...]
    missing_types: tuple[tuple[str, tuple[str, ...]], ...]
    missing_interactions: tuple[TopologyRequirement, ...]
    unused_types: tuple[tuple[str, tuple[str, ...]], ...]

    @property
    def complete(self):
        return not self.missing_interactions

    def require_complete(self):
        """Gate assigned-type coverage only; this does not certify production MD."""
        if not self.complete:
            raise MissingParametersError(self)
        return self

    def to_dict(self):
        return {"schema_version": "0.1", "parameter_set_name": self.parameter_set_name,
                "assigned_type_coverage_complete": self.complete,
                "topology_counts": dict(self.topology_counts),
                "missing_types": {kind: list(types) for kind, types in self.missing_types},
                "missing_interactions": [{"kind": row.kind, "id": row.interaction_id, "type_id": row.type_id} for row in self.missing_interactions],
                "unused_types": {kind: list(types) for kind, types in self.unused_types},
                "chemical_typing": "not_performed", "physical_accuracy": "not_assessed",
                "production_ready": False,
                "scope": "Coverage of the supplied assigned interaction records only; graph completeness, atom typing, partial charges, coordinates and physical validity are not established."}


def check_coverage(parameter_set, topology):
    """Return missing and unused types without silently inventing assignments."""
    if not isinstance(parameter_set, ParameterSet):
        parameter_set = ParameterSet.from_mapping(parameter_set)
    if not isinstance(topology, TopologyRequirements):
        topology = TopologyRequirements.from_mapping(topology)
    defined = {kind: {row.type_id for row in parameter_set.definitions if row.kind == kind} for kind in KINDS}
    used = {kind: {row.type_id for row in topology.entries if row.kind == kind} for kind in KINDS}
    missing = tuple(row for row in topology.entries if row.type_id not in defined[row.kind])
    return CoverageReport(parameter_set.name,
                          tuple((kind, sum(row.kind == kind for row in topology.entries)) for kind in KINDS),
                          tuple((kind, tuple(sorted(used[kind] - defined[kind]))) for kind in KINDS), missing,
                          tuple((kind, tuple(sorted(defined[kind] - used[kind]))) for kind in KINDS))
