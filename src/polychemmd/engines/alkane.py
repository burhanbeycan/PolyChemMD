"""Explicit isolated n-alkane LAMMPS single points and energy/force parity.

This adapter exports the supported assigned OPLS-AA numerical subset. It runs
zero integration steps, not a relaxation, equilibrium or density protocol.
"""

from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import re

from ..forcefields.contracts import ParameterSet
from ..geometry import MolecularGeometry
from .lammps import LammpsError, LammpsRunSpec, run_lammps


KCAL_TO_KJ = 4.184
_MODEL = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"
_KINDS = (("bonds", 2), ("angles", 3), ("dihedrals", 4))
_THERMO_COLUMNS = ("Step", "Temp", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "PotEng")
_COMPONENTS = {"bond": "E_bond", "angle": "E_angle", "dihedral": "E_dihed",
               "lj": "E_vdwl", "coulomb": "E_coul", "total": "PotEng"}
_PINNED_OPLS_SINE_FLOOR = 1e-5
_PINNED_OPLS_SOURCE = (
    "https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/"
    "src/MOLECULE/dihedral_opls.cpp#L171-L180"
)


class AlkaneParityError(LammpsError):
    """Actual engine output failed the declared numerical parity checks."""


def _finite(value, name, *, lower=None, positive=False):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a finite number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(value) or (positive and value <= 0) or (lower is not None and value < lower):
        raise ValueError(f"{name} is outside the supported finite range")
    return value


def _format(value):
    return format(_finite(value, "exported coefficient or coordinate"), ".17g")


def _assignment(assigned):
    from ..forcefields.alkanes import AssignedAlkane
    if not isinstance(assigned, AssignedAlkane):
        raise ValueError("assigned must be an AssignedAlkane result")
    record = assigned.to_dict()
    if record.get("kind") != "assigned_opls_aa_n_alkane" or record.get("model_id") != _MODEL:
        raise ValueError("unsupported assigned force-field model")
    geometry = MolecularGeometry.from_mapping(record["geometry"])
    parameters = ParameterSet.from_mapping(record["parameter_set"])
    if parameters.nonbonded.to_dict() != {
        "mixing_rule": "geometric", "excluded_bonded_distances": [1, 2],
        "one_four": {"lj_scale": .5, "coulomb_scale": .5, "handling": "scaled_pair_once"},
    }:
        raise ValueError("the export requires geometric OPLS mixing and one-half 1-4 scaling")
    definitions = {
        kind: {p.type_id: p for p in parameters.definitions if p.kind == kind}
        for kind in ("atoms", "bonds", "angles", "dihedrals", "impropers")
    }
    if record["interactions"].get("impropers"):
        raise ValueError("impropers are outside the supported saturated n-alkane model")
    ids = [atom.atom_id for atom in geometry.atoms]
    assignments = record["atom_assignments"]
    if not isinstance(assignments, list) or len(assignments) != len(ids):
        raise ValueError("assigned atoms must match the geometry")
    by_id = {a["atom_id"]: a for a in assignments}
    if len(by_id) != len(ids) or set(by_id) != set(ids):
        raise ValueError("atom assignments lost or duplicate geometry identities")
    for atom in geometry.atoms:
        row = by_id[atom.atom_id]
        if row["atom_map"] != atom.atom_map or row["type_id"] not in definitions["atoms"]:
            raise ValueError("invalid atom assignment or parameter reference")
        _finite(row["charge_e"], "partial charge")
    for kind, arity in _KINDS:
        rows = record["interactions"][kind]
        if not isinstance(rows, list):
            raise ValueError(f"assigned {kind} must be an array")
        identities = set()
        for row in rows:
            endpoints = row["atom_ids"]
            if (not isinstance(endpoints, list) or len(endpoints) != arity
                    or len(set(endpoints)) != arity or not set(endpoints).issubset(by_id)
                    or row["type_id"] not in definitions[kind] or row["id"] in identities):
                raise ValueError(f"invalid {kind} interaction or parameter reference")
            identities.add(row["id"])
    return record, geometry, definitions, by_id


def _opls_coefficients(values):
    """Convert canonical Fourier terms to explicit LAMMPS OPLS coefficients."""
    coefficients = [0.0] * 4
    for term in values["terms"]:
        n = term["periodicity"]
        if type(n) is not int or not 1 <= n <= 4:
            raise ValueError("LAMMPS OPLS export supports only periodicities 1 through 4")
        k = _finite(term["k_kj_mol"], "Fourier coefficient")
        phase = _finite(term["phase_rad"], "Fourier phase")
        expected = 0.0 if n % 2 else math.pi
        if k and abs(math.remainder(phase - expected, math.tau)) > 1e-12:
            raise ValueError("Fourier phase cannot be represented by the OPLS convention")
        coefficients[n-1] += 2.0 * k / KCAL_TO_KJ
    if not all(math.isfinite(v) for v in coefficients):
        raise ValueError("OPLS coefficients exceed the finite numeric range")
    return tuple(coefficients)


def _geometric_abs_sin_dihedral(points):
    """Geometric |sin(phi)| from unit planes, without an acos subtraction.

    This screen neither reproduces nor modifies the engine derivative.
    Normalizing bonds first avoids unnecessarily large intermediates.
    """
    def unit(vector):
        length = math.hypot(*vector)
        if not math.isfinite(length) or length <= 0:
            raise ValueError("undefined geometric torsion in backend diagnostic")
        return tuple(v/length for v in vector)
    def cross(a, b):
        return (a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0])
    bonds = [unit(tuple(b-a for a, b in zip(points[i], points[i+1]))) for i in range(3)]
    first = unit(cross(bonds[0], bonds[1]))
    second = unit(cross(bonds[1], bonds[2]))
    return min(1.0, math.hypot(*cross(first, second)))


def _opls_backend_diagnostics(assignment, geometry, engine_version_banner=None):
    """Record planar geometries relevant to one pinned OPLS implementation.

    Its threshold is a cited backend detail, not a universal engine rule or a
    change to the reference potential or acceptance tolerances.
    """
    positions = {atom.atom_id: point for atom, point in zip(geometry.atoms, geometry.coordinates_angstrom)}
    near_planar = []
    for row in assignment["interactions"]["dihedrals"]:
        value = _geometric_abs_sin_dihedral([positions[identity] for identity in row["atom_ids"]])
        if value < _PINNED_OPLS_SINE_FLOOR:
            near_planar.append({"interaction_id": row["id"], "atom_ids": list(row["atom_ids"]),
                                "type_id": row["type_id"], "abs_sin_phi": value})
    return {
        "kind": "near_planar_opls_backend_numerical_screen",
        "screen_method": "absolute cross product of unit adjacent-plane normals",
        "pinned_source_engine": "LAMMPS 22 Jul 2025 Update 4, dihedral_style opls",
        "pinned_source_url": _PINNED_OPLS_SOURCE,
        "observed_source_sine_threshold": _PINNED_OPLS_SINE_FLOOR,
        "threshold_scope": "observed in the cited pinned source; not guaranteed for other engine versions or builds",
        "actual_engine_version_banner": engine_version_banner,
        "near_planar_dihedrals": near_planar,
        "interpretation": (
            "The pinned OPLS implementation floors the sine used in derivative ratios near planar torsions; "
            "force differences can exceed roundoff while energies agree. This screen does not establish "
            "that another build uses the same numerical treatment."
        ),
        "reference_or_potential_modified": False,
        "parity_tolerances_modified": False,
    }


def build_alkane_input(assigned, *, cutoff_angstrom=12.0):
    """Return ``(script, data)`` without writing files or running an engine.

    Coordinates remain in their input frame. A nonperiodic orthorhombic box
    surrounds them by more than the requested unshifted LJ/Coulomb cutoff.
    """
    cutoff = _finite(cutoff_angstrom, "cutoff_angstrom", positive=True)
    record, geometry, definitions, atom_assignments = _assignment(assigned)
    type_maps = {kind: {name: i+1 for i, name in enumerate(sorted(rows))}
                 for kind, rows in definitions.items()}
    atom_ids = {atom.atom_id: i+1 for i, atom in enumerate(geometry.atoms)}
    points = geometry.coordinates_angstrom
    margin = cutoff + 5.0
    bounds = [(min(p[axis] for p in points) - margin,
               max(p[axis] for p in points) + margin) for axis in range(3)]
    if any(not math.isfinite(v) for pair in bounds for v in pair):
        raise ValueError("single-point box bounds exceed the finite numeric range")
    lines = ["PolyChemMD isolated OPLS-AA n-alkane single point", "",
             f"{len(points)} atoms"]
    for kind, _ in _KINDS:
        lines.append(f"{len(record['interactions'][kind])} {kind}")
    lines.extend(["0 impropers", "", f"{len(type_maps['atoms'])} atom types"])
    for kind, _ in _KINDS:
        singular = {"bonds": "bond", "angles": "angle", "dihedrals": "dihedral"}[kind]
        lines.append(f"{len(type_maps[kind])} {singular} types")
    lines.extend(["0 improper types", ""])
    for axis, (lower, upper) in zip("xyz", bounds):
        lines.append(f"{_format(lower)} {_format(upper)} {axis}lo {axis}hi")

    lines.extend(["", "Masses", ""])
    for name, type_number in type_maps["atoms"].items():
        values = definitions["atoms"][name].values
        lines.append(f"{type_number} {_format(values['mass_g_mol'])} # {name}")
    lines.extend(["", "Pair Coeffs # lj/cut/coul/cut", ""])
    for name, type_number in type_maps["atoms"].items():
        values = definitions["atoms"][name].values
        lines.append(f"{type_number} {_format(values['epsilon_kj_mol']/KCAL_TO_KJ)} "
                     f"{_format(values['sigma_angstrom'])} # {name}")
    lines.extend(["", "Bond Coeffs # harmonic", ""])
    for name, type_number in type_maps["bonds"].items():
        values = definitions["bonds"][name].values
        lines.append(f"{type_number} {_format(values['k_kj_mol_angstrom2']/KCAL_TO_KJ)} "
                     f"{_format(values['r0_angstrom'])} # {name}")
    lines.extend(["", "Angle Coeffs # harmonic", ""])
    for name, type_number in type_maps["angles"].items():
        values = definitions["angles"][name].values
        lines.append(f"{type_number} {_format(values['k_kj_mol_rad2']/KCAL_TO_KJ)} "
                     f"{_format(math.degrees(values['theta0_rad']))} # {name}")
    lines.extend(["", "Dihedral Coeffs # opls", ""])
    for name, type_number in type_maps["dihedrals"].items():
        coefficients = _opls_coefficients(definitions["dihedrals"][name].values)
        lines.append(f"{type_number} " + " ".join(_format(k) for k in coefficients) + f" # {name}")
    lines.extend(["", "Atoms # full", ""])
    for atom, point in zip(geometry.atoms, points):
        row = atom_assignments[atom.atom_id]
        lines.append(f"{atom_ids[atom.atom_id]} 1 {type_maps['atoms'][row['type_id']]} "
                     f"{_format(row['charge_e'])} " + " ".join(_format(v) for v in point))
    for kind, _ in _KINDS:
        lines.extend(["", kind.capitalize(), ""])
        for i, row in enumerate(record["interactions"][kind], 1):
            lines.append(f"{i} {type_maps[kind][row['type_id']]} "
                         + " ".join(str(atom_ids[identity]) for identity in row["atom_ids"]))
    data = "\n".join(lines) + "\n"
    script = "\n".join((
        "# Isolated assigned-force-field energy/force check; no MD integration.",
        "clear", "units real", "atom_style full", "boundary f f f",
        "bond_style harmonic", "angle_style harmonic", "dihedral_style opls",
        f"pair_style lj/cut/coul/cut {_format(cutoff)}",
        "pair_modify mix geometric shift no tail no", "dielectric 1.0",
        "special_bonds lj/coul 0.0 0.0 0.5", "read_data system.data",
        "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes",
        "thermo 1", "thermo_style custom step temp ebond eangle edihed evdwl ecoul pe",
        "thermo_modify norm no format float %.17g", "run 0",
        "write_dump all custom forces.dump id x y z fx fy fz modify sort id format float %.17g",
        'print "POLYCHEMMD_RUN_COMPLETE"', "",
    ))
    return script, data


def prepare_alkane_singlepoint(assigned, *, run_directory, executable=None,
                              timeout_seconds=30.0, cutoff_angstrom=12.0):
    """Create a bounded run specification; run nothing and create no files."""
    script, data = build_alkane_input(assigned, cutoff_angstrom=cutoff_angstrom)
    return LammpsRunSpec(input_script=script, run_directory=run_directory,
                         executable=executable, timeout_seconds=timeout_seconds,
                         expected_artifacts=("log.lammps", "forces.dump"),
                         required_thermo_columns=_THERMO_COLUMNS,
                         input_files=(("system.data", data),), units_style="real")


def parse_alkane_force_dump(text, expected_coordinates):
    """Parse exactly one step-zero custom dump and verify atom/coordinate identity.

    Forces are raw kcal/(mol angstrom). Rows may arrive out of order, but IDs
    must be exactly 1..N with no duplicates or lost atoms.
    """
    if not isinstance(text, str):
        raise ValueError("force dump must be text")
    expected = []
    for point in expected_coordinates:
        if not isinstance(point, (tuple, list)) or len(point) != 3:
            raise ValueError("expected coordinates must be triples")
        expected.append(tuple(_finite(v, "expected coordinate") for v in point))
    if not expected:
        raise ValueError("expected coordinates must not be empty")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    count = len(expected)
    if len(lines) != count + 9:
        raise ValueError("force dump must contain exactly one complete frame")
    if lines[0] != "ITEM: TIMESTEP" or lines[1] != "0" or lines[2] != "ITEM: NUMBER OF ATOMS":
        raise ValueError("force dump must be a fresh step-zero frame")
    if lines[3] != str(count) or lines[4] != "ITEM: BOX BOUNDS ff ff ff":
        raise ValueError("force dump atom count or nonperiodic boundary differs from the request")
    for row in lines[5:8]:
        fields = row.split()
        if len(fields) != 2:
            raise ValueError("unsupported force-dump box representation")
        try:
            lower, upper = (float(v) for v in fields)
        except (ValueError, OverflowError) as exc:
            raise ValueError("invalid force-dump box bounds") from exc
        if not math.isfinite(lower) or not math.isfinite(upper) or not lower < upper:
            raise ValueError("force-dump box bounds must be ordered and finite")
    if lines[8] != "ITEM: ATOMS id x y z fx fy fz":
        raise ValueError("unsupported force-dump columns")
    rows = {}
    for row in lines[9:]:
        fields = row.split()
        if len(fields) != 7 or re.fullmatch(r"[1-9][0-9]*", fields[0]) is None:
            raise ValueError("malformed force-dump atom row")
        identity = int(fields[0])
        if identity in rows or not 1 <= identity <= count:
            raise ValueError("force-dump atom IDs are duplicated or outside 1..N")
        try:
            values = tuple(float(v) for v in fields[1:])
        except (ValueError, OverflowError) as exc:
            raise ValueError("nonnumeric force-dump atom value") from exc
        if not all(math.isfinite(v) for v in values):
            raise ValueError("nonfinite force-dump atom value")
        if any(not math.isclose(x, y, rel_tol=1e-12, abs_tol=1e-10)
               for x, y in zip(values[:3], expected[identity-1])):
            raise ValueError("force-dump coordinates differ from the assigned atom geometry")
        rows[identity] = values[3:]
    if set(rows) != set(range(1, count+1)):
        raise ValueError("force dump lost one or more requested atom IDs")
    return tuple(rows[i] for i in range(1, count+1))


def verify_alkane_singlepoint(assigned, *, run_directory, executable=None,
                             timeout_seconds=30.0, cutoff_angstrom=12.0,
                             energy_abs_tol_kj_mol=1e-6,
                             force_abs_tol_kj_mol_angstrom=1e-5):
    """Execute actual LAMMPS and compare all components/forces with the reference.

    Declared absolute tolerances are checked before any engine execution. A
    missing/failed engine, malformed output or numerical mismatch raises an
    explicit error; none is converted into a successful parity report.
    """
    from ..forcefields.reference import evaluate_alkane
    energy_tolerance = _finite(energy_abs_tol_kj_mol, "energy tolerance", lower=0)
    force_tolerance = _finite(force_abs_tol_kj_mol_angstrom, "force tolerance", lower=0)
    assignment, geometry, _, _ = _assignment(assigned)
    reference = evaluate_alkane(assigned, cutoff_angstrom=cutoff_angstrom)
    spec = prepare_alkane_singlepoint(assigned, run_directory=run_directory,
                                     executable=executable, timeout_seconds=timeout_seconds,
                                     cutoff_angstrom=cutoff_angstrom)
    run = run_lammps(spec).to_dict()
    report = {
        "schema_version": "1", "status": "failed", "evidence_level": "engine_energy_force_parity",
        "model_id": assignment["model_id"], "boundary": "isolated_nonperiodic",
        "cutoff_angstrom": float(cutoff_angstrom),
        "nonbonded_policy": "unshifted finite-cutoff LJ and Coulomb; geometric mixing; scaled 1-4 pairs",
        "tolerances": {"energy_abs_kj_mol": energy_tolerance,
                       "force_component_abs_kj_mol_angstrom": force_tolerance},
        "reference": reference, "engine_run": run,
        "backend_numerical_diagnostics": _opls_backend_diagnostics(
            assignment, geometry, run.get("engine_version_banner")),
        "atom_ids": [a.atom_id for a in geometry.atoms],
        "engine_atom_identity": [{"engine_atom_id": i+1, "source_atom_id": atom.atom_id,
                                  "source_atom_map": atom.atom_map}
                                 for i, atom in enumerate(geometry.atoms)],
        "source_graph_sha256": geometry.source_graph_sha256,
        "assignment_canonical_sha256": hashlib.sha256(json.dumps(
            assignment, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
        "source_snapshot": assignment["source_snapshot"],
        "input_files": {"input.lammps_sha256": hashlib.sha256(spec.input_script.encode()).hexdigest(),
                        "system.data_sha256": hashlib.sha256(dict(spec.input_files)["system.data"].encode()).hexdigest()},
        "md_integration_steps": 0, "trajectory_production_performed": False,
        "equilibrium_status": "not_assessed", "force_field_physical_accuracy": "not_assessed",
        "material_property_validation": "not_performed",
    }
    try:
        blocks = run["thermo_blocks"]
        if len(blocks) != 1 or len(blocks[0]["rows"]) != 1:
            raise ValueError("single point requires exactly one thermo block and row")
        thermo = dict(zip(blocks[0]["columns"], blocks[0]["rows"][0]))
        if thermo["Step"] != 0:
            raise ValueError("single-point thermo step must be zero")
        energies = {name: thermo[column] * KCAL_TO_KJ for name, column in _COMPONENTS.items()}
        raw_forces = parse_alkane_force_dump(
            (spec.run_directory/"forces.dump").read_text(encoding="utf-8"), geometry.coordinates_angstrom)
        forces = [[value*KCAL_TO_KJ for value in row] for row in raw_forces]
        ref_energies, ref_forces = reference["energy_components_kj_mol"], reference["forces_kj_mol_angstrom"]
        if len(ref_forces) != len(forces) or any(len(row) != 3 for row in ref_forces):
            raise ValueError("reference forces do not match the assigned atom sequence")
        errors = {name: abs(value - ref_energies[name]) for name, value in energies.items()}
        force_errors = [abs(x-y) for actual, expected in zip(forces, ref_forces)
                        for x, y in zip(actual, expected)]
        if not all(math.isfinite(v) for v in [*errors.values(), *force_errors]):
            raise ValueError("energy/force differences are not finite")
        maximum_force_error = max(force_errors)
        rms_force_error = (maximum_force_error * math.sqrt(math.fsum(
            (v/maximum_force_error)**2 for v in force_errors)/len(force_errors))
                           if maximum_force_error else 0.0)
        report.update(engine_energy_components_kj_mol=energies,
                      engine_forces_kj_mol_angstrom=forces,
                      energy_component_abs_errors_kj_mol=errors,
                      max_force_component_abs_error_kj_mol_angstrom=maximum_force_error,
                      rms_force_component_error_kj_mol_angstrom=rms_force_error,
                      atom_and_coordinate_identity="passed")
        if any(error > energy_tolerance for error in errors.values()) or max(force_errors) > force_tolerance:
            raise ValueError("actual LAMMPS energy/force parity exceeds the predeclared tolerance")
        report["status"] = "passed"
    except (OSError, ValueError, KeyError, TypeError) as exc:
        report["failure_reason"] = str(exc)
        (spec.run_directory/"parity_record.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        raise AlkaneParityError(str(exc), record=report) from exc
    (spec.run_directory/"parity_record.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    return report
