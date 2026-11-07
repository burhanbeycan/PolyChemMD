"""Periodic OPLS-AA/Ewald input export and actual single-point comparison.

This defines a finite numerical periodic model. Run-zero agreement is neither
a production trajectory nor evidence of material-property accuracy.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
from types import SimpleNamespace

from ..forcefields.contracts import ParameterSet
from .alkane import (
    AlkaneParityError, KCAL_TO_KJ, _finite, _format, _opls_coefficients,
    _opls_backend_diagnostics,
)
from .lammps import LammpsError, LammpsRunSpec, run_lammps


_MODEL = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"
_KINDS = ("bonds", "angles", "dihedrals")
_THERMO = ("Step", "Temp", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "E_long", "PotEng")
_COMPONENTS = {"bond": "E_bond", "angle": "E_angle", "dihedral": "E_dihed",
               "lj": "E_vdwl", "coulomb_realspace": "E_coul",
               "coulomb_kspace": "E_long", "total": "PotEng"}


def _record(assigned):
    from ..forcefields.periodic import AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell
    if not isinstance(assigned, (AssignedPeriodicAlkaneCell, CompactAssignedPeriodicAlkaneCell)):
        raise ValueError("assigned must be a validated periodic alkane assignment")
    record = assigned.to_dict()
    if record.get("kind") not in ("assigned_periodic_opls_aa_n_alkane_cell",
            "compact_assigned_periodic_opls_aa_n_alkane_cell") or record.get("model_id") != _MODEL:
        raise ValueError("unsupported periodic assigned model")
    if not 0 < len(record["atoms"]) <= 8192:
        raise ValueError("the bounded periodic engine export supports 1 to 8192 atoms")
    return record


def _settings(record, cutoff, alpha, nmax):
    cutoff = _finite(cutoff, "cutoff_angstrom", positive=True)
    alpha = _finite(alpha, "ewald_alpha_per_angstrom", positive=True)
    if type(nmax) is not int or not 1 <= nmax <= 20:
        raise ValueError("reciprocal_nmax must be an integer in [1, 20]")
    box = record["box"]
    if (box["origin_angstrom"] != [0.0, 0.0, 0.0]
            or any(v is not True for v in box["periodic"])):
        raise ValueError("the periodic export requires a zero-origin fully periodic box")
    lengths = tuple(_finite(v, "box length", positive=True) for v in box["lengths_angstrom"])
    if len(lengths) != 3 or cutoff >= min(lengths)/2:
        raise ValueError("cutoff must be strictly smaller than half the shortest box length")
    atoms = record["atoms"]
    if not atoms or abs(math.fsum(a["charge_e"] for a in atoms)) > 1e-12:
        raise ValueError("the periodic Ewald model requires a nonempty neutral cell")
    coordinates = {a["cell_atom_id"]: a["coordinates_angstrom"] for a in atoms}
    # Every short topology pair must refer to its bonded molecule image, not a
    # different nearest image. This also makes the real-space Ewald special-pair
    # correction complete inside the chosen cutoff.
    pairs = record["special_pairs"] if record["schema_version"] == "2" else record["nonbonded_pairs"]
    for pair in pairs:
        distance = pair["topological_distance"]
        if not pair["same_molecule"] or distance is None or distance > 3:
            continue
        first, second = (coordinates[a] for a in pair["atom_ids"])
        delta = tuple(x-y for x, y in zip(first, second))
        if any(not math.isfinite(v) or abs(v) >= length/2 for v, length in zip(delta, lengths)):
            raise ValueError("a short bonded-topology pair has an unsupported periodic image")
        if math.hypot(*delta) >= cutoff:
            raise ValueError("all excluded/scaled short topology pairs must lie inside the real-space cutoff")
    return cutoff, alpha, nmax, lengths


def build_periodic_alkane_setup(assigned, *, cutoff_angstrom=10.0,
                                ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    """Return scientific setup and data, with no thermo/run/output commands.

    Protocol builders can append their own bounded operations without slicing
    arbitrary single-point text. No files are written and no engine is run.
    """
    record = _record(assigned)
    cutoff, alpha, nmax, lengths = _settings(record, cutoff_angstrom,
                                           ewald_alpha_per_angstrom, reciprocal_nmax)
    parameters = ParameterSet.from_mapping(record["parameter_set"])
    if parameters.nonbonded.to_dict() != {"mixing_rule": "geometric",
            "excluded_bonded_distances": [1, 2],
            "one_four": {"lj_scale": .5, "coulomb_scale": .5, "handling": "scaled_pair_once"}}:
        raise ValueError("periodic export requires the exact supported OPLS pair conventions")
    definitions = {kind: {p.type_id: p for p in parameters.definitions if p.kind == kind}
                   for kind in ("atoms", *_KINDS, "impropers")}
    if record["interactions"]["impropers"]:
        raise ValueError("impropers are outside the supported n-alkane model")
    type_maps = {kind: {name: i+1 for i, name in enumerate(sorted(rows))}
                 for kind, rows in definitions.items()}
    atoms = record["atoms"]
    ids = {a["cell_atom_id"]: i+1 for i, a in enumerate(atoms)}
    lines = ["PolyChemMD periodic OPLS-AA Ewald numerical model", "", f"{len(atoms)} atoms"]
    for kind in _KINDS:
        lines.append(f"{len(record['interactions'][kind])} {kind}")
    lines.extend(["0 impropers", "", f"{len(type_maps['atoms'])} atom types"])
    for kind in _KINDS:
        singular = {"bonds": "bond", "angles": "angle", "dihedrals": "dihedral"}[kind]
        lines.append(f"{len(type_maps[kind])} {singular} types")
    lines.extend(["0 improper types", ""])
    for axis, length in zip("xyz", lengths):
        lines.append(f"0 {_format(length)} {axis}lo {axis}hi")
    sections = [
        ("Masses", "atoms", lambda p: [p["mass_g_mol"]]),
        ("Pair Coeffs # lj/cut/coul/long", "atoms", lambda p: [p["epsilon_kj_mol"]/KCAL_TO_KJ, p["sigma_angstrom"]]),
        ("Bond Coeffs # harmonic", "bonds", lambda p: [p["k_kj_mol_angstrom2"]/KCAL_TO_KJ, p["r0_angstrom"]]),
        ("Angle Coeffs # harmonic", "angles", lambda p: [p["k_kj_mol_rad2"]/KCAL_TO_KJ, math.degrees(p["theta0_rad"])]),
        ("Dihedral Coeffs # opls", "dihedrals", _opls_coefficients),
    ]
    for title, kind, converter in sections:
        lines.extend(["", title, ""])
        for name, number in type_maps[kind].items():
            values = converter(definitions[kind][name].values)
            lines.append(f"{number} " + " ".join(_format(v) for v in values) + f" # {name}")
    lines.extend(["", "Atoms # full", ""])
    for atom in atoms:
        images = atom["image_flags"]
        if len(images) != 3 or any(type(v) is not int for v in images):
            raise ValueError("periodic atom image flags must be three integers")
        points = atom["wrapped_coordinates_angstrom"]
        if len(points) != 3 or any(not 0 <= x < length for x, length in zip(points, lengths)):
            raise ValueError("wrapped atom positions must be inside the periodic cell")
        lines.append(f"{ids[atom['cell_atom_id']]} {atom['molecule_index']+1} "
                     f"{type_maps['atoms'][atom['type_id']]} {_format(atom['charge_e'])} "
                     + " ".join(_format(v) for v in points) + " "
                     + " ".join(str(v) for v in images))
    for kind in _KINDS:
        lines.extend(["", kind.capitalize(), ""])
        for number, row in enumerate(record["interactions"][kind], 1):
            lines.append(f"{number} {type_maps[kind][row['type_id']]} "
                         + " ".join(str(ids[a]) for a in row["atom_ids"]))
    data = "\n".join(lines) + "\n"
    setup = "\n".join((
        "# Periodic assigned OPLS-AA model with explicit finite Ewald settings.",
        "clear", "units real", "atom_style full", "boundary p p p",
        "bond_style harmonic", "angle_style harmonic", "dihedral_style opls",
        f"pair_style lj/cut/coul/long {_format(cutoff)}",
        "pair_modify mix geometric shift no tail no table 0", "dielectric 1.0",
        "special_bonds lj/coul 0.0 0.0 0.5", "read_data system.data",
        "kspace_style ewald 1e-8",
        f"kspace_modify gewald {_format(alpha)} kmax/ewald {nmax} {nmax} {nmax}",
        "neighbor 2.0 bin", "neigh_modify every 1 delay 0 check yes", "",
    ))
    return setup, data


def build_periodic_alkane_input(assigned, *, cutoff_angstrom=10.0,
                               ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    """Return a complete run-zero input script and periodic data file."""
    setup, data = build_periodic_alkane_setup(assigned, cutoff_angstrom=cutoff_angstrom,
                    ewald_alpha_per_angstrom=ewald_alpha_per_angstrom,
                    reciprocal_nmax=reciprocal_nmax)
    return _singlepoint_script(setup), data


def _singlepoint_script(setup, *, native_readback=False):
    commands = [
        "thermo 1", "thermo_style custom step temp ebond eangle edihed evdwl ecoul elong pe",
        "thermo_modify norm no format float %.17g", "run 0",
        "write_dump all custom forces.dump id mol x y z ix iy iz fx fy fz modify sort id format float %.17g",
    ]
    if native_readback:
        commands.append("write_data singlepoint.data")
    commands.extend(['print "POLYCHEMMD_RUN_COMPLETE"', ""])
    return setup + "\n".join(commands)


def prepare_periodic_alkane_singlepoint(assigned, *, run_directory, executable=None,
        timeout_seconds=30.0, cutoff_angstrom=10.0,
        ewald_alpha_per_angstrom=.30, reciprocal_nmax=16, native_readback=False):
    if type(native_readback) is not bool:
        raise ValueError("native_readback must be a Boolean")
    setup, data = build_periodic_alkane_setup(assigned, cutoff_angstrom=cutoff_angstrom,
                    ewald_alpha_per_angstrom=ewald_alpha_per_angstrom,
                    reciprocal_nmax=reciprocal_nmax)
    script = _singlepoint_script(setup, native_readback=native_readback)
    expected = ("log.lammps", "forces.dump", "singlepoint.data") if native_readback else ("log.lammps", "forces.dump")
    return LammpsRunSpec(input_script=script, run_directory=run_directory,
        executable=executable, timeout_seconds=timeout_seconds,
        expected_artifacts=expected, required_thermo_columns=_THERMO,
        input_files=(("system.data", data),), units_style="real")


def parse_periodic_alkane_force_dump(text, expected_atoms, box_lengths_angstrom):
    """Parse one step-zero periodic frame with exact molecule/image identity.

    The returned forces are raw LAMMPS real units, ordered by engine atom ID.
    """
    if not isinstance(text, str) or not isinstance(expected_atoms, (tuple, list)) or not expected_atoms:
        raise ValueError("a text dump and nonempty expected atom records are required")
    box = tuple(_finite(v, "expected box length", positive=True) for v in box_lengths_angstrom)
    if len(box) != 3:
        raise ValueError("expected periodic box must have three lengths")
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    count = len(expected_atoms)
    if (len(lines) != count+9 or lines[:3] != ["ITEM: TIMESTEP", "0", "ITEM: NUMBER OF ATOMS"]
            or lines[3] != str(count) or lines[4] != "ITEM: BOX BOUNDS pp pp pp"
            or lines[8] != "ITEM: ATOMS id mol x y z ix iy iz fx fy fz"):
        raise ValueError("periodic force dump must be one exact step-zero frame with the requested fields")
    for line, length in zip(lines[5:8], box):
        fields = line.split()
        if len(fields) != 2:
            raise ValueError("unsupported periodic force-dump box representation")
        try:
            lower, upper = (float(v) for v in fields)
        except (ValueError, OverflowError) as exc:
            raise ValueError("nonnumeric periodic dump box") from exc
        if not math.isclose(lower, 0.0, rel_tol=0, abs_tol=1e-10) or not math.isclose(upper, length, rel_tol=1e-12, abs_tol=1e-10):
            raise ValueError("periodic dump box differs from the assigned cell")
    rows = {}
    integer = re.compile(r"-?(?:0|[1-9][0-9]*)")
    for line in lines[9:]:
        fields = line.split()
        if len(fields) != 11 or any(integer.fullmatch(fields[i]) is None for i in (0, 1, 5, 6, 7)):
            raise ValueError("malformed periodic force-dump atom row")
        identity, molecule = int(fields[0]), int(fields[1])
        if not 1 <= identity <= count or identity in rows:
            raise ValueError("duplicate or out-of-range periodic engine atom ID")
        expected = expected_atoms[identity-1]
        if molecule != expected["molecule_index"]+1 or tuple(int(v) for v in fields[5:8]) != tuple(expected["image_flags"]):
            raise ValueError("periodic dump molecule or image identity changed")
        try:
            xyz = tuple(float(v) for v in fields[2:5])
            forces = tuple(float(v) for v in fields[8:11])
        except (ValueError, OverflowError) as exc:
            raise ValueError("nonnumeric periodic position or force") from exc
        if not all(math.isfinite(v) for v in (*xyz, *forces)):
            raise ValueError("nonfinite periodic position or force")
        if any(not math.isclose(a, b, rel_tol=1e-12, abs_tol=1e-10)
               for a, b in zip(xyz, expected["wrapped_coordinates_angstrom"])):
            raise ValueError("periodic dump coordinates differ from the assigned wrapped positions")
        rows[identity] = forces
    if set(rows) != set(range(1, count+1)):
        raise ValueError("periodic dump lost one or more atoms")
    return tuple(rows[i] for i in range(1, count+1))


class PeriodicSinglepointError(LammpsError):
    """An actual single-point process or its output identity/unit audit failed."""


def run_periodic_alkane_singlepoint(assigned, *, run_directory, executable=None,
        timeout_seconds=30.0, cutoff_angstrom=10.0,
        ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    """Execute a bounded periodic single point without claiming reference parity.

    This accepts the compact representation for larger cells. It does not call
    the deliberately small independent Ewald reference. Scientific sensitivity
    comparisons can consume these energies/forces, but a completed run alone
    establishes neither numerical convergence nor material-property accuracy.
    """
    timeout = _finite(timeout_seconds, "timeout_seconds", positive=True)
    if timeout > 120:
        raise ValueError("a bounded periodic single point permits at most 120 seconds")
    record = _record(assigned)
    spec = prepare_periodic_alkane_singlepoint(assigned, run_directory=run_directory,
        executable=executable, timeout_seconds=timeout, cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom, reciprocal_nmax=reciprocal_nmax,
        native_readback=True)
    report = {
        "schema_version": "1", "kind": "periodic_alkane_engine_singlepoint",
        "status": "failed", "evidence_level": "engine_execution_and_state_audit",
        "model_id": record["model_id"], "assignment_schema_version": record["schema_version"],
        "assignment_canonical_sha256": hashlib.sha256(json.dumps(record,
            sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
        "input_files": {"input.lammps_sha256": hashlib.sha256(spec.input_script.encode()).hexdigest(),
            "system.data_sha256": hashlib.sha256(dict(spec.input_files)["system.data"].encode()).hexdigest()},
        "ewald_settings": {"cutoff_angstrom": float(cutoff_angstrom),
            "alpha_per_angstrom": float(ewald_alpha_per_angstrom), "reciprocal_nmax": reciprocal_nmax,
            "relative_engine_accuracy_request": 1e-8, "surface_boundary": "tin_foil",
            "realspace_coulomb_table": False, "lj_tail_correction": False},
        "atom_ids": [atom["cell_atom_id"] for atom in record["atoms"]],
        "atom_count": len(record["atoms"]), "molecule_count": len(record["molecules"]),
        "simulation_performed": False, "independent_reference_status": "not_performed",
        "native_parameter_readback_status": "not_performed_by_singlepoint_adapter",
        "parameter_identity_basis": "validated frozen assignment and exact generated-input hashes; native coefficient/charge read-back is a separate audit",
        "numerical_convergence_status": "not_established",
        "md_integration_steps": 0, "trajectory_production_performed": False,
        "equilibrium_status": "not_assessed", "bulk_production_model_validated": False,
        "material_property_validation": "not_performed",
        "energy_units": "kJ/mol", "force_units": "kJ/(mol angstrom)",
    }

    def persist():
        (spec.run_directory/"periodic_singlepoint_record.json").write_text(
            json.dumps(report, indent=2, allow_nan=False)+"\n", encoding="utf-8")

    try:
        run = run_lammps(spec).to_dict()
    except LammpsError as exc:
        report.update(failure_phase="engine", failure_reason=str(exc), engine_run=exc.record,
            simulation_performed=bool(exc.record and exc.record.get("simulation_performed")))
        if exc.record is not None and spec.run_directory.is_dir():
            persist()
        raise PeriodicSinglepointError(str(exc), record=report) from exc
    report.update(engine_run=run, simulation_performed=True)
    try:
        blocks = run["thermo_blocks"]
        if (len(blocks) != 1 or tuple(blocks[0]["columns"]) != _THERMO
                or len(blocks[0]["rows"]) != 1):
            raise ValueError("periodic single point requires exactly the prescribed thermo block and row")
        row = blocks[0]["rows"][0]
        if len(row) != len(_THERMO) or not all(math.isfinite(v) for v in row):
            raise ValueError("periodic single-point thermo row is malformed or nonfinite")
        thermo = dict(zip(_THERMO, row))
        if thermo["Step"] != 0:
            raise ValueError("periodic single point must execute step zero")
        energies = {name: thermo[column]*KCAL_TO_KJ for name,column in _COMPONENTS.items()}
        raw = parse_periodic_alkane_force_dump((spec.run_directory/"forces.dump").read_text(),
            record["atoms"], record["box"]["lengths_angstrom"])
        forces = [[value*KCAL_TO_KJ for value in force] for force in raw]
        if not all(math.isfinite(v) for v in [*energies.values(), *(v for force in forces for v in force)]):
            raise ValueError("canonical single-point energies or forces exceed the finite numeric range")
        component_sum = math.fsum(value for name,value in energies.items() if name != "total")
        if not math.isclose(component_sum, energies["total"], rel_tol=1e-12, abs_tol=1e-8):
            raise ValueError("single-point energy components do not sum to the actual potential energy")
        header = (spec.run_directory/"singlepoint.data").read_text().splitlines()[0]
        if "units = real" not in header:
            raise ValueError("actual native data header does not confirm LAMMPS real units")
        report.update(status="completed_singlepoint", energy_components_kj_mol=energies,
            forces_kj_mol_angstrom=forces,
            atom_coordinate_molecule_image_identity="passed",
            independent_output_unit_check={"basis": "actual LAMMPS write_data header",
                "header": header, "confirmed_units_style": "real"},
            energy_component_sum_abs_difference_kj_mol=abs(component_sum-energies["total"]))
    except (OSError, ValueError, KeyError, TypeError, OverflowError, IndexError) as exc:
        report.update(failure_phase="singlepoint_output_audit", failure_reason=str(exc))
        persist()
        raise PeriodicSinglepointError(str(exc), record=report) from exc
    persist()
    return report


def verify_periodic_alkane_singlepoint(assigned, *, run_directory, executable=None,
        timeout_seconds=30.0, cutoff_angstrom=10.0, ewald_alpha_per_angstrom=.30,
        reciprocal_nmax=16, energy_abs_tol_kj_mol=1e-4,
        force_abs_tol_kj_mol_angstrom=1e-4):
    """Compare actual periodic LAMMPS energies and forces with the reference."""
    from ..forcefields.periodic_reference import evaluate_periodic_alkanes
    e_tol = _finite(energy_abs_tol_kj_mol, "energy tolerance", positive=True)
    f_tol = _finite(force_abs_tol_kj_mol_angstrom, "force tolerance", positive=True)
    record = _record(assigned)
    reference = evaluate_periodic_alkanes(assigned, cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom, reciprocal_nmax=reciprocal_nmax)
    spec = prepare_periodic_alkane_singlepoint(assigned, run_directory=run_directory,
        executable=executable, timeout_seconds=timeout_seconds, cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom, reciprocal_nmax=reciprocal_nmax)
    run = run_lammps(spec).to_dict()
    atoms = record["atoms"]
    geometry = SimpleNamespace(atoms=[SimpleNamespace(atom_id=a["cell_atom_id"]) for a in atoms],
                               coordinates_angstrom=[a["coordinates_angstrom"] for a in atoms])
    report = {
        "schema_version": "1", "status": "failed", "evidence_level": "periodic_engine_energy_force_parity",
        "model_id": record["model_id"], "boundary": "orthorhombic_3d_periodic",
        "ewald_settings": {"cutoff_angstrom": float(cutoff_angstrom),
            "alpha_per_angstrom": float(ewald_alpha_per_angstrom), "reciprocal_nmax": reciprocal_nmax,
            "relative_engine_accuracy_request": 1e-8, "surface_boundary": "tin_foil",
            "realspace_coulomb_table": False, "lj_tail_correction": False,
            "reciprocal_truncation_note": "reference uses a complete integer cube; pinned LAMMPS Ewald screens a sphere. Matching nmax does not promise identical finite mode lists; declared parity tolerances still apply."},
        "tolerances": {"energy_abs_kj_mol": e_tol, "force_component_abs_kj_mol_angstrom": f_tol},
        "reference": reference, "engine_run": run,
        "atom_ids": [a["cell_atom_id"] for a in atoms],
        "engine_atom_identity": [{"engine_atom_id": i+1, "cell_atom_id": a["cell_atom_id"],
            "source_atom_id": a["source_atom_id"], "source_atom_map": a["source_atom_map"],
            "molecule_index": a["molecule_index"]} for i, a in enumerate(atoms)],
        "backend_numerical_diagnostics": _opls_backend_diagnostics(record, geometry, run.get("engine_version_banner")),
        "realspace_backend_numerical_context": {
            "pinned_engine": "LAMMPS 22 Jul 2025 Update 4, pair_style lj/cut/coul/long",
            "source_url": "https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/KSPACE/pair_lj_cut_coul_long.cpp#L136-L141",
            "constants_source_url": "https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/ewald_const.h#L19-L25",
            "engine_realspace_function": "pinned polynomial erfc approximation with rounded force constant; pair table disabled",
            "reference_realspace_function": "Python math.erfc and 2/sqrt(pi)",
            "scope": "observed in the cited source, not guaranteed for other engine builds; finite numerical residuals remain subject to the predeclared tolerances",
            "reference_or_potential_modified": False, "parity_tolerances_modified": False,
        },
        "assignment_canonical_sha256": hashlib.sha256(json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False).encode()).hexdigest(),
        "source_snapshot": record["source_snapshot"],
        "input_files": {"input.lammps_sha256": hashlib.sha256(spec.input_script.encode()).hexdigest(),
                       "system.data_sha256": hashlib.sha256(dict(spec.input_files)["system.data"].encode()).hexdigest()},
        "md_integration_steps": 0, "trajectory_production_performed": False,
        "equilibrium_status": "not_assessed", "bulk_production_model_validated": False,
        "material_property_validation": "not_performed",
    }
    try:
        blocks = run["thermo_blocks"]
        if len(blocks) != 1 or len(blocks[0]["rows"]) != 1:
            raise ValueError("periodic single point requires one thermo block and row")
        thermo = dict(zip(blocks[0]["columns"], blocks[0]["rows"][0]))
        if thermo["Step"] != 0:
            raise ValueError("periodic single-point step must be zero")
        energies = {name: thermo[column]*KCAL_TO_KJ for name, column in _COMPONENTS.items()}
        raw = parse_periodic_alkane_force_dump((spec.run_directory/"forces.dump").read_text(), atoms, record["box"]["lengths_angstrom"])
        forces = [[v*KCAL_TO_KJ for v in row] for row in raw]
        expected_forces = reference["forces_kj_mol_angstrom"]
        if len(expected_forces) != len(forces) or any(len(row) != 3 for row in expected_forces):
            raise ValueError("periodic reference forces have incompatible shape")
        errors = {name: abs(value-reference["energy_components_kj_mol"][name]) for name, value in energies.items()}
        force_errors = [abs(a-b) for first, second in zip(forces, expected_forces) for a, b in zip(first, second)]
        if not all(math.isfinite(v) for v in (*errors.values(), *force_errors)):
            raise ValueError("periodic numerical differences must be finite")
        maximum = max(force_errors)
        rms = maximum*math.sqrt(math.fsum((v/maximum)**2 for v in force_errors)/len(force_errors)) if maximum else 0.0
        report.update(engine_energy_components_kj_mol=energies, engine_forces_kj_mol_angstrom=forces,
            energy_component_abs_errors_kj_mol=errors, max_force_component_abs_error_kj_mol_angstrom=maximum,
            rms_force_component_error_kj_mol_angstrom=rms, atom_coordinate_molecule_image_identity="passed")
        if max(errors.values()) > e_tol or maximum > f_tol:
            raise ValueError("actual periodic energy/force differences exceed predeclared tolerances")
        report["status"] = "passed"
    except (ValueError, KeyError, TypeError, OSError) as exc:
        report["failure_reason"] = str(exc)
        (spec.run_directory/"periodic_parity_record.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
        raise AlkaneParityError(str(exc), record=report) from exc
    (spec.run_directory/"periodic_parity_record.json").write_text(json.dumps(report, indent=2, allow_nan=False)+"\n")
    return report
