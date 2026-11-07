"""Bounded static virial-pressure and energy evaluations for periodic alkanes.

The explicit pressure compute excludes the kinetic term. A single constructed
configuration, even with zero velocities, is not a thermodynamic pressure or
density prediction. Existing periodic single-point/parity scripts are unchanged.
"""

import hashlib
import json
import math
import re

from ..analysis.density import density_from_mass_volume
from .alkane import KCAL_TO_KJ, _finite
from .lammps import LammpsError, LammpsRunSpec, plan_lammps, run_lammps
from .periodic_alkane import _record, build_periodic_alkane_setup


BAR_PER_ATM = 1.01325
_NATIVE_REAL_DENSITY_FACTOR = 1.0/0.602214129
_THERMO = ("Step", "Temp", "KinEng", "Press", "Pxx", "Pyy", "Pzz", "Pxy", "Pxz", "Pyz",
           "Volume", "Density", "E_bond", "E_angle", "E_dihed", "E_vdwl", "E_coul", "E_long", "PotEng")
_DUMP = ("id", "mol", "q", "x", "y", "z", "ix", "iy", "iz", "fx", "fy", "fz", "vx", "vy", "vz")
_COMPONENTS = {"bond":"E_bond", "angle":"E_angle", "dihedral":"E_dihed", "lj":"E_vdwl",
               "coulomb_realspace":"E_coul", "coulomb_kspace":"E_long", "total":"PotEng"}
_GEOMETRY_ATOM_FIELDS = ("cell_atom_id", "molecule_index", "charge_e", "coordinates_angstrom",
                         "wrapped_coordinates_angstrom", "image_flags")


class PeriodicThermoError(LammpsError):
    """The engine or current-run static thermodynamic output audit failed."""


def _canonical(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _sha(raw):
    return hashlib.sha256(raw).hexdigest()


def _geometry_digest(record):
    state = {"box":record["box"], "atoms":[{key:atom[key] for key in _GEOMETRY_ATOM_FIELDS}
                                                for atom in record["atoms"]]}
    return _sha(_canonical(state))


def _prepare(assigned, *, run_directory, executable=None, timeout_seconds=60.,
             cutoff_angstrom=10., ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    timeout = _finite(timeout_seconds, "timeout_seconds", positive=True)
    if timeout > 120:
        raise ValueError("a bounded static thermo single point permits at most 120 seconds")
    record = _record(assigned)
    setup, data = build_periodic_alkane_setup(assigned, cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom, reciprocal_nmax=reciprocal_nmax)
    script = setup + "\n".join((
        "# Static full virial: no kinetic contribution or integration.",
        "velocity all set 0.0 0.0 0.0",
        "compute polychemmd_virial all pressure NULL virial",
        "thermo 1",
        "thermo_style custom step temp ke press pxx pyy pzz pxy pxz pyz vol density ebond eangle edihed evdwl ecoul elong pe",
        "thermo_modify press polychemmd_virial norm no format float %.17g",
        "run 0",
        "write_dump all custom forces.dump id mol q x y z ix iy iz fx fy fz vx vy vz modify sort id format float %.17g",
        "write_data singlepoint.data", 'print "POLYCHEMMD_RUN_COMPLETE"', ""))
    spec = LammpsRunSpec(input_script=script, run_directory=run_directory,
        executable=executable, timeout_seconds=timeout,
        expected_artifacts=("log.lammps", "forces.dump", "singlepoint.data"),
        required_thermo_columns=_THERMO, input_files=(("system.data",data),), units_style="real")
    plan = {
        "schema_version":"1", "kind":"periodic_alkane_static_thermo_singlepoint",
        "status":"planned", "simulation_performed":False, "evidence_level":"specified_engine_evaluation",
        "model_id":record["model_id"], "assignment_schema_version":record["schema_version"],
        "assignment_canonical_sha256":_sha(_canonical(record)),
        "geometry_canonical_sha256":_geometry_digest(record),
        "geometry_fingerprint_definition":{"box":"complete assignment box record",
            "ordered_atom_fields":list(_GEOMETRY_ATOM_FIELDS),
            "serialization":"UTF-8 JSON, sorted keys, compact separators, finite values"},
        "atom_ids":[atom["cell_atom_id"] for atom in record["atoms"]],
        "atom_count":len(record["atoms"]), "molecule_count":len(record["molecules"]),
        "ewald_settings":{"cutoff_angstrom":float(cutoff_angstrom),
            "alpha_per_angstrom":float(ewald_alpha_per_angstrom), "reciprocal_nmax":reciprocal_nmax,
            "relative_engine_accuracy_request":1e-8, "surface_boundary":"tin_foil",
            "realspace_coulomb_table":False, "lj_tail_correction":False},
        "energy_units":"kJ/mol", "force_units":"kJ/(mol angstrom)", "pressure_units":"bar",
        "volume_units":"angstrom^3", "density_units":"g/cm^3",
        "pressure_interpretation":"static_virial_not_equilibrium_average",
        "pressure_compute":"compute polychemmd_virial all pressure NULL virial",
        "pressure_tensor_convention":"symmetric 3x3 Cartesian matrix, positive-compression pressure; not stress",
        "pressure_conversion":"1 atm = 1.01325 bar",
        "initial_velocities":"will_be_zeroed_and_output_verified",
        "md_integration_steps":0, "trajectory_production_performed":False,
        "independent_reference_status":"not_performed", "numerical_convergence_status":"not_established",
        "native_parameter_readback_status":"not_performed_by_static_adapter",
        "parameter_identity_basis":"validated frozen assignment and generated-input hashes; native coefficients/charges are a separate read-back audit",
        "equilibrium_status":"not_assessed", "production_ready":False,
        "density_prediction_validation":"not_performed", "material_property_validation":"not_performed",
        "engine_plan":plan_lammps(spec), "input_script":script,
        "staged_inputs":{"system.data":data},
        "input_files":{"input.lammps_sha256":_sha(script.encode()), "system.data_sha256":_sha(data.encode())},
    }
    return spec, plan, record


def plan_periodic_alkane_thermo_singlepoint(assigned, *, run_directory, executable=None,
        timeout_seconds=60., cutoff_angstrom=10., ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    """Return the exact reviewable recipe without filesystem or engine access."""
    return _prepare(assigned, run_directory=run_directory, executable=executable,
        timeout_seconds=timeout_seconds, cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom, reciprocal_nmax=reciprocal_nmax)[1]


def parse_static_thermo_dump(text, expected_atoms, box_lengths_angstrom):
    """Audit one zero-velocity step-zero frame, returning canonical forces.

    Atom IDs/molecules/fixed charges/images and wrapped coordinates must match
    their assigned source. This parser does not reconstruct native parameters.
    """
    if not isinstance(text,str) or not isinstance(expected_atoms,(tuple,list)) or not expected_atoms:
        raise ValueError("static dump requires text and nonempty expected atoms")
    if not isinstance(box_lengths_angstrom,(tuple,list)) or len(box_lengths_angstrom)!=3:
        raise ValueError("static dump requires three expected box lengths")
    box = tuple(_finite(value,"expected box length",positive=True) for value in box_lengths_angstrom)
    lines = text.splitlines(); count = len(expected_atoms)
    if (len(lines)!=count+9 or lines[:3]!=["ITEM: TIMESTEP","0","ITEM: NUMBER OF ATOMS"]
            or lines[3]!=str(count) or lines[4]!="ITEM: BOX BOUNDS pp pp pp"
            or lines[8].split()!=["ITEM:","ATOMS",*_DUMP]):
        raise ValueError("static dump must contain exactly one complete step-zero frame with prescribed columns")
    for line,length in zip(lines[5:8],box):
        tokens=line.split()
        if len(tokens)!=2:
            raise ValueError("static dump box must be orthorhombic")
        lower,upper=map(float,tokens)
        if (not all(math.isfinite(value) for value in (lower,upper))
                or not math.isclose(lower,0.,rel_tol=0,abs_tol=1e-10)
                or not math.isclose(upper,length,rel_tol=1e-12,abs_tol=1e-10)):
            raise ValueError("static dump box differs from assigned geometry")
    result={};integer=re.compile(r"[+-]?\d+")
    for line in lines[9:]:
        tokens=line.split()
        if len(tokens)!=len(_DUMP) or any(integer.fullmatch(tokens[index]) is None for index in (0,1,6,7,8)):
            raise ValueError("malformed static dump atom row")
        identity,molecule=int(tokens[0]),int(tokens[1])
        if not 1<=identity<=count or identity in result:
            raise ValueError("static dump has duplicate or out-of-range atom identities")
        expected=expected_atoms[identity-1]
        charge=float(tokens[2]);xyz=[float(value) for value in tokens[3:6]]
        images=[int(value) for value in tokens[6:9]]
        forces=[float(value)*KCAL_TO_KJ for value in tokens[9:12]]
        velocities=[float(value) for value in tokens[12:15]]
        if not all(math.isfinite(value) for value in (charge,*xyz,*forces,*velocities)):
            raise ValueError("static dump has nonfinite atom values or unrepresentable canonical forces")
        if (molecule!=expected["molecule_index"]+1 or images!=expected["image_flags"]
                or not math.isclose(charge,expected["charge_e"],rel_tol=0,abs_tol=1e-12)):
            raise ValueError("static dump molecule, charge or image identity changed")
        if any(not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-10)
               for a,b in zip(xyz,expected["wrapped_coordinates_angstrom"])):
            raise ValueError("static dump coordinates differ from assigned source")
        if any(value!=0. for value in velocities):
            raise ValueError("static virial evaluation requires exactly zero output velocities")
        result[identity]=forces
    if len(result)!=count:
        raise ValueError("static dump lost assigned atoms")
    return [result[index] for index in range(1,count+1)]


def _audit_outputs(run, spec, record):
    blocks=run["thermo_blocks"]
    if (len(blocks)!=1 or tuple(blocks[0]["columns"])!=_THERMO or len(blocks[0]["rows"])!=1):
        raise ValueError("static thermo requires exactly one prescribed thermo row")
    values=blocks[0]["rows"][0]
    if len(values)!=len(_THERMO) or not all(math.isfinite(value) for value in values):
        raise ValueError("static thermo row is malformed or nonfinite")
    thermo=dict(zip(_THERMO,values))
    if thermo["Step"]!=0 or thermo["Temp"]!=0 or thermo["KinEng"]!=0:
        raise ValueError("static thermo must be step zero with zero temperature and kinetic energy")
    energies={key:thermo[column]*KCAL_TO_KJ for key,column in _COMPONENTS.items()}
    tensor_native=[[thermo["Pxx"],thermo["Pxy"],thermo["Pxz"]],
                   [thermo["Pxy"],thermo["Pyy"],thermo["Pyz"]],
                   [thermo["Pxz"],thermo["Pyz"],thermo["Pzz"]]]
    tensor=[[value*BAR_PER_ATM for value in row] for row in tensor_native]
    pressure=thermo["Press"]*BAR_PER_ATM
    if not all(math.isfinite(value) for value in (*energies.values(),pressure,*(v for row in tensor for v in row))):
        raise ValueError("canonical static thermo values exceed the finite numeric range")
    component_sum=math.fsum(value for key,value in energies.items() if key!="total")
    if not math.isclose(component_sum,energies["total"],rel_tol=1e-12,abs_tol=1e-8):
        raise ValueError("static energy components do not sum to potential energy")
    tensor_mean=math.fsum(tensor[index][index]/3 for index in range(3))
    if not math.isclose(pressure,tensor_mean,rel_tol=1e-12,abs_tol=1e-8):
        raise ValueError("static scalar pressure differs from pressure-tensor trace divided by three")
    volume=thermo["Volume"]
    expected_volume=math.prod(record["box"]["lengths_angstrom"])
    if not math.isfinite(volume) or volume<=0 or not math.isclose(volume,expected_volume,rel_tol=1e-12,abs_tol=1e-8):
        raise ValueError("static thermo volume differs from assigned box")
    masses=[atom["mass_g_mol"] for atom in record["atoms"]]
    native_density=math.fsum(masses)/volume*_NATIVE_REAL_DENSITY_FACTOR
    if thermo["Density"]<=0 or not math.isclose(thermo["Density"],native_density,rel_tol=1e-10,abs_tol=1e-12):
        raise ValueError("static native density disagrees with the pinned real-unit mass/volume convention")
    density=density_from_mass_volume(masses,volume)
    header=(spec.run_directory/"singlepoint.data").read_text(encoding="utf-8").splitlines()[0]
    if "units = real" not in header:
        raise ValueError("static native data header does not confirm real units")
    forces=parse_static_thermo_dump((spec.run_directory/"forces.dump").read_text(encoding="utf-8"),
                                   record["atoms"],record["box"]["lengths_angstrom"])
    return {"energy_components_kj_mol":energies,"forces_kj_mol_angstrom":forces,
        "pressure_bar":pressure,"pressure_tensor_bar":tensor,"volume_angstrom3":volume,
        "density_g_cm3":density,"engine_density_g_cm3":thermo["Density"],
        "temperature_k":0.,"kinetic_energy_kj_mol":0.,"initial_velocities":"zeroed_and_output_verified",
        "atom_coordinate_molecule_charge_image_identity":"passed",
        "energy_component_sum_abs_difference_kj_mol":abs(component_sum-energies["total"]),
        "pressure_trace_abs_difference_bar":abs(pressure-tensor_mean),
        "independent_output_unit_check":{"basis":"actual LAMMPS write_data header","header":header,
                                         "confirmed_units_style":"real"},
        "density_conversion_audit":{"canonical_basis":"exact SI Avogadro constant, 6.02214076e23 mol^-1",
            "native_real_factor":_NATIVE_REAL_DENSITY_FACTOR,"expected_engine_density_g_cm3":native_density,
            "source":"https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/update.cpp#L175"}}


def run_periodic_alkane_thermo_singlepoint(assigned, *, run_directory, executable=None,
        timeout_seconds=60., cutoff_angstrom=10., ewald_alpha_per_angstrom=.30, reciprocal_nmax=16):
    """Run one static virial evaluation, preserving explicit failures and limits."""
    spec,plan,record=_prepare(assigned,run_directory=run_directory,executable=executable,
        timeout_seconds=timeout_seconds,cutoff_angstrom=cutoff_angstrom,
        ewald_alpha_per_angstrom=ewald_alpha_per_angstrom,reciprocal_nmax=reciprocal_nmax)
    report={key:value for key,value in plan.items() if key not in ("input_script","staged_inputs")}
    report.update(status="failed",evidence_level="engine_execution_and_static_state_audit")

    def persist():
        (spec.run_directory/"periodic_thermo_record.json").write_text(
            json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")

    try:
        run=run_lammps(spec).to_dict()
    except LammpsError as exc:
        report.update(failure_phase="engine",failure_reason=str(exc),engine_run=exc.record,
            simulation_performed=bool(exc.record and exc.record.get("simulation_performed")))
        if exc.record is not None and spec.run_directory.is_dir():
            persist()
        raise PeriodicThermoError(str(exc),record=report) from exc
    report.update(engine_run=run,simulation_performed=True)
    try:
        audited=_audit_outputs(run,spec,record)
        report.update(audited,status="completed_static_thermo_singlepoint")
    except (OSError,ValueError,KeyError,TypeError,OverflowError,IndexError) as exc:
        report.update(failure_phase="static_thermo_output_audit",failure_reason=str(exc))
        persist()
        raise PeriodicThermoError(str(exc),record=report) from exc
    persist()
    return report
