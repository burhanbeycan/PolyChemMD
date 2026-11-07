"""Bounded continuation of audited native periodic-alkane NPT checkpoints.

A native restart preserves an evolved dynamical state. It must never be fed
back into the rigid starting-cell contract. Compatibility is deliberately
restricted to an unchanged engine build/platform/configuration and one process.
Successful continuation is execution evidence, not equilibrium certification.
"""

from copy import deepcopy
from collections import Counter
from dataclasses import asdict, dataclass, fields
import hashlib
import json
import math
from pathlib import Path
import re

from ..analysis.density import density_from_mass_volume
from ..engines.lammps import LammpsError, LammpsRunSpec, _finite, plan_lammps, run_lammps
from .alkane import (BAR_PER_ATM, KCAL_TO_KJ, PeriodicAlkaneProtocolSpec,
    PeriodicProtocolError, _PINNED_REAL_MASS_DENSITY_FACTOR, _THERMO_COLUMNS,
    _canonical, _integer, _positive, _sha, parse_periodic_final_dump)

_CHECKPOINT_FIELDS = {"schema_version", "kind", "model_id", "assignment_canonical_sha256",
    "artifact", "state_snapshot", "units_style", "step", "elapsed_time_fs", "atom_count",
    "box", "protocol_spec", "integrator", "producer_engine_identity", "engine_version_banner",
    "parent_checkpoint_canonical_sha256", "checkpoint_chain_index", "equilibrium_status",
    "production_ready", "checkpoint_canonical_sha256"}
_INTEGRATOR = {"fix_id":"polychemmd_npt","group_id":"all","style":"npt",
               "thermostat_damping_fs":100.0,"barostat_damping_fs":1000.0,
               "target_schedule":"constant_temperature_and_isotropic_pressure"}


@dataclass(frozen=True)
class PeriodicAlkaneContinuationSpec:
    steps: int = 100
    timeout_seconds: float = 60.0
    thermo_stride_steps: int = 10

    def __post_init__(self):
        _integer(self.steps,"steps",1,10000)
        _integer(self.thermo_stride_steps,"thermo_stride_steps",1,10000)
        timeout = _positive(self.timeout_seconds,"timeout_seconds")
        if timeout > 120:
            raise ValueError("a bounded continuation cannot exceed 120 seconds")
        object.__setattr__(self,"timeout_seconds",timeout)

    def to_dict(self):
        return asdict(self)

    @classmethod
    def from_mapping(cls,data):
        expected = {field.name for field in fields(cls)}
        if not isinstance(data,dict) or not set(data).issubset(expected):
            raise ValueError(f"continuation fields must be drawn from {sorted(expected)}")
        return cls(**data)


def _digest_record(checkpoint):
    return _sha(_canonical({k:v for k,v in checkpoint.items() if k != "checkpoint_canonical_sha256"}))


def _artifact(directory,name):
    path = directory/name
    if not path.is_file() or path.is_symlink() or not path.resolve().is_relative_to(directory.resolve()):
        raise ValueError(f"checkpoint artifact missing or unsafe: {name}")
    if not 0 < path.stat().st_size <= 64*1024*1024:
        raise ValueError(f"checkpoint artifact exceeds the bounded 64 MiB read budget: {name}")
    raw = path.read_bytes()
    if not raw:
        raise ValueError(f"checkpoint artifact is empty: {name}")
    return raw,{"file_name":name,"sha256":_sha(raw),"size_bytes":len(raw)}


def checkpoint_manifest(assignment,protocol,engine_run,directory,*,step,parent=None):
    """Create an engine-produced checkpoint manifest after current-run audits.

    The native artifact contains velocities, images and active NPT fix state;
    the separate text snapshot makes state and atom identities reviewable.
    Binary checksums are not claimed to be portable across platforms/builds.
    """
    directory = Path(directory).resolve()
    _,binary = _artifact(directory,"checkpoint.restart")
    text,snapshot = _artifact(directory,"checkpoint.dump")
    state = parse_periodic_final_dump(text.decode("utf-8"),assignment,step,include_velocities=True)
    if not engine_run.get("engine_identity"):
        raise ValueError("checkpoint requires the actual producer engine identity")
    checkpoint = {
        "schema_version":"1","kind":"periodic_alkane_npt_checkpoint",
        "model_id":assignment["model_id"],
        "assignment_canonical_sha256":_sha(_canonical(assignment)),
        "artifact":binary,"state_snapshot":snapshot,"units_style":"real",
        "step":step,"elapsed_time_fs":step*protocol.timestep_fs,
        "atom_count":len(assignment["atoms"]),"box":state["box"],
        "protocol_spec":protocol.to_dict(),"integrator":deepcopy(_INTEGRATOR),
        "producer_engine_identity":deepcopy(engine_run["engine_identity"]),
        "engine_version_banner":engine_run.get("engine_version_banner"),
        "parent_checkpoint_canonical_sha256":parent["checkpoint_canonical_sha256"] if parent else None,
        "checkpoint_chain_index":parent["checkpoint_chain_index"]+1 if parent else 0,
        "equilibrium_status":"insufficient_sampling","production_ready":False,
    }
    checkpoint["checkpoint_canonical_sha256"] = _digest_record(checkpoint)
    return checkpoint


def _validate_checkpoint(checkpoint,assignment):
    if not isinstance(checkpoint,dict) or set(checkpoint) != _CHECKPOINT_FIELDS:
        raise ValueError("checkpoint fields differ from the exact native NPT checkpoint schema")
    if (checkpoint["schema_version"] != "1" or checkpoint["kind"] != "periodic_alkane_npt_checkpoint"
            or checkpoint["units_style"] != "real" or checkpoint["integrator"] != _INTEGRATOR):
        raise ValueError("unsupported checkpoint schema, units or integrator")
    if (checkpoint["model_id"] != assignment["model_id"]
            or checkpoint["assignment_canonical_sha256"] != _sha(_canonical(assignment))
            or type(checkpoint["atom_count"]) is not int
            or checkpoint["atom_count"] != len(assignment["atoms"])):
        raise ValueError("checkpoint assignment/model/atom identity differs from the supplied source")
    if checkpoint["checkpoint_canonical_sha256"] != _digest_record(checkpoint):
        raise ValueError("checkpoint manifest checksum changed")
    _integer(checkpoint["step"],"checkpoint step",1,2**53-10001)
    _integer(checkpoint["checkpoint_chain_index"],"checkpoint_chain_index",0,1000000)
    parent = checkpoint["parent_checkpoint_canonical_sha256"]
    if (checkpoint["checkpoint_chain_index"] == 0) != (parent is None):
        raise ValueError("checkpoint parent identity/chain index is inconsistent")
    if parent is not None and (not isinstance(parent,str) or not re.fullmatch("[0-9a-f]{64}",parent)):
        raise ValueError("checkpoint parent identity must be a SHA256")
    if checkpoint["production_ready"] is not False or checkpoint["equilibrium_status"] != "insufficient_sampling":
        raise ValueError("native continuation cannot certify equilibrium or production suitability")
    protocol = PeriodicAlkaneProtocolSpec.from_mapping(checkpoint["protocol_spec"])
    if protocol.to_dict() != checkpoint["protocol_spec"]:
        raise ValueError("checkpoint must preserve all original protocol controls explicitly")
    elapsed = _positive(checkpoint["elapsed_time_fs"],"elapsed_time_fs")
    if elapsed != checkpoint["step"]*protocol.timestep_fs:
        raise ValueError("checkpoint step and physical time disagree")
    box = checkpoint["box"]
    if (not isinstance(box,dict) or set(box) != {"origin_angstrom","lengths_angstrom","periodic"}
            or any(not isinstance(box.get(field),list) or len(box[field]) != 3
                   for field in ("origin_angstrom","lengths_angstrom","periodic"))
            or any(v is not True for v in box["periodic"])):
        raise ValueError("checkpoint box must be orthorhombic and fully periodic")
    for value in [*box["origin_angstrom"],*box["lengths_angstrom"]]:
        _finite(value,"checkpoint box component")
    if min(box["lengths_angstrom"]) <= 2*protocol.cutoff_angstrom:
        raise ValueError("checkpoint box violates the unique minimum-image cutoff bound")
    identity = checkpoint["producer_engine_identity"]
    if (not isinstance(identity,dict) or not identity
            or not isinstance(checkpoint["engine_version_banner"],str) or not checkpoint["engine_version_banner"]):
        raise ValueError("checkpoint must identify its actual native engine configuration")
    for key,name in (("artifact","checkpoint.restart"),("state_snapshot","checkpoint.dump")):
        item = checkpoint[key]
        if (not isinstance(item,dict) or set(item) != {"file_name","sha256","size_bytes"}
                or item["file_name"] != name or not isinstance(item["sha256"],str)
                or not re.fullmatch("[0-9a-f]{64}",item["sha256"])):
            raise ValueError("checkpoint artifact metadata is malformed")
        _integer(item["size_bytes"],"checkpoint artifact size",1,64*1024*1024)
    return protocol


def _build(assigned,checkpoint,spec,*,run_directory,executable,binary=None):
    from ..engines.periodic_alkane import _record, _neighbor_capacity_commands
    if not isinstance(spec,PeriodicAlkaneContinuationSpec):
        raise ValueError("spec must be a PeriodicAlkaneContinuationSpec")
    assignment = _record(assigned)
    checkpoint = deepcopy(checkpoint)
    protocol = _validate_checkpoint(checkpoint,assignment)
    target = format(protocol.temperature_k,".17g")
    pressure = format(protocol.pressure_bar/BAR_PER_ATM,".17g")
    cutoff = protocol.cutoff_angstrom
    script = "\n".join((
        "# Native NPT continuation: no minimization, reseeding or timestep reset.",
        "clear","read_restart source.restart",
        "# Kspace/neighbor/fixes are not restored by read_restart; redeclare them.",
        "kspace_style ewald 1e-8",
        f"kspace_modify gewald {protocol.ewald_alpha_per_angstrom:.17g} kmax/ewald {protocol.reciprocal_nmax} {protocol.reciprocal_nmax} {protocol.reciprocal_nmax}",
        "pair_modify mix geometric shift no tail no table 0","dielectric 1.0",
        "neighbor 2.0 bin","neigh_modify every 1 delay 0 check yes",
        *_neighbor_capacity_commands(len(assignment["atoms"])),
        f"timestep {protocol.timestep_fs:.17g}",
        f"fix polychemmd_npt all npt temp {target} {target} 100.0 iso {pressure} {pressure} 1000.0",
        f'variable polychemmd_box_invalid equal "(lx<={2*cutoff:.17g})||(ly<={2*cutoff:.17g})||(lz<={2*cutoff:.17g})"',
        "fix polychemmd_boxguard all halt 1 v_polychemmd_box_invalid != 0 error hard message yes",
        f"fix polychemmd_bondguard all halt 1 bondmax >= {cutoff/3:.17g} error hard message yes",
        f"thermo {spec.thermo_stride_steps}",
        "thermo_style custom step temp press vol density pe etotal ebond eangle edihed evdwl ecoul elong lx ly lz",
        "thermo_modify norm no format float %.17g",
        'print "POLYCHEMMD_RESTART_INITIAL_BEGIN"',"run 0",
        'print "POLYCHEMMD_RESTART_INITIAL_END"',
        "write_data resumed_initial.data",
        "write_dump all custom resumed_initial.dump id mol q x y z ix iy iz fx fy fz vx vy vz modify sort id format float %.17g",
        'print "POLYCHEMMD_CONTINUATION_BEGIN"',f"run {spec.steps}",
        'print "POLYCHEMMD_CONTINUATION_END"',
        "write_restart checkpoint.restart",
        "write_dump all custom checkpoint.dump id mol q x y z ix iy iz fx fy fz vx vy vz modify sort id format float %.17g",
        "write_data prepared.data",
        'print "POLYCHEMMD_RUN_COMPLETE"',""))
    run_spec = LammpsRunSpec(script,run_directory,executable=executable,
        timeout_seconds=spec.timeout_seconds,units_style="real",
        expected_artifacts=("log.lammps","resumed_initial.data","resumed_initial.dump","prepared.data","checkpoint.restart","checkpoint.dump"),
        required_thermo_columns=_THERMO_COLUMNS,
        input_binary_files=(("source.restart",binary),) if binary is not None else (),
        expected_engine_identity=checkpoint["producer_engine_identity"])
    report = {
        "schema_version":"1","kind":"periodic_alkane_npt_continuation",
        "status":"planned","simulation_performed":False,"input_script":script,
        "model_id":assignment["model_id"],
        "assignment_canonical_sha256":checkpoint["assignment_canonical_sha256"],
        "source_checkpoint_canonical_sha256":checkpoint["checkpoint_canonical_sha256"],
        "source_checkpoint_artifact_sha256":checkpoint["artifact"]["sha256"],
        "continuation_spec":spec.to_dict(),"protocol_spec":protocol.to_dict(),
        "start_step":checkpoint["step"],"end_step_requested":checkpoint["step"]+spec.steps,
        "duration_fs_requested":spec.steps*protocol.timestep_fs,
        "physical_time_fs_requested":(checkpoint["step"]+spec.steps)*protocol.timestep_fs,
        "engine_plan":plan_lammps(run_spec),"equilibrium_status":"not_assessed",
        "production_ready":False,"density_prediction_validation":"not_performed",
        "state_interpretation":"evolved native MD state; original rigid packing remains immutable",
        "native_restart_compatibility":"same engine build/platform/configuration/direct process layout only; linked-system libraries beyond recorded native package hashes are not exhaustively fingerprinted",
        "parameter_readback_precision":{
            "binary":"unchanged native checkpoint SHA256; no coefficient rounding for dynamics",
            "text_coefficients":"pinned Pair/Bond/Angle/Dihedral write_data uses C %g; compare the six-significant-digit serialization",
            "angle_writer_source":"https://github.com/lammps/lammps/blob/stable_22Jul2025_update4/src/MOLECULE/angle_harmonic.cpp",
            "restored_force_component_abs_tolerance_kj_mol_angstrom":1e-7},
    }
    return run_spec,report,assignment,protocol,checkpoint


def plan_periodic_alkane_continuation(assigned,checkpoint,spec,*,checkpoint_directory,
                                     run_directory,executable=None):
    """Validate the recipe/manifest without reading files or invoking an engine.

    Binary presence/hash and producer-engine matching are execution-time gates;
    a plan therefore explicitly reports them as not yet checked.
    """
    result = _build(assigned,checkpoint,spec,run_directory=run_directory,executable=executable)[1]
    result["checkpoint_artifact_validation"] = "deferred_until_execution"
    result["checkpoint_directory"] = str(Path(checkpoint_directory).expanduser().resolve())
    return result


def _sections(text):
    names = {"Masses","Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs",
             "Atoms","Velocities","Bonds","Angles","Dihedrals","Impropers"}
    sections,current = {},None
    for line in text.splitlines():
        clean = line.split("#",1)[0].strip()
        if clean in names:
            current=clean;sections[current]={};continue
        if not clean or current is None:
            continue
        tokens=clean.split()
        if not re.fullmatch(r"\d+",tokens[0]):
            raise ValueError("unsupported text data section in checkpoint audit")
        key=int(tokens[0])
        if key in sections[current]:
            raise ValueError("duplicate text data row identity in checkpoint audit")
        values=[float(v) for v in tokens[1:]]
        if not all(math.isfinite(v) for v in values):
            raise ValueError("nonfinite data coefficient/topology value")
        sections[current][key]=values
    return sections


def _audit_data_identity(text,assigned,protocol):
    from ..engines.periodic_alkane import build_periodic_alkane_setup
    if "units = real" not in text.splitlines()[0]:
        raise ValueError("native restart read-back does not confirm real units")
    _,expected = build_periodic_alkane_setup(assigned,cutoff_angstrom=protocol.cutoff_angstrom,
        ewald_alpha_per_angstrom=protocol.ewald_alpha_per_angstrom,
        reciprocal_nmax=protocol.reciprocal_nmax)
    reference,actual = _sections(expected),_sections(text)
    for name,rows in reference.items():
        if name in {"Bonds","Angles","Dihedrals"}:
            # write_data renumbers interactions in local-atom traversal order;
            # those row numbers have no molecular identity. Preserve the full
            # type/endpoint multiset, including multiplicity. Full reversal is
            # equivalent for these harmonic bonds/angles and OPLS proper terms.
            def topology_multiset(table):
                canonical=[]
                for row in table.values():
                    if any(not value.is_integer() or value < 1 for value in row):
                        raise ValueError("native interaction type/atom IDs must be positive integers")
                    endpoints=tuple(int(value) for value in row[1:])
                    canonical.append((int(row[0]),min(endpoints,endpoints[::-1])))
                return Counter(canonical)
            if name not in actual or topology_multiset(rows) != topology_multiset(actual[name]):
                raise ValueError(f"native checkpoint read-back changed {name} topology/type/multiplicity")
            continue
        if name not in actual or set(actual[name]) != set(rows):
            raise ValueError(f"native checkpoint read-back changed {name} identities")
        for key,values in rows.items():
            observed=actual[name][key]
            # The Atoms record retains id/molecule/type/charge; coordinates and
            # image flags are evolved state, checked independently in dumps.
            if name == "Atoms": values,observed=values[:3],observed[:3]
            # These pinned LAMMPS styles write their coefficients with C %g
            # (six significant digits). Audit that documented serialization,
            # rather than treating rounded text as the exact native doubles.
            # The unchanged binary is separately bound to its producer by SHA.
            if name in {"Pair Coeffs","Bond Coeffs","Angle Coeffs","Dihedral Coeffs"}:
                values=[float(format(value,".6g")) for value in values]
            if len(values) != len(observed) or any(not math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12)
                                                  for a,b in zip(values,observed)):
                raise ValueError(f"native checkpoint read-back changed {name} parameters/topology")
    if actual.get("Impropers"):
        raise ValueError("native checkpoint unexpectedly contains improper topology")


def _state_equal(source,restored):
    if source["step"] != restored["step"] or source["atom_count"] != restored["atom_count"]:
        raise ValueError("native restart changed the saved step/atom count")
    for field in ("origin_angstrom","lengths_angstrom"):
        if any(not math.isclose(a,b,rel_tol=0,abs_tol=1e-11)
               for a,b in zip(source["box"][field],restored["box"][field])):
            raise ValueError("native restart changed the saved box")
    for before,after in zip(source["atoms"],restored["atoms"]):
        for field in ("engine_id","cell_atom_id","molecule_index","charge_e","image_flags"):
            if before[field] != after[field]:
                raise ValueError("native restart changed an atom identity/charge/image")
        for field in ("wrapped_coordinates_angstrom","coordinates_angstrom","velocities_angstrom_fs"):
            if any(not math.isclose(a,b,rel_tol=0,abs_tol=1e-11) for a,b in zip(before[field],after[field])):
                raise ValueError("native restart changed saved coordinates/velocities")
        # The checkpoint snapshot retains full-precision current forces. This
        # independently checks the native parameter behavior beyond the rounded
        # %g coefficients in write_data. This tolerance is specified before the
        # successful actual restart comparison, in canonical kJ/(mol Angstrom).
        if any(not math.isclose(a,b,rel_tol=0,abs_tol=1e-7)
               for a,b in zip(before["forces_kj_mol_angstrom"],after["forces_kj_mol_angstrom"])):
            raise ValueError("native restart changed saved forces beyond 1e-7 kJ/(mol Angstrom)")


def _audit_continuation(run,assigned,assignment,protocol,checkpoint,spec,directory,source_state):
    blocks=run["thermo_blocks"]
    if len(blocks) != 2 or any(tuple(b["columns"]) != _THERMO_COLUMNS for b in blocks):
        raise ValueError("continuation requires one run-zero restore audit and one integration block")
    start,end=checkpoint["step"],checkpoint["step"]+spec.steps
    if (len(blocks[0]["rows"]) != 1 or blocks[0]["rows"][0][0] != start
            or blocks[1]["rows"][0][0] != start or blocks[1]["rows"][-1][0] != end):
        raise ValueError("native continuation failed to preserve/advance exactly its requested timestep")
    log=(directory/"log.lammps").read_text(encoding="utf-8")
    if not re.search(r"Resetting global fix info from restart file:\s*fix style: npt, fix ID: polychemmd_npt",log):
        raise ValueError("LAMMPS did not confirm restoring the active NPT fix state")
    if re.search(r"Unused restart file global fix info",log):
        raise ValueError("LAMMPS discarded saved global fix state")
    for name in ("resumed_initial.data","prepared.data"):
        _audit_data_identity((directory/name).read_text(encoding="utf-8"),assigned,protocol)
    initial=parse_periodic_final_dump((directory/"resumed_initial.dump").read_text(),assignment,start,include_velocities=True)
    _state_equal(source_state,initial)
    final=parse_periodic_final_dump((directory/"checkpoint.dump").read_text(),assignment,end,include_velocities=True)
    rows=[]
    masses=[atom["mass_g_mol"] for atom in assignment["atoms"]]
    for block in blocks:
        for values in block["rows"]:
            row=dict(zip(block["columns"],values))
            if (not all(math.isfinite(v) for v in values) or row["Volume"] <= 0
                    or min(row[axis] for axis in ("Lx","Ly","Lz")) <= 2*protocol.cutoff_angstrom):
                raise ValueError("continuation thermo is nonfinite or outside the supported box regime")
            native=math.fsum(masses)/row["Volume"]*_PINNED_REAL_MASS_DENSITY_FACTOR
            if not math.isclose(row["Density"],native,rel_tol=1e-10,abs_tol=1e-12):
                raise ValueError("continuation thermo density violates the pinned real-unit convention")
            if block is blocks[-1]:
                rows.append({"step":int(row["Step"]),"time_fs":row["Step"]*protocol.timestep_fs,
                         "temperature_k":row["Temp"],"pressure_bar":row["Press"]*BAR_PER_ATM,
                         "volume_angstrom3":row["Volume"],
                         "density_g_cm3":density_from_mass_volume(masses,row["Volume"]),
                             "engine_density_g_cm3":row["Density"]})
    last=dict(zip(blocks[-1]["columns"],blocks[-1]["rows"][-1]))
    volume=math.prod(final["box"]["lengths_angstrom"])
    if not math.isclose(volume,last["Volume"],rel_tol=1e-10,abs_tol=1e-8):
        raise ValueError("continuation final volume disagrees between dump and thermo")
    final_by_id={atom["cell_atom_id"]:atom for atom in final["atoms"]}
    for bond in assignment["interactions"]["bonds"]:
        length=math.dist(*(final_by_id[key]["coordinates_angstrom"] for key in bond["atom_ids"]))
        if not math.isfinite(length) or not 0 < length < protocol.cutoff_angstrom/3:
            raise ValueError("continuation final bond image violates the supported cutoff bound")
    final.update(temperature_k=last["Temp"],pressure_bar=last["Press"]*BAR_PER_ATM,
                 volume_angstrom3=volume,density_g_cm3=density_from_mass_volume(masses,volume),
                 engine_density_g_cm3=last["Density"],potential_energy_kj_mol=last["PotEng"]*KCAL_TO_KJ,
                 total_energy_kj_mol=last["TotEng"]*KCAL_TO_KJ)
    return final,rows


def _write_report(directory,report):
    (directory/"continuation_record.json").write_text(json.dumps(report,indent=2,allow_nan=False)+"\n",encoding="utf-8")


def run_periodic_alkane_continuation(assigned,checkpoint,spec,*,checkpoint_directory,
                                    run_directory,executable=None):
    """Hash-check, stage and continue a native checkpoint within finite resources."""
    # Manifest/assignment validation precedes reading or writing any artifacts.
    _,plan,assignment,protocol,checkpoint=_build(assigned,checkpoint,spec,
        run_directory=run_directory,executable=executable)
    source=Path(checkpoint_directory).expanduser().resolve()
    binary,metadata=_artifact(source,"checkpoint.restart")
    text,snapshot_metadata=_artifact(source,"checkpoint.dump")
    if metadata != checkpoint["artifact"] or snapshot_metadata != checkpoint["state_snapshot"]:
        raise ValueError("checkpoint binary or reviewable state snapshot checksum changed")
    source_state=parse_periodic_final_dump(text.decode("utf-8"),assignment,checkpoint["step"],include_velocities=True)
    if source_state["box"] != checkpoint["box"]:
        raise ValueError("checkpoint snapshot and manifest box disagree")
    run_spec,report,_,_,_=_build(assigned,checkpoint,spec,run_directory=run_directory,
        executable=executable,binary=binary)
    report.update(status="failed",checkpoint_artifact_validation="passed")
    try:
        run=run_lammps(run_spec).to_dict()
    except LammpsError as exc:
        report.update(failure_phase="engine",failure_reason=str(exc),engine_run=exc.record,
                      simulation_performed=bool(exc.record and exc.record.get("simulation_performed")))
        if exc.record is not None and run_spec.run_directory.is_dir():
            _write_report(run_spec.run_directory,report)
        raise PeriodicProtocolError(str(exc),record=report) from exc
    report.update(engine_run=run,simulation_performed=True)
    try:
        final,rows=_audit_continuation(run,assigned,assignment,protocol,checkpoint,spec,
                                       run_spec.run_directory,source_state)
        next_checkpoint=checkpoint_manifest(assignment,protocol,run,run_spec.run_directory,
             step=checkpoint["step"]+spec.steps,parent=checkpoint)
        report.update(status="completed_npt_continuation",final_state=final,
                      checkpoint=next_checkpoint,thermo_rows=rows,
                      equilibrium_status="insufficient_sampling",md_integration_steps_executed=spec.steps,
                      npt_fix_restart_state="engine_confirmed_restoration_same_fix_id",
                      topology_parameter_unit_readback_audit="passed",
                      initial_position_velocity_image_box_audit="passed")
    except (OSError,ValueError,KeyError,TypeError,OverflowError) as exc:
        report.update(failure_phase="continuation_output_audit",failure_reason=str(exc))
        _write_report(run_spec.run_directory,report)
        raise PeriodicProtocolError(str(exc),record=report) from exc
    _write_report(run_spec.run_directory,report)
    return report
