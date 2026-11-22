"""Dispatch scientifically distinct outputs. Optional engines load only on demand."""

from __future__ import annotations
from importlib import metadata, util
import sysconfig
from pathlib import Path
import platform
import shutil

DEFAULTS = dict(
    polymer="PEO",
    dp=8,
    seed=17,
    geometry="single",
    nx=6,
    ny=4,
    gap_A=22,
    water_height_A=22,
    oh_top=0.5,
    oh_bottom=0.5,
    hydrolysis=1.0,
    acetate_pattern="random",
    pattern="random",
    conformer_ff="MMFF94",
    chains=1,
    salt_M=0.15,
    water_density=0.997,
    packing_engine="builtin",
)

DESIGN_FIELDS = set(DEFAULTS) | {
    "side_dp",
    "graft_every",
    "water_count",
    "min_distance_A",
    "acetate_dyads",
    "acetate_dyad_fraction",
    "acetate_positions",
    "surface_correlation_top",
    "surface_correlation_bottom",
    "surface_seed",
    "sequence_seed",
    "tacticity_seed",
    "conformer_seed",
    "vacuum_padding_A",
}


def versions():
    result = {"python": platform.python_version()}
    for package in (
        "polychemmd",
        "rdkit",
        "ase",
        "numpy",
        "scipy",
        "moltemplate",
        "molcrafts-molpack",
        "packmol",
        "lammps",
        "MDAnalysis",
    ):
        try:
            result[package] = metadata.version(package)
        except metadata.PackageNotFoundError:
            continue
    return result


def status():
    v = versions()
    tools = {}
    for name, module in (
        ("rdkit", "rdkit"),
        ("ase", "ase"),
        ("mdanalysis", "MDAnalysis"),
        ("molpack", "molpack"),
    ):
        pkg = {"mdanalysis": "MDAnalysis", "molpack": "molcrafts-molpack"}.get(
            name, name
        )
        tools[name] = {
            "available": util.find_spec(module) is not None,
            "version": v.get(pkg),
        }
    for name, exe in (("packmol", "packmol"), ("lammps", "lmp")):
        tools[name] = {
            "available": bool(
                shutil.which(exe)
                or any(
                    (Path(sysconfig.get_path("scripts")) / (exe + ext)).is_file()
                    for ext in ("", ".exe")
                )
            ),
            "version": v.get(name),
        }
    from .radon import status as radon_status

    tools["radonpy"] = radon_status()
    tools["builtin"] = {"available": True, "version": "1"}
    tools["slit"] = {
        "available": tools["rdkit"]["available"] and tools["ase"]["available"]
    }
    tools["moltemplate"] = {
        "available": "moltemplate" in v,
        "version": v.get("moltemplate"),
    }
    from .iff_sources import status as iff_status

    return {
        "tools": tools,
        "defaults": DEFAULTS,
        "software": v,
        "profiles": iff_status(),
    }


def execute(kind, config, out, store):
    keys = {
        "geometry": DESIGN_FIELDS,
        "packing": DESIGN_FIELDS | {"packing_seed"},
        "topology": {"source_job", "cutoff_A"},
        "packing_benchmark": {"cases", "engines", "seeds", "evaluate_forces"},
        "pva_ensemble": {
            "mn_g_mol",
            "dispersity",
            "hydrolysis",
            "acetate_dyad_fraction",
            "dp_min",
            "dp_max",
            "n_chains",
            "seed",
        },
        "radonpy": {"polymer", "dp", "seed"},
        "trajectory": {
            "source_job",
            "source_job_id",
            "bin_width_A",
            "contact_cutoff_A",
        },
        "experiment": {"assay", "reference", "rows"},
    }
    if kind in keys and set(config) - keys[kind]:
        raise ValueError(
            "Unsupported recipe field(s): "
            + ", ".join(sorted(set(config) - keys[kind]))
        )
    if kind in {"geometry", "packing"}:
        from .geometry import build_geometry

        full = {**DEFAULTS, **config}
        if config.get("polymer") == "bare" and "chains" not in config:
            full["chains"] = 0
        result = build_geometry(full, out)
        if kind == "geometry":
            summary = {
                k: v
                for k, v in result.items()
                if k not in {"atoms", "templates", "bonds", "polymer_bonds"}
            }
            summary["atom_count"] = len(result["atoms"])
            return summary
        from .packing import pack_system

        return pack_system(full, out)
    if kind == "topology":
        from .topology_build import run_topology

        source = config.get("source_job")
        if not isinstance(source, str) or store.get(source)["kind"] != "packing":
            raise ValueError("Choose a completed packing job to generate its topology.")
        return run_topology(config, out, source_dir=store.output_dir(source))
    if kind == "packing_benchmark":
        from .packing_benchmark import run_benchmark

        return run_benchmark(config, out)
    if kind == "pva_ensemble":
        from .ensemble import run_ensemble

        return run_ensemble(config, out)
    if kind == "md":
        if (
            config.get("profile") == "registered"
            or config.get("md_profile") == "registered"
        ):
            from .register_topology import run_registered

            source = config.get("source_job")
            if not isinstance(source, str) or store.get(source)["kind"] != "packing":
                raise ValueError(
                    "Choose a completed packed structure for this topology."
                )
            return run_registered(config, out, source_dir=store.output_dir(source))
        from .md import run_md

        if "restart_job" in config:
            source = config["restart_job"]
            if not isinstance(source, str) or store.get(source)["kind"] != "md":
                raise ValueError("Choose a completed checkpoint-containing MD source.")
            return run_md(config, out, source_dir=store.output_dir(source))
        return run_md(config, out)
    if kind == "radonpy":
        from .radon import run_radon

        return run_radon(config, out)
    if kind == "experiment":
        from .experiment import run_experiment

        return run_experiment(config, out)
    if kind == "trajectory":
        from .trajectory import run_trajectory

        source = config.get("source_job", config.get("source_job_id"))
        if not isinstance(source, str) or store.get(source)["kind"] != "md":
            raise ValueError("Select a successfully completed MD job.")
        return run_trajectory(config, out, source_dir=store.output_dir(source))
    raise ValueError("Unknown kind")
