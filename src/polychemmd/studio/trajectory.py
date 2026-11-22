"""Analyse a Studio MD job in the instantaneous basal-plane coordinate system.

LAMMPS restricted-triclinic bounds follow the LAMMPS dump specification:
https://docs.lammps.org/Howto_triclinic.html . Cell vectors are rows, in Å.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import itertools
import json
import math
from dataclasses import dataclass
from pathlib import Path


def _number(value, name, low, high):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be a number")
    if not math.isfinite(value) or not low <= value <= high:
        raise ValueError(f"{name} must be finite and between {low} and {high}")
    return float(value)


def _expect(stream, expected):
    line = stream.readline().strip()
    if line != expected:
        raise ValueError(f"Expected {expected!r}, found {line!r}")


@dataclass
class Frame:
    step: int
    ids: object
    fractional: object
    cell: object
    pbc: tuple


def frames(path):
    """Read complete native orthorhombic/restricted-triclinic LAMMPS dump frames."""
    import numpy as np

    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="ascii") as stream:
        while True:
            line = stream.readline()
            if not line:
                return
            if line.strip() != "ITEM: TIMESTEP":
                raise ValueError("Invalid or truncated LAMMPS timestep header")
            step = int(stream.readline())
            _expect(stream, "ITEM: NUMBER OF ATOMS")
            count = int(stream.readline())
            if count <= 0:
                raise ValueError("Empty trajectory frame")
            header = stream.readline().split()
            if header[:3] != ["ITEM:", "BOX", "BOUNDS"]:
                raise ValueError("Missing LAMMPS box bounds")
            tilt = header[3:6] == ["xy", "xz", "yz"]
            flags = header[6:] if tilt else header[3:]
            if len(flags) != 3 or any(
                flag not in ("pp", "ff", "ss", "mm", "fs", "sf", "fm", "mf", "sm", "ms")
                for flag in flags
            ):
                raise ValueError(
                    "Expected three LAMMPS boundary flags on orthorhombic or restricted-triclinic cell"
                )
            pbc = tuple(flag == "pp" for flag in flags)
            bounds = np.asarray([stream.readline().split() for _ in range(3)], float)
            if bounds.shape != (3, 3 if tilt else 2) or not np.isfinite(bounds).all():
                raise ValueError("Invalid box bounds")
            lo = bounds[:, 0].copy()
            hi = bounds[:, 1].copy()
            xy, xz, yz = bounds[:, 2] if tilt else (0.0, 0.0, 0.0)
            lo[0] -= min(0.0, xy, xz, xy + xz)
            hi[0] -= max(0.0, xy, xz, xy + xz)
            lo[1] -= min(0.0, yz)
            hi[1] -= max(0.0, yz)
            lx, ly, lz = hi - lo
            if min(lx, ly, lz) <= 0:
                raise ValueError("Nonpositive lattice height")
            cell = np.array([[lx, 0.0, 0.0], [xy, ly, 0.0], [xz, yz, lz]])
            header = stream.readline().split()
            if header[:2] != ["ITEM:", "ATOMS"] or "id" not in header[2:]:
                raise ValueError("Expected atom IDs and a coordinate triplet")
            cols = header[2:]
            scaled = all(col in cols for col in ("xs", "ys", "zs"))
            xyz = ("xs", "ys", "zs") if scaled else ("x", "y", "z")
            if not all(col in cols for col in xyz):
                raise ValueError("Supported coordinate columns: xs ys zs or x y z")
            ii, *jj = [cols.index(col) for col in ("id", *xyz)]
            ids = np.empty(count, np.int64)
            coords = np.empty((count, 3))
            for i in range(count):
                values = stream.readline().split()
                if len(values) != len(cols):
                    raise ValueError("Truncated or malformed atom row")
                ids[i] = int(values[ii])
                coords[i] = [float(values[j]) for j in jj]
            if (
                not np.isfinite(coords).all()
                or min(ids) <= 0
                or len(np.unique(ids)) != count
            ):
                raise ValueError("Duplicate/invalid atom IDs or nonfinite coordinates")
            frac = coords if scaled else (coords - lo) @ np.linalg.inv(cell)
            for axis, periodic in enumerate(pbc):
                if periodic:
                    frac[:, axis] %= 1.0
            yield Frame(step, ids, frac, cell, pbc)


def _group_indices(groups, ids):
    keys = (
        "reference_plane",
        "comparison_plane",
        "surface_sites",
        "water_oxygen",
        "polymer_heavy",
    )
    indices = {}
    lookup = {int(identifier): i for i, identifier in enumerate(ids)}
    for key in keys:
        values = groups.get(key, [])
        if not isinstance(values, list) or any(
            type(i) is not int or i not in lookup for i in values
        ):
            raise ValueError(f"Invalid or missing atom ID in {key}")
        if len(set(values)) != len(values):
            raise ValueError(f"Duplicate atom in {key}")
        indices[key] = [lookup[i] for i in values]
    if not indices["reference_plane"] or not indices["surface_sites"]:
        raise ValueError(
            "Nonempty reference-plane and surface-site selections required"
        )
    if set(indices["reference_plane"]) & set(indices["comparison_plane"]):
        raise ValueError("Reference and comparison planes must be disjoint")
    solvent, polymer, sites = (
        set(indices[k]) for k in ("water_oxygen", "polymer_heavy", "surface_sites")
    )
    if solvent & polymer or solvent & sites or polymer & sites:
        raise ValueError("Solvent, polymer and surface contact groups must be disjoint")
    return indices


def _plane(z, periodic):
    import numpy as np

    if not periodic:
        return float(np.mean(z))
    # Substrate atom IDs identify one actual layer; unwrap near one member and
    # reject a selection extending across unrelated layers of the periodic cell.
    dz = z - z[0]
    dz -= np.floor(dz + 0.5)
    if np.ptp(dz) >= 0.5:
        raise ValueError("Reference-plane selection spans half the periodic box")
    return float((z[0] + np.mean(dz)) % 1)


def _contact_mask(frame, sites, probes, cutoff):
    import numpy as np
    from scipy.spatial import cKDTree

    if not probes:
        return np.zeros(0, dtype=bool)
    inv = np.linalg.inv(frame.cell)
    # If a Cartesian displacement is <= r, its fractional component i is
    # <= r*||inv(cell)[:,i]||. Enumerate all possible translated neighbours.
    shifts = [
        range(
            -math.ceil(cutoff * np.linalg.norm(inv[:, i])),
            math.ceil(cutoff * np.linalg.norm(inv[:, i])) + 1,
        )
        if on
        else (0,)
        for i, on in enumerate(frame.pbc)
    ]
    if math.prod(map(len, shifts)) > 343:
        raise ValueError("Contact cutoff too large for this cell")
    site_xyz = frame.fractional[sites] @ frame.cell
    offsets = np.array(list(itertools.product(*shifts))) @ frame.cell
    tree = cKDTree((site_xyz[None, :, :] + offsets[:, None, :]).reshape(-1, 3))
    dist, _ = tree.query(frame.fractional[probes] @ frame.cell)
    return dist <= cutoff


def _write_table(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analyze_trajectory(dump, groups, out, *, bin_width_A=0.5, contact_cutoff_A=3.5):
    """Analyse discrete dumped frames. Results are descriptive, without equilibration claims."""
    import numpy as np

    width = _number(bin_width_A, "bin_width_A", 0.05, 50)
    cutoff = _number(contact_cutoff_A, "contact_cutoff_A", 0.1, 20)
    restraint = groups.get("plane_restraint", "unknown")
    if restraint not in ("fixed", "mobile", "unknown"):
        raise ValueError("Invalid plane restraint")
    records = []
    ids_previous = None
    indices = None
    pbc_reference = None
    last_step = None
    for frame in frames(dump):
        if last_step is not None and frame.step <= last_step:
            raise ValueError("Timesteps must increase")
        last_step = frame.step
        if pbc_reference is None:
            pbc_reference = frame.pbc
        elif frame.pbc != pbc_reference:
            raise ValueError("Boundary convention changed during the trajectory")
        if ids_previous is None or not np.array_equal(ids_previous, frame.ids):
            if ids_previous is not None and not np.array_equal(
                np.sort(ids_previous), np.sort(frame.ids)
            ):
                raise ValueError("Frame atom IDs changed")
            indices = _group_indices(groups, frame.ids)
            ids_previous = frame.ids.copy()
        water, polymer, sites = (
            indices[k] for k in ("water_oxygen", "polymer_heavy", "surface_sites")
        )
        lz = float(frame.cell[2, 2])
        if not records:
            bins = max(2, 2 * math.ceil(lz / (2 * width)))
            if bins > 4096:
                raise ValueError("Requested density bin count is excessive")
            density_sum = np.zeros(bins)
            density_sq_sum = np.zeros(bins)
            mean_bin_volume = 0.0
            height_sum = 0.0
        frac_z = frame.fractional[:, 2]
        plane = _plane(frac_z[indices["reference_plane"]], frame.pbc[2])
        # With free z boundaries there is no physical wrap: measure from box floor.
        z_water = ((frac_z[water] - plane + 0.5) % 1) if frame.pbc[2] else frac_z[water]
        if np.any(z_water < 0) or np.any(z_water > 1):
            raise ValueError("Density selection extends beyond nonperiodic cell bounds")
        counts, _ = np.histogram(z_water, bins=bins, range=(0, 1))
        volume = float(np.linalg.det(frame.cell)) / bins
        density = counts / volume
        density_sum += density
        density_sq_sum += density**2
        mean_bin_volume += volume
        height_sum += lz
        in_contact = _contact_mask(frame, sites, water + polymer, cutoff)
        row = {
            "step": frame.step,
            "water_O": len(water),
            "water_O_in_contact": int(in_contact[: len(water)].sum()),
            "polymer_heavy": len(polymer),
            "polymer_heavy_in_contact": int(in_contact[len(water) :].sum()),
        }
        if indices["comparison_plane"]:
            other = _plane(frac_z[indices["comparison_plane"]], frame.pbc[2])
            dz = abs(other - plane)
            shortest = min(dz, 1 - dz) if frame.pbc[2] else dz
            row["plane_separation_A"] = shortest * lz
            if frame.pbc[2]:
                row["complementary_separation_A"] = (1 - shortest) * lz
        records.append(row)
    if not records:
        raise ValueError("Empty trajectory")
    n = len(records)
    mean = density_sum / n
    sd = np.sqrt(np.maximum(density_sq_sum / n - mean**2, 0))
    shift = 0.5 if pbc_reference[2] else 0.0
    profile = [
        {
            "bin": i,
            "fractional_center": (i + 0.5) / bins - shift,
            "mean_normal_position_A": ((i + 0.5) / bins - shift) * height_sum / n,
            "water_O_number_density_A-3": float(mean[i]),
            "frame_sd_A-3": float(sd[i]),
            "mean_bin_volume_A3": mean_bin_volume / n,
        }
        for i in range(bins)
    ]
    out.mkdir(exist_ok=True, parents=True)
    _write_table(out / "density.csv", profile)
    _write_table(out / "contacts.csv", records)
    result = {
        "n_frames": n,
        "first_step": records[0]["step"],
        "last_step": records[-1]["step"],
        "atom_count": len(ids_previous),
        "density_bins": bins,
        "requested_bin_width_A": width,
        "contact_cutoff_A": cutoff,
        "plane_restraint": restraint,
        "density_reference": "periodic reference-plane center"
        if pbc_reference[2]
        else "nonperiodic lower cell boundary",
        "selection_note": groups.get("selection_note", ""),
        "density_unit": "water oxygen atoms per cubic angstrom",
        "plane_distance_observable": "nearest center-to-center distance along basal-plane normal"
        if "plane_separation_A" in records[0]
        else None,
        "evidence": "Finite-trajectory descriptive statistics. Contact denotes a geometric cutoff. Fixed layers encode the input geometry. XRD comparison requires a matched mobile, equilibrated laminate. Frame SD describes temporal spread.",
    }
    for label in ("water_O", "polymer_heavy"):
        denom = sum(row[label] for row in records)
        result["mean_" + label + "_contact_fraction"] = (
            sum(row[label + "_in_contact"] for row in records) / denom
            if denom
            else None
        )
    if "plane_separation_A" in records[0]:
        result["mean_plane_separation_nm"] = (
            float(np.mean([row["plane_separation_A"] for row in records])) / 10
        )
    return result


def run_trajectory(config, out, *, source_dir):
    """Consume only the file pair produced in a resolved local canonical MD job."""
    directory = Path(source_dir).resolve(strict=True)
    if directory == Path(out).resolve():
        raise ValueError("Trajectory analysis requires a separate output directory")
    metadata = directory / "trajectory-groups.json"
    candidates = [
        directory / subdir / name
        for subdir in ("", "engine")
        for name in ("trajectory.dump", "trajectory.dump.gz")
    ]
    dumps = [p for p in candidates if p.is_file()]
    if len(dumps) != 1 or not metadata.is_file():
        raise ValueError(
            "Source MD job lacks the unique trajectory and stable atom selections"
        )
    dump = dumps[0]
    if any(
        not p.resolve().is_relative_to(directory)
        or p.is_symlink()
        or p.parent.is_symlink()
        for p in (dump, metadata)
    ):
        raise ValueError("Source data must be native artifacts of the selected job")
    result = analyze_trajectory(
        dump,
        json.loads(metadata.read_text(encoding="utf-8")),
        Path(out),
        bin_width_A=config.get("bin_width_A", 0.5),
        contact_cutoff_A=config.get("contact_cutoff_A", 3.5),
    )
    for key, path in (("source_dump_sha256", dump), ("selection_sha256", metadata)):
        with path.open("rb") as handle:
            result[key] = hashlib.file_digest(handle, "sha256").hexdigest()
    result["source_job"] = config.get("source_job", config.get("source_job_id"))
    return result
