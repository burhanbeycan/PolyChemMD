"""Rigid placement of explicit molecules in an XY-periodic, Z-confined slab cell.

Coordinates use angstroms. Distinct backends feed the same independent checks.
No atom typing, electrostatic charge assignment or equilibration occurs here.
"""

from __future__ import annotations

import hashlib
from importlib import metadata
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import sysconfig

AVOGADRO = 6.02214076e23
WATER_MOLAR_MASS = 18.01528
READ_TOL_A = 0.002
ENGINE_BORDER_A = 0.15
ENGINE_DISTANCE_BUFFER_A = 0.035


def _number(c, key, default, minimum, maximum, integer=False):
    value = c.get(key, default)
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, float))
        or not math.isfinite(value)
    ):
        raise ValueError(f"{key} must be a finite number")
    if integer and not isinstance(value, int):
        raise ValueError(f"{key} must be an integer")
    if not minimum <= value <= maximum:
        raise ValueError(f"{key} must be between {minimum} and {maximum}")
    return value


def _write_xyz(path, atoms, coords):
    with path.open("w", encoding="utf-8") as f:
        f.write(f"{len(atoms)}\nPolyChemMD geometry; angstrom\n")
        for a, (x, y, z) in zip(atoms, coords, strict=True):
            f.write(f"{a['element']} {x:.9f} {y:.9f} {z:.9f}\n")


def _pos(atoms):
    import numpy as np

    return np.array([[a[t] for t in ("x", "y", "z")] for a in atoms], dtype=float)


def _simple_atom(element, position, role, charge=0):
    return dict(
        element=element,
        x=position[0],
        y=position[1],
        z=position[2],
        role=role,
        formal_charge=charge,
    )


def _species(shape, chains, waters, pairs):
    # Experimental gas-phase water geometry, used only as a geometric template.
    # The selected water force field would set its own equilibrium geometry.
    angle = math.radians(104.52)
    arm = 0.9572
    water = [
        _simple_atom(e, v, "water")
        for e, v in [
            ("O", (0, 0, 0)),
            ("H", (arm, 0, 0)),
            ("H", (arm * math.cos(angle), arm * math.sin(angle), 0)),
        ]
    ]
    specs = [
        ("slab", shape["templates"]["slab"], shape.get("slab_bonds", []), 1),
        ("polymer", shape["templates"]["polymer"], shape["polymer_bonds"], chains),
        ("sodium", [_simple_atom("Na", (0, 0, 0), "sodium", 1)], [], pairs),
        ("chloride", [_simple_atom("Cl", (0, 0, 0), "chloride", -1)], [], pairs),
        ("water", water, [[0, 1, 1], [0, 2, 1]], waters),
    ]
    return [spec for spec in specs if spec[3]]


def _template_order(specs):
    atoms, bonds, groups = [], [], []
    for name, template, links, count in specs:
        if not template:
            raise ValueError(f"Missing {name} template")
        for copy in range(count):
            offset = len(atoms)
            groups.append((name, offset, len(template)))
            for i, atom in enumerate(template):
                item = dict(atom)
                if name != "slab":
                    item["molecule_id"] = f"{name}:{copy + 1}"
                item["template_atom_index"] = i
                item["id"] = len(atoms)
                atoms.append(item)
            for i, j, order in links:
                bonds.append([offset + i, offset + j, order])
    return atoms, bonds, groups


def packmol_path():
    """Resolve the canonical CLI from this Python environment or PATH."""
    location = Path(sysconfig.get_path("scripts")) / (
        "packmol.exe" if sys.platform == "win32" else "packmol"
    )
    return str(location) if location.is_file() else shutil.which("packmol")


def _external(name, out, specs, box, interval, tol, seed):
    import numpy as np

    work = out / "packing_inputs"
    work.mkdir(exist_ok=True)
    # The native engines take a 3D PBC box. Padding z beyond two cells prevents
    # a fictitious z image, while x/y use the real periodic lengths.
    lines = [
        f"tolerance {tol + ENGINE_DISTANCE_BUFFER_A:.7f}",
        "filetype xyz",
        "output engine.xyz",
        f"pbc {box[0]:.9f} {box[1]:.9f} {2 * box[2] + 10:.9f}",
        f"seed {seed}",
        "nloop 250",
    ]
    for species, atoms, _, count in specs:
        _write_xyz(work / f"{species}.xyz", atoms, _pos(atoms))
        lines.extend([f"structure {species}.xyz", f" number {count}"])
        if species == "slab":
            lines.append(" fixed 0 0 0 0 0 0")
        else:
            lines.append(
                f" inside box {ENGINE_BORDER_A} {ENGINE_BORDER_A} {interval[0] + ENGINE_BORDER_A:.9f} "
                f"{box[0] - ENGINE_BORDER_A:.9f} {box[1] - ENGINE_BORDER_A:.9f} {interval[1] - ENGINE_BORDER_A:.9f}"
            )
        lines.append("end structure")
    inp = work / "system.inp"
    inp.write_text("\n".join(lines) + "\n")
    if name == "molpack":
        try:
            import molpack
        except ImportError as err:
            raise RuntimeError(
                "MolPack unavailable. Install molcrafts-molpack==0.4.0 (Python 3.12+)."
            ) from err
        job = molpack.load_script(str(inp.resolve()))
        state = job.packer.with_precision(1e-5).run(job.targets, max_loops=job.nloop)
        if not state.converged:
            raise RuntimeError(
                f"MolPack did not converge (fdist={state.fdist}, frest={state.frest})"
            )
        return (
            list(state.elements),
            np.asarray(state.positions),
            {
                "package": "molcrafts-molpack",
                "version": molpack.version,
                "converged": True,
                "fdist": state.fdist,
                "frest": state.frest,
            },
        )
    binary = packmol_path()
    if not binary:
        raise RuntimeError(
            "Packmol unavailable. Install the canonical Fortran CLI (tested 21.2.0)."
        )
    output = work / "engine.xyz"
    output.unlink(missing_ok=True)
    # 21.x uses stdin; pass no untrusted executable name or shell fragment.
    with inp.open() as input_stream:
        proc = subprocess.run(
            [binary],
            stdin=input_stream,
            text=True,
            cwd=work,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=180,
        )
    (work / "packmol.log").write_text(proc.stdout, encoding="utf-8")
    if proc.returncode or not output.is_file() or "Success!" not in proc.stdout:
        raise RuntimeError(
            "Packmol did not converge; inspect packing_inputs/packmol.log"
        )
    rows = output.read_text().splitlines()
    data = [line.split() for line in rows[2:] if line.strip()]
    if len(data) != int(rows[0]):
        raise ValueError("Packmol XYZ line count mismatch")
    banner = next(
        (line.strip() for line in proc.stdout.splitlines() if "Version " in line),
        "unknown",
    )
    return (
        [r[0] for r in data],
        np.array([[float(v) for v in r[1:4]] for r in data]),
        {
            "package": "packmol",
            "version": banner,
            "executable_sha256": hashlib.sha256(Path(binary).read_bytes()).hexdigest(),
        },
    )


def _periodic_tree(coords, box):
    import numpy as np
    from scipy.spatial import cKDTree

    # Shift z into the centre of a long virtual cell to remove z image contacts.
    supercell = np.array([box[0], box[1], 4 * box[2] + 10])
    offset = np.array([0.0, 0.0, box[2] + 2])
    return cKDTree((coords + offset) % supercell, boxsize=supercell), offset, supercell


def _self_clear(points, box, tol):
    import numpy as np
    from scipy.spatial import cKDTree

    if len(points) < 2 or all(np.ptp(points, axis=0)[:2] < box[:2] - tol):
        return True
    tree = cKDTree(points)
    for ix, iy in ((1, 0), (0, 1), (1, 1), (1, -1)):
        if tree.query(points + [ix * box[0], iy * box[1], 0])[0].min() < tol:
            return False
    return True


def _check_rigid_fit(specs, box, bounds):
    """Rotation-invariant necessary conditions for the unwrapped envelope used by all backends."""
    import numpy as np
    from scipy.spatial.distance import pdist

    sides = np.array([box[0], box[1], bounds[1] - bounds[0]]) + 2 * READ_TOL_A
    for species, atoms, _, _ in specs[1:]:
        if len(atoms) < 2:
            continue
        points = _pos(atoms)
        if np.max(pdist(points)) > np.linalg.norm(sides):
            raise ValueError(
                f"{species}: rigid template diameter exceeds accessible box diagonal"
            )
        # Along any axis spanning width w, the variance is <= w**2/4.
        # The smallest eigenvalue lower-bounds the variance in *every* direction.
        centered = points - points.mean(axis=0)
        min_variance = max(
            0.0, np.linalg.eigvalsh(centered.T @ centered / len(atoms))[0]
        )
        if 2 * math.sqrt(min_variance) > min(sides):
            raise ValueError(
                f"{species}: rigid template minimum width exceeds the narrowest box dimension"
            )


def _slit_pose(template, axes, rng, lo, hi, attempt):
    """Find a feasible translation interval for a rigid orientation proposal in a slit."""
    import numpy as np
    from scipy.spatial.transform import Rotation

    if attempt % 2:
        rotated = Rotation.random(random_state=rng).apply(template)
    else:
        # Prefer an orientation with the largest-spread direction in the basal plane.
        # Rotate around the layer normal and permit either face to point to the surface.
        xy_angle = rng.uniform(-np.pi, np.pi)
        tilted = template @ axes
        if rng.integers(2):
            tilted *= [1.0, -1.0, -1.0]
        rotated = Rotation.from_rotvec([0, 0, xy_angle]).apply(tilted)
    start = lo - rotated.min(axis=0)
    stop = hi - rotated.max(axis=0)
    if np.any(start > stop):
        return None
    return rotated + rng.uniform(start, stop)


class _InsertionExhausted(RuntimeError):
    def __init__(self, species, copy, count, attempts):
        self.attempts = attempts
        super().__init__(
            f"Rigid insertion exhausted at {species} {copy}/{count}. "
            "Use a larger accessible cell, a different packer or a lower initial density."
        )


def _builtin(specs, box, bounds, tol, seed, *, slit=False, matched_constraints=False):
    """Sequential rigid insertion with a randomized solvent grid followed by rejection.

    The grid is an initialization device. Close-to-packed or oversized recipes
    may fail; failed jobs are never reported as converged liquid structures.
    """
    import numpy as np
    from scipy.spatial.transform import Rotation

    rng = np.random.default_rng(seed)
    border = ENGINE_BORDER_A if matched_constraints else 0.001
    lo = np.array([border, border, bounds[0] + border])
    hi = np.array([box[0] - border, box[1] - border, bounds[1] - border])
    points = [_pos(specs[0][1])]
    tree, shift, cell = _periodic_tree(points[0], box)
    pending = np.empty((0, 3))
    target = tol + (ENGINE_DISTANCE_BUFFER_A if matched_constraints else 0.0001)

    # Rebuild a balanced periodic tree in batches, check recent insertions directly.
    def acceptable(candidate):
        if np.any(candidate < lo) or np.any(candidate > hi):
            return False
        if tree.query((candidate + shift) % cell)[0].min() < target:
            return False
        if len(pending):
            delta = candidate[:, None, :] - pending[None, :, :]
            delta[:, :, :2] -= box[:2] * np.rint(delta[:, :, :2] / box[:2])
            if np.any(np.einsum("ijk,ijk->ij", delta, delta) < target**2):
                return False
        return True

    attempts = 0
    polymer_proposals = 0
    for name, atoms, _, count in specs[1:]:
        template = _pos(atoms)
        template -= template.mean(axis=0)
        axes = None
        if slit and name == "polymer":
            _, _, rows = np.linalg.svd(template, full_matrices=False)
            # Columns map from the molecular frame to principal-axis coordinates.
            axes = rows.T.copy() if rows.shape == (3, 3) else np.eye(3)
            axes[:, -1] *= np.sign(np.linalg.det(axes))
        grid = iter(())
        if name == "water":
            # A staggered shape-agnostic candidate grid with randomly ordered
            # sites. Geometric exclusions are applied before accepting each water.
            dims = np.maximum(1, np.floor((hi - lo) / (tol + 0.74)).astype(int))
            sites = (np.indices(dims).reshape(3, -1).T + 0.5) * (hi - lo) / dims + lo
            grid = iter(sites[rng.permutation(len(sites))])
        for i in range(count):
            for trial in range(3500):
                attempts += 1
                if name == "polymer":
                    polymer_proposals += 1
                if axes is not None:
                    candidate = _slit_pose(template, axes, rng, lo, hi, trial)
                    if (
                        candidate is not None
                        and acceptable(candidate)
                        and _self_clear(candidate, box, target)
                    ):
                        break
                    continue
                centre = next(grid, None) if name == "water" else None
                if centre is None:
                    centre = rng.uniform(lo, hi)
                for _orientation in range(4 if name == "water" else 1):
                    candidate = (
                        Rotation.random(random_state=rng).apply(template) + centre
                    )
                    if acceptable(candidate) and _self_clear(candidate, box, target):
                        break
                else:
                    continue
                break
            else:
                raise _InsertionExhausted(name, i + 1, count, attempts)
            points.append(candidate)
            pending = np.concatenate((pending, candidate))
            if len(pending) >= 60:
                tree, shift, cell = _periodic_tree(np.concatenate(points), box)
                pending = np.empty((0, 3))
    return np.concatenate(points), {
        "package": "polychemmd.studio",
        "method": "envelope-conditioned rigid insertion v2"
        if slit
        else "seeded rigid insertion v1",
        "polymer_orientation_proposals": polymer_proposals,
        "orientation_policy": (
            "Alternating in-plane principal-axis and isotropic proposals; feasible centers sampled uniformly."
            if slit
            else "Isotropic orientation proposals; center sampled in the box."
        ),
        "placement_attempts": attempts,
        "numpy": np.__version__,
        "scipy": metadata.version("scipy"),
    }


def _run_builtin(specs, box, bounds, tol, seed, *, slit, matched_constraints):
    import numpy as np

    # Retry confined packings from scratch after a dead end, including new water grid order.
    # Stop at the first audited candidate, without ranking any candidate by contacts or energy.
    limit = 3 if slit else 1
    failures = []
    efforts = 0
    for number in range(limit):
        placement_seed = (
            seed
            if not number
            else int(np.random.SeedSequence([seed, number]).generate_state(1)[0])
        )
        try:
            coords, info = _builtin(
                specs,
                box,
                bounds,
                tol,
                placement_seed,
                slit=slit,
                matched_constraints=matched_constraints,
            )
        except _InsertionExhausted as exc:
            efforts += exc.attempts
            failures.append(
                {
                    "seed": placement_seed,
                    "reason": str(exc),
                    "placement_attempts": exc.attempts,
                }
            )
            if number + 1 == limit:
                if not slit:
                    raise
                raise RuntimeError(
                    f"{limit} confined rigid starts exhausted; final attempt: {exc}"
                ) from exc
        else:
            info.update(
                rejected_starts=failures,
                successful_placement_seed=placement_seed,
                total_placement_attempts=efforts + info["placement_attempts"],
            )
            return coords, info
    raise AssertionError("positive restart budget expected")


def _audit(atoms, coords, groups, specs, box, bounds, tol):
    """Independent minimum-image and rigid-identity checks on accepted coordinates."""
    import numpy as np
    from scipy.spatial.distance import pdist

    if coords.shape != (len(atoms), 3) or not np.isfinite(coords).all():
        raise ValueError("Invalid packed coordinates")
    slab = _pos(specs[0][1])
    slab_size = len(slab)
    if np.max(np.abs(coords[:slab_size] - slab)) > READ_TOL_A:
        raise ValueError("The packing engine changed fixed surface coordinates")
    mobile = coords[slab_size:]
    lo = np.array([0.0, 0.0, bounds[0]])
    hi = np.array([box[0], box[1], bounds[1]])
    if len(mobile) and (
        np.any(mobile < lo - READ_TOL_A) or np.any(mobile > hi + READ_TOL_A)
    ):
        raise ValueError("Packing violates the accessible region")
    refs = {
        name: pdist(_pos(template)) for name, template, _, _ in specs if name != "slab"
    }
    centered = {
        name: _pos(template) - _pos(template).mean(axis=0)
        for name, template, _, _ in specs
        if name != "slab" and len(template) >= 4
    }
    mol = np.empty(len(mobile), dtype=int)
    for idx, (name, start, n) in enumerate(groups[1:]):
        s = coords[start : start + n]
        # All intramolecular distances test template geometry and atom ordering.
        if not np.allclose(pdist(s), refs[name], rtol=0, atol=READ_TOL_A):
            raise ValueError(f"Packing changed the intra-molecular geometry of {name}")
        if name in centered:
            # Distances alone also accept a mirror image of a stereogenic template.
            ref = centered[name]
            observed = s - s.mean(axis=0)
            left, _, right = np.linalg.svd(ref.T @ observed)
            left[:, -1] *= np.sign(np.linalg.det(left @ right))
            residual = np.linalg.norm(ref @ left @ right - observed, axis=1).max()
            if residual > 2 * READ_TOL_A:
                raise ValueError(
                    f"Packing changed the handedness or rigid coordinates of {name}"
                )
        if not _self_clear(s, box, tol - READ_TOL_A):
            raise ValueError(f"{name} overlaps its periodic image")
        mol[start - slab_size : start - slab_size + n] = idx
    closest = None
    if len(mobile):
        tree, offset, periodic = _periodic_tree(slab, box)
        closest = float(tree.query((mobile + offset) % periodic)[0].min())
        # Only pairs capable of lowering the present minimum need checking.
        # Bound query radius for sparse/huge cells, expanding until a different
        # molecule is found, or the surface minimum itself bounds the search.
        tree, _, _ = _periodic_tree(mobile, box)
        radius = min(closest, max(tol + 0.1, 3.5))
        while len(groups) > 2:
            pairs = tree.query_pairs(radius, output_type="ndarray")
            pairs = pairs[mol[pairs[:, 0]] != mol[pairs[:, 1]]]
            if len(pairs):
                delta = mobile[pairs[:, 0]] - mobile[pairs[:, 1]]
                delta[:, :2] -= box[:2] * np.rint(delta[:, :2] / box[:2])
                closest = min(closest, float(np.linalg.norm(delta, axis=1).min()))
                break
            if radius >= closest:
                break
            radius = min(2 * radius, closest)
        if closest < tol - READ_TOL_A:
            raise ValueError(
                f"Minimum XY-image intermolecular distance {closest:.5f} Å below {tol} Å"
            )
    return {
        "minimum_mobile_to_other_A": closest,
        "fixed_surface_invariant": True,
        "rigid_molecule_identity": True,
        "proper_rigid_transforms": True,
        "mobile_self_image_check": True,
        "pbc": [True, True, False],
        "read_tolerance_A": READ_TOL_A,
    }


def pack_system(config: dict, out: Path, *, matched_constraints=False) -> dict:
    """Pack one job directory after geometry.build_geometry; return an audit manifest."""
    import numpy as np

    out = Path(out)
    source = out / "structure.json"
    shape = json.loads(source.read_text(encoding="utf-8"))
    box = np.array(shape["box_A"], dtype=float)
    bounds = np.array(shape["accessible_z_A"], dtype=float)
    if box.shape != (3,) or not np.isfinite(box).all() or np.any(box <= 0):
        raise ValueError("Invalid cell")
    if bounds.shape != (2,) or not 0 <= bounds[0] < bounds[1] <= box[2]:
        raise ValueError("Invalid accessible region")
    volume = float(box[0] * box[1] * (bounds[1] - bounds[0]))
    litre = volume * 1.0e-27
    salt = _number(config, "salt_M", 0.15, 0.0, 5.0)
    density = _number(config, "water_density", 1.0, 0.0, 2.0)
    # Round to nearest whole neutral salt pair; the small-cell realized
    # molarity and the one-pair increment are more meaningful than decimals.
    pairs = math.floor(salt * AVOGADRO * litre + 0.5)
    waters = (
        _number(config, "water_count", 0, 0, 12000, integer=True)
        if config.get("water_count") is not None
        else math.floor(density * volume * 1e-24 / WATER_MOLAR_MASS * AVOGADRO + 0.5)
    )
    chains = _number(
        config, "chains", shape.get("requested_chains", 1), 0, 96, integer=True
    )
    if chains != shape.get("requested_chains", chains):
        raise ValueError("Geometry and packing chain counts disagree")
    if not shape["templates"]["polymer"] and chains:
        raise ValueError("Positive chain count requires a polymer template")
    tol = _number(config, "min_distance_A", 1.8, 0.75, 4.0)
    seed = _number(
        config, "packing_seed", config.get("seed", 42), 1, 2147483646, integer=True
    )
    name = config.get("packing_engine", "builtin")
    if name not in ("builtin", "slit", "molpack", "packmol"):
        raise ValueError("Unknown packing engine")
    specs = _species(shape, chains, waters, pairs)
    natom = sum(len(template) * count for _, template, _, count in specs)
    if natom > 50000 or waters > 12000:
        raise ValueError(
            "Starting cell exceeds 50,000 atoms or 12,000 waters; reduce the cell"
        )
    _check_rigid_fit(specs, box, bounds)
    atoms, bonds, groups = _template_order(specs)
    if len(specs) == 1:
        coords, provenance = (
            _pos(specs[0][1]),
            {"package": "none", "method": "fixed surface only"},
        )
    elif name in ("builtin", "slit"):
        coords, provenance = _run_builtin(
            specs,
            box,
            bounds,
            tol,
            seed,
            slit=name == "slit",
            matched_constraints=matched_constraints,
        )
    else:
        elements, coords, provenance = _external(
            name, out, specs, box, bounds, tol, seed
        )
        if [s.capitalize() for s in elements] != [
            a["element"].capitalize() for a in atoms
        ]:
            raise ValueError("Packing engine changed atom order or element identity")
    audit = _audit(atoms, coords, groups, specs, box, bounds, tol)
    for a, (x, y, z) in zip(atoms, coords, strict=True):
        a.update(x=float(x), y=float(y), z=float(z))
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    packed = {
        "geometry_input_sha256": source_sha,
        "units": "angstrom",
        "box_A": box.tolist(),
        "pbc": [True, True, False],
        "geometry": shape.get("geometry"),
        "accessible_z_A": bounds.tolist(),
        "atoms": atoms,
        "bonds": bonds,
    }
    _write_xyz(out / "packed.xyz", atoms, coords)
    (out / "packed.json").write_text(
        json.dumps(packed, separators=(",", ":"), allow_nan=False)
    )
    # Geometric loading counts every chain in the fluid, independent of contact.
    faces = 2 if shape.get("geometry") == "slit" else 1
    exposed_area_nm2 = float(faces * box[0] * box[1] * 0.01)
    molar_mass = shape.get("polymer", {}).get("molar_mass_g_mol") if chains else 0.0
    nominal_loading = (
        None
        if molar_mass is None
        else chains * molar_mass * 1e21 / (AVOGADRO * exposed_area_nm2)
    )
    report = {
        "stage": "rigid_starting_geometry",
        "engine": name,
        "engine_info": provenance,
        "seed": seed,
        "atom_count": len(atoms),
        "chains": chains,
        "water_count": waters,
        "chain_molarity_M": chains / (AVOGADRO * litre),
        "exposed_faces": faces,
        "total_exposed_area_nm2": exposed_area_nm2,
        "nominal_chains_per_nm2": chains / exposed_area_nm2,
        "nominal_polymer_mass_per_exposed_area_mg_m2": nominal_loading,
        "area_loading_basis": "All chains supplied to the fluid per exposed planar area; adsorption has not been measured",
        "ion_pairs": pairs,
        "salt_M_requested": salt,
        "salt_M_realized": pairs / (AVOGADRO * litre),
        "salt_M_per_pair": 1 / (AVOGADRO * litre),
        "accessible_volume_A3": volume,
        "volume_basis": "Geometric accessible solvent envelope, before mobile solute excluded volume",
        "water_count_basis": "explicit"
        if config.get("water_count") is not None
        else "nominal density",
        "water_density_realized_g_cm3": waters
        * WATER_MOLAR_MASS
        / (AVOGADRO * volume * 1e-24),
        "minimum_distance_target_A": tol,
        "audit": audit,
        "geometry_input_sha256": source_sha,
        "packed_xyz_sha256": hashlib.sha256(
            (out / "packed.xyz").read_bytes()
        ).hexdigest(),
    }
    (out / "packing_report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report
