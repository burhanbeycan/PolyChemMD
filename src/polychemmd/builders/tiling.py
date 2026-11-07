"""Deterministic molecular-motif tiling with bounded periodic contact checks.

Tiling creates an ordered geometric starting configuration. It does not assign
a crystalline or liquid phase, relax coordinates, or generate independent
material samples. The contact floor is a geometric exclusion, not a force field.
"""

from itertools import combinations, product
import json
import math
import platform

from .. import __version__
from ..analysis.density import density_from_mass_volume
from ..geometry import MolecularGeometry
from .packing import (PackedCellResult, PackingError, _fits_self_image_bound,
                      _minimum_image_distance, _positive, _rotate, _wrap)


_IDENTITY = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))
_HARD_ATOM_LIMIT = 8192


def _integer(value, name, *, maximum=None):
    if type(value) is not int or value < 1 or maximum is not None and value > maximum:
        raise ValueError(f"{name} must be a positive integer" +
                         (f" <= {maximum}" if maximum is not None else ""))
    return value


def _triple(value, name, converter):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} must contain exactly three components")
    return tuple(converter(item, name) for item in value)


def _matrix(value):
    if value is None:
        return _IDENTITY
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError("rotation_matrix must contain three finite rows")
    rows = []
    for row in value:
        if not isinstance(row, (list, tuple)) or len(row) != 3:
            raise ValueError("rotation_matrix must contain three finite rows")
        result = []
        for v in row:
            if isinstance(v, bool) or not isinstance(v, (int, float)):
                raise ValueError("rotation_matrix must contain finite numeric values")
            try:
                v = float(v)
            except (ValueError, OverflowError) as exc:
                raise ValueError("rotation_matrix must be representable") from exc
            if not math.isfinite(v) or abs(v) > 1.0 + 1e-10:
                raise ValueError("rotation_matrix must be finite and orthonormal")
            result.append(v)
        rows.append(tuple(result))
    for i, row in enumerate(rows):
        for j, other in enumerate(rows):
            if not math.isclose(math.fsum(a*b for a, b in zip(row, other)),
                                float(i == j), rel_tol=0, abs_tol=1e-10):
                raise ValueError("rotation_matrix must be orthonormal")
    a, b, c = rows
    determinant = (a[0]*(b[1]*c[2]-b[2]*c[1]) -
                   a[1]*(b[0]*c[2]-b[2]*c[0]) +
                   a[2]*(b[0]*c[1]-b[1]*c[0]))
    if not math.isclose(determinant, 1., rel_tol=0, abs_tol=1e-10):
        raise ValueError("rotation_matrix must be a proper rotation with determinant +1")
    return tuple(rows)


class _Budget:
    def __init__(self, limit):
        self.limit, self.used = limit, 0

    def take(self):
        if self.used >= self.limit:
            raise PackingError("max_distance_checks exhausted; no tiled result accepted")
        self.used += 1


def _contact_audit(molecules, box, cutoff, budget):
    """Check all potentially forbidden pairs using periodic linked bins.

    Bin widths are at least cutoff. Thus a pair closer than cutoff must lie in
    the same or neighboring bins. Wrapped neighboring keys are deduplicated
    when an axis has one or two bins. No global nearest-distance claim follows.
    """
    count = sum(len(points) for points in molecules)
    # Compare after division by the bounded atom count, avoiding an overflowing
    # L/cutoff ratio for extremely small, finite positive contact floors.
    counts = []
    for length in box:
        n = count if length/count >= cutoff else max(1, math.floor(length/cutoff))
        # A rounded-up quotient at an integer boundary must not create bins
        # narrower than the geometric contact floor.
        if n > 1 and length/n < cutoff:
            n -= 1
        counts.append(n)
    bins_per_axis = tuple(counts)
    bins = {}
    minimum = None
    checks = 0
    for molecule_index, points in enumerate(molecules):
        for point in points:
            key = tuple(min(n-1, math.floor(value/length*n))
                        for value, length, n in zip(point, box, bins_per_axis))
            neighbors = {tuple((v+d) % n for v, d, n in zip(key, delta, bins_per_axis))
                         for delta in product((-1, 0, 1), repeat=3)}
            for neighbor in sorted(neighbors):
                for previous_molecule, previous in bins.get(neighbor, ()):
                    if previous_molecule == molecule_index:
                        continue
                    budget.take(); checks += 1
                    distance = _minimum_image_distance(point, previous, box)
                    if not math.isfinite(distance):
                        raise ValueError("periodic contact distance is unrepresentable")
                    minimum = distance if minimum is None else min(minimum, distance)
                    if distance < cutoff:
                        raise PackingError("tiled geometry violates the periodic intermolecular contact floor")
            bins.setdefault(key, []).append((molecule_index, point))
    return bins_per_axis, checks, minimum


def tile_molecular_geometry(
    geometry, *, repeats=(2, 6, 6), spacing_angstrom=(13., 5.1, 5.1),
    rotation_matrix=None, min_distance_angstrom=1.5, max_atoms=8192,
    max_distance_checks=5000000,
):
    """Place one rigid molecule at every center of an orthorhombic motif grid.

    The axis spacings and proper rotation are explicit construction inputs.
    At most 8192 atoms are supported in this bounded builder. Its work budget
    includes every actual intra-molecular rigidity and contact distance check.
    Failed construction returns no partial result and writes no files.
    """
    if not isinstance(geometry, MolecularGeometry):
        raise ValueError("geometry must be a MolecularGeometry")
    # Validate the structural source again, independently of provenance flags.
    geometry = MolecularGeometry.from_mapping(geometry.to_dict())
    repeats = _triple(repeats, "repeats", _integer)
    spacing = _triple(spacing_angstrom, "spacing_angstrom", _positive)
    _integer(max_atoms, "max_atoms", maximum=_HARD_ATOM_LIMIT)
    _integer(max_distance_checks, "max_distance_checks")
    atom_count = len(geometry.atoms)*math.prod(repeats)
    if atom_count > max_atoms:
        raise ValueError("total tiled atom count exceeds max_atoms")
    if math.prod(repeats)*len(geometry.atoms)*(len(geometry.atoms)-1)//2 > max_distance_checks:
        raise PackingError("max_distance_checks cannot cover mandatory intra-molecular rigidity checks")
    matrix = _matrix(rotation_matrix)
    cutoff = _positive(min_distance_angstrom, "min_distance_angstrom")
    try:
        box = tuple(n*s for n, s in zip(repeats, spacing))
        volume = math.prod(box)
    except OverflowError as exc:
        raise ValueError("tiled box is not representable") from exc
    if any(not math.isfinite(v) or v <= 0 for v in box) or not math.isfinite(volume) or volume <= 0:
        raise ValueError("tiled box and volume must be positive and representable")
    if any(length < cutoff for length in box):
        raise PackingError("tiled box is smaller than the geometric self-image cutoff")
    positions = geometry.coordinates_angstrom
    if len(set(positions)) != len(positions):
        raise ValueError("source geometry contains coincident atoms")
    try:
        center = tuple(math.fsum(point[axis]/len(positions) for point in positions) for axis in range(3))
        centered = tuple(tuple(v-c for v, c in zip(point, center)) for point in positions)
    except (OverflowError, ValueError) as exc:
        raise ValueError("source molecular centering is unrepresentable") from exc
    if any(not math.isfinite(v) for point in centered for v in point):
        raise ValueError("source molecular extent is unrepresentable")
    rotated = tuple(_rotate(point, matrix) for point in centered)
    if not _fits_self_image_bound(rotated, box, cutoff):
        raise PackingError("source extent violates the conservative nonzero self-image bound")
    budget = _Budget(max_distance_checks)
    molecule_records, placed = [], []
    maximum_pair_error = 0.
    for index, grid in enumerate(product(*(range(n) for n in repeats))):
        translation = tuple((v+.5)*s for v, s in zip(grid, spacing))
        unwrapped = tuple(tuple(v+t for v, t in zip(point, translation)) for point in rotated)
        if any(not math.isfinite(v) for point in unwrapped for v in point):
            raise ValueError("tiled coordinates are unrepresentable")
        wrapped_images = tuple(_wrap(point, box) for point in unwrapped)
        wrapped = tuple(item[0] for item in wrapped_images)
        images = tuple(item[1] for item in wrapped_images)
        for first, second in combinations(range(len(positions)), 2):
            budget.take()
            before, after = math.dist(positions[first], positions[second]), math.dist(unwrapped[first], unwrapped[second])
            if not math.isfinite(before) or before <= 0 or not math.isfinite(after):
                raise ValueError("source/tiled intra-molecular distance is degenerate or unrepresentable")
            error = abs(after-before)
            maximum_pair_error = max(maximum_pair_error, error)
            if error > 1e-8:
                raise PackingError("tiled coordinates lost source molecular pair distances")
        placed.append((unwrapped, wrapped, images))
        molecule_records.append({
            "molecule_index": index, "source_geometry": geometry.to_dict(),
            "source_graph_sha256": geometry.source_graph_sha256,
            "placement": {"rotation_matrix": [list(row) for row in matrix],
                          "source_centroid_angstrom": list(center),
                          "translated_centroid_angstrom": list(translation),
                          "motif_grid_index": list(grid), "attempts": 1,
                          "self_image_extent_rejections": 0, "intermolecular_contact_rejections": 0}})
    bins, contact_checks, minimum = _contact_audit([row[1] for row in placed], box, cutoff, budget)
    atoms, bonds = [], []
    for molecule_index, (unwrapped, wrapped, images) in enumerate(placed):
        maps = {}
        for index, atom in enumerate(geometry.atoms):
            mapping = len(atoms)+1; maps[atom.atom_id] = mapping
            atoms.append({"cell_atom_id": f"molecule:{molecule_index}/{atom.atom_id}",
                          "atom_map": mapping, "source_atom_id": atom.atom_id,
                          "source_atom_map": atom.atom_map, "molecule_index": molecule_index,
                          "element": atom.element, "isotope": atom.isotope, "formal_charge": atom.formal_charge,
                          "mass_g_mol": atom.mass_g_mol, "lineage": atom.lineage,
                          "coordinates_angstrom": list(unwrapped[index]),
                          "wrapped_coordinates_angstrom": list(wrapped[index]), "image_flags": list(images[index])})
        for bond in geometry.bonds:
            bonds.append({"molecule_index": molecule_index,
                          "cell_atom_ids": [f"molecule:{molecule_index}/{name}" for name in bond.atom_ids],
                          "atom_maps": [maps[name] for name in bond.atom_ids],
                          "order": bond.order, "aromatic": bond.aromatic, "stereo": bond.stereo})
    masses = [atom["mass_g_mol"] for atom in atoms]
    try:
        mass = math.fsum(masses)
    except OverflowError as exc:
        raise ValueError("total tiled molar mass is unrepresentable") from exc
    if not math.isfinite(mass):
        raise ValueError("total tiled molar mass is unrepresentable")
    density = density_from_mass_volume(masses, volume)
    payload = {
        "schema_version": "1", "kind": "orthorhombic_rigid_packed_cell",
        "evidence_level": "geometry_constructed", "simulation_performed": False,
        "coordinate_units": "angstrom",
        "coordinate_representation": "unwrapped molecule coordinates plus wrapped coordinates and image flags",
        "box": {"origin_angstrom": [0., 0., 0.], "lengths_angstrom": list(box),
                "periodic": [True, True, True], "volume_angstrom3": volume},
        "molecule_count": len(placed), "atom_count": atom_count, "bond_count": len(bonds),
        "atoms": atoms, "bonds": bonds, "molecules": molecule_records,
        "formal_charge": sum(atom["formal_charge"] for atom in atoms),
        "total_molar_mass_g_mol": mass, "construction_density_g_cm3": density,
        "density_interpretation": "mass divided by requested motif-grid volume; constructed, not predicted",
        "packing": {"method": "deterministic_rigid_molecular_motif_tiling",
                    "seed": None, "random_generator": None, "stochastic_sampling_performed": False,
                    "repeats": list(repeats), "spacing_angstrom": list(spacing),
                    "rotation_matrix": [list(row) for row in matrix],
                    "min_distance_angstrom": cutoff, "cutoff_model": "geometric exclusion, not a force field",
                    "max_atoms": max_atoms, "max_distance_checks": max_distance_checks,
                    "distance_checks_used": budget.used, "contact_distance_checks_used": contact_checks,
                    "contact_audit": "periodic linked bins; all potentially sub-cutoff intermolecular pairs inspected",
                    "linked_bin_counts": list(bins),
                    "contact_distance_lower_bound_angstrom": cutoff,
                    "minimum_checked_candidate_distance_angstrom": minimum,
                    "minimum_inter_molecular_distance_angstrom": None,
                    "global_minimum_distance_status": "not_measured; candidate minimum is not a global minimum",
                    "self_image_policy": "rotated molecule axis extent <= global box length minus cutoff"},
        "validation": {"intermolecular_minimum_image_audit": True,
                       "self_periodic_image_exclusion_bound": True, "proper_rotations_only": True,
                       "maximum_intramolecular_pair_absolute_error_angstrom": maximum_pair_error,
                       "source_identity_and_lineage_preserved": True, "force_field_assigned": False},
        "construction_interpretation": "ordered repeated molecular motif; not an identified crystal, liquid or amorphous phase",
        "independent_material_cell_count": 1,
        "repeated_molecules_are_independent_material_cells": False,
        "equilibrium_status": "not_assessed", "phase": "not_assigned", "material_validation": "not_performed",
        "software": {"polychemmd": __version__, "python": platform.python_version(),
                     "python_implementation": platform.python_implementation()},
    }
    return PackedCellResult(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False))
