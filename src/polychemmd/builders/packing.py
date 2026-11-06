"""Bounded rigid-molecule packing in an orthorhombic periodic cell.

This independent constructor places supplied geometries; it does not relax them,
assign a force field, create an equilibrium amorphous phase or predict density.
The all-atom distance cutoff is an explicit geometric exclusion, not a van der
Waals model. Only proper rotations are used, preserving a supplied chirality.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from numbers import Real
import platform
import random
import statistics
from typing import Sequence

from ..analysis.density import density_from_mass_volume
from ..geometry import MolecularGeometry
from .. import __version__


class PackingError(RuntimeError):
    """A bounded construction failed; no partly packed result is returned."""


def _positive(value, name):
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a positive finite real number")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be finite") from exc
    if not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be positive and finite")
    return value


def _integer(value, name, minimum):
    if type(value) is not int or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def _minimum_image_distance(first, second, box):
    """Nearest-image distance for wrapped points in an orthorhombic box."""
    displacement = []
    for x, y, length in zip(first, second, box):
        delta = x - y
        displacement.append(delta - length * math.floor(delta / length + 0.5))
    return math.hypot(*displacement)


def _rotation_matrix(rng):
    # Four independent centred normal draws, normalised on S^3, induce Haar
    # measure on SO(3). q and -q describe the same proper rotation.
    for _ in range(32):
        q = tuple(rng.gauss(0.0, 1.0) for _ in range(4))
        norm = math.hypot(*q)
        if norm > 0:
            break
    else:
        raise PackingError("rotation generator did not produce a nonzero quaternion")
    w, x, y, z = (v / norm for v in q)
    return (
        (1 - 2*(y*y + z*z), 2*(x*y - w*z), 2*(x*z + w*y)),
        (2*(x*y + w*z), 1 - 2*(x*x + z*z), 2*(y*z - w*x)),
        (2*(x*z - w*y), 2*(y*z + w*x), 1 - 2*(x*x + y*y)),
    )


def _rotate(point, matrix):
    try:
        rotated = tuple(math.fsum(a*b for a, b in zip(row, point)) for row in matrix)
    except (OverflowError, ValueError) as exc:
        raise ValueError("rotated coordinates exceed the finite numeric range") from exc
    if any(not math.isfinite(value) for value in rotated):
        raise ValueError("rotated coordinates exceed the finite numeric range")
    return rotated


def _wrap(point, box):
    wrapped, images = [], []
    for value, length in zip(point, box):
        try:
            image = math.floor(value / length)
            position = value - image * length
        except (OverflowError, ValueError) as exc:
            raise ValueError("wrapped coordinates or image flags exceed the finite numeric range") from exc
        # Floating rounding exactly at the upper face must not produce L.
        if position >= length:
            position, image = 0.0, image + 1
        elif position < 0:
            position, image = position + length, image - 1
        if not math.isfinite(position):
            raise ValueError("wrapped coordinates exceed the finite numeric range")
        wrapped.append(position)
        images.append(image)
    return tuple(wrapped), tuple(images)


def _fits_self_image_bound(points, box, cutoff):
    """Sufficient, conservative condition excluding all nonzero self images."""
    return all(max(p[axis] for p in points) - min(p[axis] for p in points)
               <= length - cutoff for axis, length in enumerate(box))


class _DistanceBudget:
    def __init__(self, limit):
        self.limit, self.used = limit, 0

    def distance(self, first, second, box):
        if self.used >= self.limit:
            raise PackingError("max_distance_checks exhausted; no packed result was accepted")
        self.used += 1
        return _minimum_image_distance(first, second, box)


@dataclass(frozen=True)
class PackedCellResult:
    """Immutable serialized construction record; ``to_dict`` returns a copy."""

    _payload_json: str

    def to_dict(self):
        return json.loads(self._payload_json)


def pack_cell(
    geometries: Sequence[MolecularGeometry],
    *,
    box_lengths_angstrom: Sequence[float],
    seed: int,
    min_distance_angstrom: float = 1.5,
    max_attempts_per_molecule: int = 500,
    max_atoms: int = 10000,
    max_distance_checks: int = 5000000,
) -> PackedCellResult:
    """Place supplied whole molecules without changing their internal geometry.

    A local RNG, Gaussian-unit-quaternion rotations and uniformly translated
    molecular centroids make this reproducible for the same inputs and Python
    implementation. Sequential rejection packing is not an equilibrium sampler.
    Every trial must meet the conservative self-image extent bound and all
    accepted intermolecular pairs must be at least ``min_distance_angstrom``.
    Attempt, atom and distance-check budgets fail explicitly, with no result.
    """
    if not isinstance(geometries, Sequence) or isinstance(geometries, (str, bytes)):
        raise ValueError("geometries must be a nonempty sequence of MolecularGeometry")
    geometries = tuple(geometries)
    if not geometries or any(not isinstance(g, MolecularGeometry) for g in geometries):
        raise ValueError("geometries must be nonempty and contain only MolecularGeometry")
    if (not isinstance(box_lengths_angstrom, Sequence)
            or isinstance(box_lengths_angstrom, (str, bytes)) or len(box_lengths_angstrom) != 3):
        raise ValueError("box_lengths_angstrom requires exactly three lengths")
    try:
        box = tuple(_positive(v, "box length") for v in box_lengths_angstrom)
    except TypeError as exc:
        raise ValueError("box_lengths_angstrom requires exactly three lengths") from exc
    if len(box) != 3:
        raise ValueError("box_lengths_angstrom requires exactly three lengths")
    _integer(seed, "seed", 0)
    _integer(max_attempts_per_molecule, "max_attempts_per_molecule", 1)
    _integer(max_atoms, "max_atoms", 1)
    _integer(max_distance_checks, "max_distance_checks", 1)
    cutoff = _positive(min_distance_angstrom, "min_distance_angstrom")
    if any(length < cutoff for length in box):
        raise PackingError("cell is smaller than the cutoff to a molecule's own periodic image")
    volume = math.prod(box)
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError("box volume is outside the positive finite numeric range")
    atom_count = sum(len(g.atoms) for g in geometries)
    if atom_count > max_atoms:
        raise ValueError("total atom count exceeds max_atoms")
    masses = tuple(atom.mass_g_mol for g in geometries for atom in g.atoms)
    try:
        total_mass = math.fsum(masses)
    except OverflowError as exc:
        raise ValueError("total molar mass is outside the finite numeric range") from exc
    if not math.isfinite(total_mass):
        raise ValueError("total molar mass is outside the finite numeric range")
    constructed_density = density_from_mass_volume(masses, volume)
    rng, budget = random.Random(seed), _DistanceBudget(max_distance_checks)
    placed, molecule_records = [], []

    for molecule_index, geometry in enumerate(geometries):
        positions = geometry.coordinates_angstrom
        if len(set(positions)) != len(positions):
            raise ValueError("a supplied molecular geometry contains coincident atoms")
        try:
            center = tuple(float(statistics.mean(p[axis] for p in positions)) for axis in range(3))
            centered = tuple(tuple(v-c for v, c in zip(p, center)) for p in positions)
        except (OverflowError, ValueError) as exc:
            raise ValueError("molecular centering exceeds the finite numeric range") from exc
        if any(not math.isfinite(v) for p in centered for v in p):
            raise ValueError("molecular extent is outside the finite numeric range")
        accepted = None
        self_bound_rejections = 0
        contact_rejections = 0
        for attempt in range(1, max_attempts_per_molecule + 1):
            matrix = _rotation_matrix(rng)
            rotated = tuple(_rotate(p, matrix) for p in centered)
            if not _fits_self_image_bound(rotated, box, cutoff):
                self_bound_rejections += 1
                continue
            translation = tuple(rng.random() * length for length in box)
            unwrapped = tuple(tuple(a+b for a, b in zip(p, translation)) for p in rotated)
            if any(not math.isfinite(v) for p in unwrapped for v in p):
                raise ValueError("transformed coordinates exceed the finite numeric range")
            wrapped_and_images = tuple(_wrap(p, box) for p in unwrapped)
            wrapped = tuple(pair[0] for pair in wrapped_and_images)
            images = tuple(pair[1] for pair in wrapped_and_images)
            collides = False
            for earlier in placed:
                for point in wrapped:
                    if any(budget.distance(point, old, box) < cutoff for old in earlier[1]):
                        collides = True
                        break
                if collides:
                    break
            if collides:
                contact_rejections += 1
                continue
            accepted = (unwrapped, wrapped, images)
            break
        if accepted is None:
            raise PackingError(
                f"molecule {molecule_index} exhausted {max_attempts_per_molecule} attempts "
                f"(self-image extent rejections={self_bound_rejections}, contact rejections={contact_rejections}); "
                "no partial packing returned"
            )
        placed.append(accepted)
        molecule_records.append({
            "molecule_index": molecule_index,
            "source_geometry": geometry.to_dict(),
            "source_graph_sha256": geometry.source_graph_sha256,
            "placement": {
                "rotation_matrix": [list(row) for row in matrix],
                "source_centroid_angstrom": list(center),
                "translated_centroid_angstrom": list(translation),
                "attempts": attempt,
                "self_image_extent_rejections": self_bound_rejections,
                "intermolecular_contact_rejections": contact_rejections,
            },
        })

    # An independent full-pair audit rechecks the actual accepted coordinates.
    minimum = None
    for i, current in enumerate(placed):
        for earlier in placed[:i]:
            for point in current[1]:
                for other in earlier[1]:
                    distance = budget.distance(point, other, box)
                    if distance < cutoff:
                        raise PackingError("internal final audit found an intermolecular contact")
                    minimum = distance if minimum is None else min(minimum, distance)

    atoms, bonds = [], []
    maximum_bond_relative_error = 0.0
    for molecule_index, (geometry, coordinates) in enumerate(zip(geometries, placed)):
        unwrapped, wrapped, images = coordinates
        indices = {atom.atom_id: i for i, atom in enumerate(geometry.atoms)}
        global_maps = {}
        for i, atom in enumerate(geometry.atoms):
            global_map = len(atoms) + 1
            global_maps[atom.atom_id] = global_map
            atoms.append({
                "cell_atom_id": f"molecule:{molecule_index}/{atom.atom_id}",
                "atom_map": global_map,
                "source_atom_id": atom.atom_id, "source_atom_map": atom.atom_map,
                "molecule_index": molecule_index, "element": atom.element,
                "isotope": atom.isotope, "formal_charge": atom.formal_charge,
                "mass_g_mol": atom.mass_g_mol, "lineage": json.loads(atom.lineage_json),
                "coordinates_angstrom": list(unwrapped[i]),
                "wrapped_coordinates_angstrom": list(wrapped[i]),
                "image_flags": list(images[i]),
            })
        for bond in geometry.bonds:
            first, second = (indices[name] for name in bond.atom_ids)
            before = math.dist(geometry.coordinates_angstrom[first], geometry.coordinates_angstrom[second])
            after = math.dist(unwrapped[first], unwrapped[second])
            if before <= 0 or not math.isfinite(before) or not math.isfinite(after):
                raise ValueError("supplied bond geometry is degenerate or unrepresentable")
            relative_error = abs(after-before) / before
            maximum_bond_relative_error = max(maximum_bond_relative_error, relative_error)
            if relative_error > 1e-10:
                raise PackingError("rigid placement failed bond-length preservation")
            bonds.append({
                "molecule_index": molecule_index,
                "cell_atom_ids": [f"molecule:{molecule_index}/{name}" for name in bond.atom_ids],
                "atom_maps": [global_maps[name] for name in bond.atom_ids],
                "order": bond.order, "aromatic": bond.aromatic, "stereo": bond.stereo,
            })

    payload = {
        "schema_version": "1", "kind": "orthorhombic_rigid_packed_cell",
        "evidence_level": "geometry_constructed", "simulation_performed": False,
        "coordinate_units": "angstrom",
        "coordinate_representation": "unwrapped molecule coordinates plus wrapped coordinates and image flags",
        "box": {"origin_angstrom": [0.0, 0.0, 0.0], "lengths_angstrom": list(box),
                "periodic": [True, True, True], "volume_angstrom3": volume},
        "molecule_count": len(geometries), "atom_count": atom_count,
        "bond_count": len(bonds), "atoms": atoms, "bonds": bonds,
        "molecules": molecule_records,
        "formal_charge": sum(atom.formal_charge for g in geometries for atom in g.atoms),
        "total_molar_mass_g_mol": total_mass,
        "construction_density_g_cm3": constructed_density,
        "density_interpretation": "mass divided by user-specified cell volume; constructed, not predicted",
        "packing": {
            "method": "sequential_rigid_rejection_orthorhombic_periodic",
            "seed": seed, "random_generator": "Python random.Random",
            "rotation_sampling": "normalised four-Gaussian quaternion; proper uniform SO(3) rotations",
            "translation_sampling": "molecular centroid uniform in each box axis",
            "min_distance_angstrom": cutoff, "cutoff_model": "all-atom geometric exclusion, not van der Waals radii",
            "max_attempts_per_molecule": max_attempts_per_molecule,
            "max_atoms": max_atoms, "max_distance_checks": max_distance_checks,
            "distance_checks_used": budget.used,
            "minimum_inter_molecular_distance_angstrom": minimum,
            "self_image_policy": "every rotated molecular axis extent <= box length minus cutoff",
        },
        "validation": {"intermolecular_minimum_image_audit": True,
                       "self_periodic_image_exclusion_bound": True,
                       "proper_rotations_only": True,
                       "maximum_bond_length_relative_error": maximum_bond_relative_error,
                       "source_identity_and_lineage_preserved": True,
                       "force_field_assigned": False},
        "equilibrium_status": "not_assessed", "phase": "not_assigned",
        "material_validation": "not_performed",
        "software": {"polychemmd": __version__, "python": platform.python_version(),
                     "python_implementation": platform.python_implementation()},
    }
    return PackedCellResult(json.dumps(payload, sort_keys=True, separators=(",", ":"), allow_nan=False))
