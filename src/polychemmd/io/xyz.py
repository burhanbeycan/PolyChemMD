"""Extended XYZ coordinate views, with authoritative metadata kept in JSON.

These functions export coordinates for visualization. They do not emit engine
restart files, force fields, bonds, isotope masses or scientific-validation
certificates. Row ordering matches the companion input record exactly.
"""

import math

from ..geometry import MolecularGeometry, _ELEMENTS


def _number(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must contain finite numbers, not booleans")
    try:
        value = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{name} is outside the finite numeric range") from exc
    if not math.isfinite(value):
        raise ValueError(f"{name} must contain only finite numbers")
    return value


def _triple(value, name):
    if not isinstance(value, (list, tuple)) or len(value) != 3:
        raise ValueError(f"{name} requires exactly three components")
    return tuple(_number(item, name) for item in value)


def _rows(elements, coordinates, header):
    lines = [str(len(elements)), header]
    for element, point in zip(elements, coordinates):
        if not isinstance(element, str) or element not in _ELEMENTS:
            raise ValueError("each exported species must be a chemical element symbol")
        point = _triple(point, "coordinates")
        lines.append(element + " " + " ".join(format(value, ".17g") for value in point))
    return "\n".join(lines) + "\n"


def _header(periodic, representation):
    pbc = "T T T" if periodic else "F F F"
    return (f'Properties=species:S:1:pos:R:3 pbc="{pbc}" units=angstrom '
            f'coordinates={representation} visualization_only=T '
            'companion_JSON="topology isotopes atom_identity"')


def molecular_xyz(geometry: MolecularGeometry) -> str:
    """Return a nonperiodic coordinate view in the geometry's atom order."""
    if not isinstance(geometry, MolecularGeometry):
        raise ValueError("geometry must be a MolecularGeometry record")
    return _rows(tuple(atom.element for atom in geometry.atoms),
                 geometry.coordinates_angstrom, _header(False, "unwrapped"))


def cell_xyz(cell_record: dict, *, wrapped: bool = False) -> str:
    """Export the selected cell coordinates without trusting validation claims.

    This validates only the coordinate-view fields it consumes. Unwrapped
    molecular coordinates are the default; wrapped coordinates place atom
    positions inside the cell but may visually split a molecule across faces.
    """
    if not isinstance(wrapped, bool):
        raise ValueError("wrapped must be a boolean")
    if not isinstance(cell_record, dict):
        raise ValueError("cell_record must be a packed-cell mapping")
    if cell_record.get("kind") != "orthorhombic_rigid_packed_cell":
        raise ValueError("unsupported cell kind for XYZ export")
    if cell_record.get("coordinate_units") != "angstrom":
        raise ValueError("XYZ export requires explicit angstrom coordinates")
    atoms, count = cell_record.get("atoms"), cell_record.get("atom_count")
    if isinstance(count, bool) or not isinstance(count, int) or count <= 0:
        raise ValueError("atom_count must be a positive integer")
    if not isinstance(atoms, list) or len(atoms) != count:
        raise ValueError("atom_count must match a nonempty atom array")
    box = cell_record.get("box")
    if not isinstance(box, dict):
        raise ValueError("packed cell requires a box mapping")
    lengths = _triple(box.get("lengths_angstrom"), "box lengths")
    if any(value <= 0 for value in lengths):
        raise ValueError("box lengths must be positive")
    periodic = box.get("periodic")
    if (not isinstance(periodic, (list, tuple)) or len(periodic) != 3
            or any(value is not True for value in periodic)):
        raise ValueError("packed-cell XYZ export requires three periodic axes")
    if "origin_angstrom" in box and any(_triple(box["origin_angstrom"], "box origin")):
        raise ValueError("packed-cell XYZ export supports a zero-origin box only")
    field = "wrapped_coordinates_angstrom" if wrapped else "coordinates_angstrom"
    elements, coordinates = [], []
    for atom in atoms:
        if not isinstance(atom, dict) or "element" not in atom or field not in atom:
            raise ValueError("each atom requires its element and selected coordinate field")
        elements.append(atom["element"])
        point = _triple(atom[field], "coordinates")
        if wrapped and any(not 0 <= value < length for value, length in zip(point, lengths)):
            raise ValueError("wrapped coordinates must lie within [0, box length)")
        coordinates.append(point)
    lattice = (lengths[0], 0., 0., 0., lengths[1], 0., 0., 0., lengths[2])
    header = (_header(True, "wrapped" if wrapped else "unwrapped")
              + ' Lattice="' + " ".join(format(value, ".17g") for value in lattice) + '"')
    return _rows(elements, coordinates, header)
