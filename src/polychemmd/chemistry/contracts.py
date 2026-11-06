"""Version 1 contracts for already transformed, two-ended repeat units.

This module does not import RDKit until a chemical object is constructed. A
repeat unit is a connectivity specification, not a free monomer or a proposed
polymerization reaction. Atom indices refer to RDKit's parsed input order.
"""

from collections.abc import Mapping
from dataclasses import asdict, dataclass


class ChemistryDependencyError(RuntimeError):
    """A requested chemistry operation needs the optional RDKit dependency."""


class UnsupportedChemistryError(ValueError):
    """The input is outside the explicitly implemented chemistry boundary."""


def _rdkit():
    try:
        from rdkit import Chem
    except ImportError as exc:
        raise ChemistryDependencyError(
            "Chemical graph operations require the optional RDKit dependency; "
            "install PolyChemMD's chemistry extra or install rdkit."
        ) from exc
    return Chem


def _positive_int(value, name):
    if type(value) is not int or value < 1:
        raise ValueError(f"{name} must be a positive integer")
    return value


def _text(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{name} must be a non-empty string")
    return value


def _mapping_keys(data, required, name):
    if not isinstance(data, Mapping):
        raise ValueError(f"{name} must be an object")
    unknown, missing = set(data) - required, required - set(data)
    if unknown or missing:
        raise ValueError(
            f"{name}: unknown keys={sorted(unknown, key=str)}; "
            f"missing keys={sorted(missing)}"
        )


@dataclass(frozen=True)
class RepeatUnitSpec:
    """A single connected repeat graph with mapped head/tail dummy atoms.

    Both attachments must be neutral, terminal dummies with single bonds. Maps
    on ordinary atoms are optional; every supplied nonzero map must be unique.
    The immutable object retains the exact input SMILES and attachment roles.
    Construction performs graph/valence validation, without making a material,
    geometry, reaction, protonation, or polymerization-mechanism claim.
    """

    label: str
    smiles: str
    head_map: int
    tail_map: int

    def __post_init__(self):
        _text(self.label, "label")
        _text(self.smiles, "smiles")
        _positive_int(self.head_map, "head_map")
        _positive_int(self.tail_map, "tail_map")
        if self.head_map == self.tail_map:
            raise ValueError("head_map and tail_map must be distinct")
        self.to_mol()

    @classmethod
    def from_mapping(cls, data):
        fields = {"label", "smiles", "head_map", "tail_map"}
        _mapping_keys(data, fields, "repeat unit")
        return cls(**data)

    def to_dict(self):
        """Return the lossless input contract, without inserting material state."""
        return asdict(self)

    def to_mol(self):
        """Parse and validate, returning a fresh RDKit Mol with attachment dummies.

        Explicit H atoms in the source SMILES are retained, so source atom
        indices and isotope/atom-map lineage are not silently discarded.
        """
        Chem = _rdkit()
        parameters = Chem.SmilesParserParams()
        parameters.removeHs = False
        parameters.sanitize = False
        parameters.parseName = False
        parameters.allowCXSMILES = False
        mol = Chem.MolFromSmiles(self.smiles, parameters)
        if mol is None:
            raise ValueError("repeat-unit SMILES is invalid or fails valence sanitization")
        requested_chirality = {a.GetIdx(): a.GetChiralTag() for a in mol.GetAtoms()}
        try:
            Chem.SanitizeMol(mol)
            Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        except (ValueError, RuntimeError) as exc:
            raise ValueError("repeat-unit SMILES is invalid or fails valence sanitization") from exc
        if len(Chem.GetMolFrags(mol)) != 1:
            raise ValueError("repeat unit must be a single connected graph; fragments are unsupported")
        dummies = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        if len(dummies) != 2:
            raise ValueError("repeat unit requires exactly two attachment dummy atoms; free monomers are not repeat units")
        nonzero_maps = [a.GetAtomMapNum() for a in mol.GetAtoms() if a.GetAtomMapNum()]
        if len(nonzero_maps) != len(set(nonzero_maps)):
            raise ValueError("all nonzero atom maps must be unique within the repeat unit")
        if {a.GetAtomMapNum() for a in dummies} != {self.head_map, self.tail_map}:
            raise ValueError("the two attachment dummies must carry the specified head_map and tail_map")
        if mol.GetNumAtoms() == 2:
            raise ValueError("repeat unit must contain a non-dummy atom")
        for atom in mol.GetAtoms():
            if atom.GetNumRadicalElectrons():
                raise UnsupportedChemistryError("radical repeat units are not supported")
        dummy_indices = {a.GetIdx() for a in dummies}
        for dummy in dummies:
            if (dummy.GetDegree() != 1 or dummy.GetFormalCharge() != 0
                    or dummy.GetIsotope() != 0 or dummy.GetTotalNumHs() != 0):
                raise ValueError("attachment dummies must be neutral terminal atoms without isotope or hydrogen labels")
            neighbor = dummy.GetNeighbors()[0]
            if neighbor.GetAtomicNum() == 0:
                raise ValueError("each attachment must connect directly to a non-dummy atom")
            bond = mol.GetBondBetweenAtoms(dummy.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() != Chem.BondType.SINGLE:
                raise UnsupportedChemistryError("only single-bond attachments are supported")
            if neighbor.GetAtomicNum() == 1:
                raise UnsupportedChemistryError("a hydrogen cannot serve as a repeat-unit attachment anchor")
            if requested_chirality[neighbor.GetIdx()] != Chem.ChiralType.CHI_UNSPECIFIED:
                raise UnsupportedChemistryError(
                    "stereochemistry at an attachment anchor is unsupported; "
                    "explicit transformation-aware stereochemistry is required"
                )
        for bond in mol.GetBonds():
            if dummy_indices.intersection(bond.GetStereoAtoms()):
                raise UnsupportedChemistryError(
                    "double-bond stereochemistry defined by an attachment dummy is unsupported"
                )
        # Only tetrahedral point chirality and ordinary cis/trans bond stereo
        # have a preservation contract in this first builder.
        for atom in mol.GetAtoms():
            if requested_chirality[atom.GetIdx()] not in (
                Chem.ChiralType.CHI_UNSPECIFIED,
                Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
            ):
                raise UnsupportedChemistryError("non-tetrahedral atom stereochemistry is unsupported")
            if (requested_chirality[atom.GetIdx()] != Chem.ChiralType.CHI_UNSPECIFIED
                    and atom.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED):
                raise ValueError("specified atom stereochemistry is not valid for the repeat graph")
        return mol

    def attachment_indices(self, mol=None):
        """Return dummy and anchor indices, using the unmodified source graph."""
        if mol is None:
            mol = self.to_mol()
        mapped = {a.GetAtomMapNum(): a for a in mol.GetAtoms() if a.GetAtomicNum() == 0}
        head, tail = mapped[self.head_map], mapped[self.tail_map]
        return {
            "head_dummy": head.GetIdx(), "tail_dummy": tail.GetIdx(),
            "head_anchor": head.GetNeighbors()[0].GetIdx(),
            "tail_anchor": tail.GetNeighbors()[0].GetIdx(),
        }

    def structural_roles(self, mol=None):
        """Identify a unique acyclic attachment path, without guessing ring paths.

        Sidegroups may contain rings. If a head-to-tail path traverses a ring,
        the possible backbone routes are ambiguous: all ordinary atoms receive
        ``unassigned`` and the ambiguity is recorded explicitly.
        """
        Chem = _rdkit()
        if mol is None:
            mol = self.to_mol()
        indices = self.attachment_indices(mol)
        head, tail = indices["head_anchor"], indices["tail_anchor"]
        path = (head,) if head == tail else Chem.GetShortestPath(mol, head, tail)
        ambiguous = any(
            mol.GetBondBetweenAtoms(a, b).IsInRing()
            for a, b in zip(path, path[1:])
        )
        roles = {}
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 0:
                role = "attachment"
            elif ambiguous:
                role = "unassigned"
            else:
                role = "backbone" if atom.GetIdx() in path else "sidechain"
            roles[atom.GetIdx()] = role
        # Explicit source H and hydrogens completed later must give the same
        # descriptor grouping for equivalent hydrogen notations.
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1:
                roles[atom.GetIdx()] = roles[atom.GetNeighbors()[0].GetIdx()]
        return {
            "method": "ambiguous_cyclic_path" if ambiguous else "unique_acyclic_attachment_path",
            "backbone_source_atom_indices": [] if ambiguous else list(path),
            "atom_roles": roles,
        }


@dataclass(frozen=True)
class AtomLineage:
    """A graph atom's stable origin, separate from its current index and map.

    Source indices and repeat indices are zero-based. Output atom maps are
    positive and unique only within one assembled molecule; atom_id encodes the
    copy/source relationship and does not depend on SMILES serialization order.
    """

    atom_id: str
    atom_index: int
    atom_map: int
    origin: str
    repeat_index: int | None
    source_atom_index: int | None
    source_atom_map: int | None
    parent_atom_id: str | None
    terminus: str | None
    structural_role: str

    def __post_init__(self):
        _text(self.atom_id, "atom_id")
        if type(self.atom_index) is not int or self.atom_index < 0:
            raise ValueError("atom_index must be a nonnegative integer")
        _positive_int(self.atom_map, "atom_map")
        if self.origin not in {"repeat", "cap", "added_hydrogen"}:
            raise ValueError("unsupported atom lineage origin")
        for name in ("repeat_index", "source_atom_index"):
            value = getattr(self, name)
            if value is not None and (type(value) is not int or value < 0):
                raise ValueError(f"{name} must be a nonnegative integer or null")
        if self.source_atom_map is not None:
            _positive_int(self.source_atom_map, "source_atom_map")
        if self.parent_atom_id is not None:
            _text(self.parent_atom_id, "parent_atom_id")
        if self.terminus not in {None, "head", "tail"}:
            raise ValueError("terminus must be head, tail, or null")
        if self.structural_role not in {"backbone", "sidechain", "terminus", "unassigned"}:
            raise ValueError("unsupported structural_role")
        if self.origin == "repeat":
            if self.repeat_index is None or self.source_atom_index is None or self.parent_atom_id is not None or self.terminus is not None:
                raise ValueError("repeat lineage requires copy/source indices and no cap or parent")
        elif self.origin == "cap":
            if self.terminus is None or self.repeat_index is not None or self.source_atom_index is not None or self.source_atom_map is not None or self.parent_atom_id is not None:
                raise ValueError("cap lineage requires a terminus and no repeat/source/parent")
        elif self.parent_atom_id is None or self.source_atom_index is not None or self.source_atom_map is not None:
            raise ValueError("added hydrogen lineage requires its parent and no source atom")

    @classmethod
    def from_mapping(cls, data):
        fields = set(cls.__dataclass_fields__)
        _mapping_keys(data, fields, "atom lineage")
        return cls(**data)

    def to_dict(self):
        return asdict(self)
