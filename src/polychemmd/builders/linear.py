"""Deterministic single-bond head-to-tail assembly of finite repeat-unit graphs.

The code assembles a supplied repeat unit, adds explicit H or methyl termini,
and validates atom balance. It does not transform a monomer, generate 3-D
coordinates, choose tacticity, or claim any simulation/physical observable.
"""

from collections import Counter
from dataclasses import dataclass, field
from typing import Any

from ..chemistry.contracts import (
    AtomLineage, RepeatUnitSpec, UnsupportedChemistryError, _positive_int, _rdkit,
)


def _isotope_counts(mol):
    return Counter((a.GetAtomicNum(), a.GetIsotope()) for a in mol.GetAtoms() if a.GetAtomicNum())


@dataclass(frozen=True)
class LinearChainResult:
    """Validated connectivity result. ``mol`` returns a defensive RDKit copy."""

    repeat_unit: RepeatUnitSpec
    degree_of_polymerization: int
    head_cap: str
    tail_cap: str
    canonical_smiles: str
    mapped_smiles: str
    molecular_formula: str
    formal_charge: int
    lineage: tuple[AtomLineage, ...]
    backbone_assignment_method: str
    backbone_source_atom_indices: tuple[int, ...]
    rdkit_version: str
    _mol: Any = field(repr=False, compare=False)

    @property
    def mol(self):
        return _rdkit().Mol(self._mol)

    def to_dict(self):
        counts = Counter(a.GetSymbol() for a in self._mol.GetAtoms())
        atoms = []
        for atom, lineage in zip(self._mol.GetAtoms(), self.lineage):
            atoms.append({
                **lineage.to_dict(), "element": atom.GetSymbol(),
                "isotope": atom.GetIsotope(), "formal_charge": atom.GetFormalCharge(),
                "chiral_tag": str(atom.GetChiralTag()),
            })
        bonds = [{
            "atom_maps": [b.GetBeginAtom().GetAtomMapNum(), b.GetEndAtom().GetAtomMapNum()],
            "order": b.GetBondTypeAsDouble(), "aromatic": b.GetIsAromatic(),
            "stereo": str(b.GetStereo()),
        } for b in self._mol.GetBonds()]
        return {
            "schema_version": "1", "kind": "finite_linear_chain_graph",
            "operation": "repeat_unit_assembly", "evidence_level": "graph_validated",
            "repeat_unit": self.repeat_unit.to_dict(),
            "degree_of_polymerization": self.degree_of_polymerization,
            "caps": {"head": self.head_cap, "tail": self.tail_cap},
            "linkage": {"direction": "tail_to_head", "bond_order": 1},
            "canonical_smiles": self.canonical_smiles,
            "mapped_smiles": self.mapped_smiles,
            "molecular_formula": self.molecular_formula,
            "formal_charge": self.formal_charge,
            "atom_count": self._mol.GetNumAtoms(),
            "heavy_atom_count": self._mol.GetNumHeavyAtoms(),
            "element_counts": dict(sorted(counts.items())),
            "hydrogen_representation": "all_explicit",
            "backbone_assignment": {
                "method": self.backbone_assignment_method,
                "source_atom_indices": list(self.backbone_source_atom_indices),
            },
            "atoms": atoms, "bonds": bonds,
            "validation": {
                "connected": True, "no_dummy_atoms": True,
                "atom_and_isotope_balance": True, "formal_charge_balance": True,
                "rdkit_sanitization": True,
                "geometry_generated": False, "force_field_assigned": False,
            },
            "software": {"rdkit": self.rdkit_version},
        }


def build_linear(
    repeat_unit: RepeatUnitSpec,
    degree_of_polymerization: int,
    head_cap: str = "H",
    tail_cap: str = "H",
    *,
    max_degree_of_polymerization: int = 1000,
    max_atoms: int = 100000,
) -> LinearChainResult:
    """Assemble identical copies with tail(i)–head(i+1) single bonds.

    The two neutral cap choices are exactly ``H`` and ``methyl``. All hydrogens
    are explicit in the returned graph. Safety limits can be raised explicitly;
    they are checked against the complete predicted atom count before assembly.
    """
    if not isinstance(repeat_unit, RepeatUnitSpec):
        raise ValueError("repeat_unit must be a validated RepeatUnitSpec")
    dp = _positive_int(degree_of_polymerization, "degree_of_polymerization")
    _positive_int(max_degree_of_polymerization, "max_degree_of_polymerization")
    _positive_int(max_atoms, "max_atoms")
    if dp > max_degree_of_polymerization:
        raise ValueError("degree_of_polymerization exceeds the configured safety limit")
    for name, value in (("head_cap", head_cap), ("tail_cap", tail_cap)):
        if not isinstance(value, str) or value not in {"H", "methyl"}:
            raise UnsupportedChemistryError(f"{name} must be H or methyl; other cap chemistry is unsupported")

    Chem = _rdkit()
    from rdkit import __version__ as rdkit_version
    from rdkit.Chem import rdMolDescriptors

    source = repeat_unit.to_mol()
    attachments = repeat_unit.attachment_indices(source)
    roles = repeat_unit.structural_roles(source)
    source_explicit = Chem.AddHs(source)
    expected_counts = _isotope_counts(source_explicit)
    expected_counts = Counter({key: count * dp for key, count in expected_counts.items()})
    for cap in (head_cap, tail_cap):
        if cap == "H":
            expected_counts[(1, 0)] += 1
        else:
            expected_counts[(6, 0)] += 1
            expected_counts[(1, 0)] += 3
    if sum(expected_counts.values()) > max_atoms:
        raise ValueError("predicted explicit atom count exceeds the configured max_atoms safety limit")

    # Keep full copied graphs until all junctions/caps have been inserted. Atom
    # properties carry source identities through RDKit's index renumbering.
    combined = None
    for repeat_index in range(dp):
        copy = Chem.Mol(source)
        for atom in copy.GetAtoms():
            atom.SetIntProp("_polychemmd_repeat_index", repeat_index)
            atom.SetIntProp("_polychemmd_source_index", atom.GetIdx())
            atom.SetIntProp("_polychemmd_source_map", atom.GetAtomMapNum())
        combined = copy if combined is None else Chem.CombineMols(combined, copy)
    working = Chem.RWMol(combined)
    source_size = source.GetNumAtoms()
    for repeat_index in range(dp - 1):
        working.AddBond(
            repeat_index * source_size + attachments["tail_anchor"],
            (repeat_index + 1) * source_size + attachments["head_anchor"],
            Chem.BondType.SINGLE,
        )
    for terminus, cap, anchor in (
        ("head", head_cap, attachments["head_anchor"]),
        ("tail", tail_cap, (dp - 1) * source_size + attachments["tail_anchor"]),
    ):
        atom = Chem.Atom(1 if cap == "H" else 6)
        atom.SetProp("_polychemmd_terminus", terminus)
        cap_index = working.AddAtom(atom)
        working.AddBond(anchor, cap_index, Chem.BondType.SINGLE)
    for atom_index in sorted(
        (i * source_size + attachments[end] for i in range(dp) for end in ("head_dummy", "tail_dummy")),
        reverse=True,
    ):
        working.RemoveAtom(atom_index)
    mol = working.GetMol()
    try:
        Chem.SanitizeMol(mol)
        mol = Chem.AddHs(mol)
        Chem.SanitizeMol(mol)
        for atom in mol.GetAtoms():
            atom.SetAtomMapNum(0)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    except (ValueError, RuntimeError) as exc:
        raise ValueError("assembled chain fails chemical valence/stereochemistry sanitization") from exc
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetAtomicNum() == 0 for a in mol.GetAtoms()):
        raise ValueError("assembled chain must be connected and contain no attachment dummies")
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("assembled chain contains unsupported radicals")
    if _isotope_counts(mol) != expected_counts:
        raise ValueError("assembled chain failed exact atom/isotope balance")
    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    expected_charge = dp * sum(a.GetFormalCharge() for a in source.GetAtoms())
    if charge != expected_charge:
        raise ValueError("assembled chain failed formal-charge balance")
    for atom in mol.GetAtoms():
        if atom.HasProp("_polychemmd_source_index"):
            original = source.GetAtomWithIdx(atom.GetIntProp("_polychemmd_source_index"))
            if (original.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
                    and atom.GetChiralTag() != original.GetChiralTag()):
                raise UnsupportedChemistryError(
                    "specified internal stereochemistry is lost or changes under this assembly/capping; "
                    "use a supported stereochemically distinct graph"
                )

    # Existing atom origins are independent of output SMILES canonicalization.
    # Added H atoms are appended by AddHs and inherit their parent atom's origin.
    lineage = []
    for atom in mol.GetAtoms():
        if atom.HasProp("_polychemmd_repeat_index"):
            repeat_index = atom.GetIntProp("_polychemmd_repeat_index")
            source_index = atom.GetIntProp("_polychemmd_source_index")
            original_map = atom.GetIntProp("_polychemmd_source_map")
            item = AtomLineage(
                f"repeat:{repeat_index}:atom:{source_index}", atom.GetIdx(), atom.GetIdx() + 1,
                "repeat", repeat_index, source_index, original_map or None,
                None, None, roles["atom_roles"][source_index],
            )
        elif atom.HasProp("_polychemmd_terminus"):
            terminus = atom.GetProp("_polychemmd_terminus")
            item = AtomLineage(
                f"cap:{terminus}:{atom.GetSymbol()}", atom.GetIdx(), atom.GetIdx() + 1,
                "cap", None, None, None, None, terminus, "terminus",
            )
        else:
            # Hydrogens added by AddHs have exactly one pre-existing parent.
            if atom.GetAtomicNum() != 1 or atom.GetDegree() != 1:
                raise ValueError("unexpected untracked atom after hydrogen completion")
            parent = atom.GetNeighbors()[0]
            parent_record = lineage[parent.GetIdx()]
            hydrogens = [a.GetIdx() for a in parent.GetNeighbors() if not a.HasProp("_polychemmd_repeat_index") and not a.HasProp("_polychemmd_terminus")]
            hydrogen_index = sorted(hydrogens).index(atom.GetIdx())
            item = AtomLineage(
                f"{parent_record.atom_id}:H:{hydrogen_index}", atom.GetIdx(), atom.GetIdx() + 1,
                "added_hydrogen", parent_record.repeat_index, None, None,
                parent_record.atom_id, parent_record.terminus, parent_record.structural_role,
            )
        atom.SetAtomMapNum(item.atom_map)
        atom.SetProp("_polychemmd_atom_id", item.atom_id)
        lineage.append(item)
    if len({a.atom_id for a in lineage}) != mol.GetNumAtoms():
        raise ValueError("assembled chain contains duplicate atom identities")

    # Mapped graph and canonical molecular identity serve different purposes.
    # Clear maps before removing H so cap-origin maps cannot inhibit removal.
    identity = Chem.Mol(mol)
    for atom in identity.GetAtoms():
        atom.SetAtomMapNum(0)
    identity = Chem.RemoveHs(identity)
    canonical_smiles = Chem.MolToSmiles(identity, canonical=True, isomericSmiles=True)
    mapped_smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True, allHsExplicit=True)
    return LinearChainResult(
        repeat_unit=repeat_unit, degree_of_polymerization=dp,
        head_cap=head_cap, tail_cap=tail_cap,
        canonical_smiles=canonical_smiles, mapped_smiles=mapped_smiles,
        molecular_formula=rdMolDescriptors.CalcMolFormula(mol, separateIsotopes=True), formal_charge=charge,
        lineage=tuple(lineage), backbone_assignment_method=roles["method"],
        backbone_source_atom_indices=tuple(roles["backbone_source_atom_indices"]),
        rdkit_version=rdkit_version, _mol=mol,
    )
