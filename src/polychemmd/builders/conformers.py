"""Seeded distance-geometry coordinates for already validated finite graphs.

This layer uses RDKit's ETKDGv3 backend, with independently written lineage,
geometry and stereochemical checks. It performs no force-field assignment,
energy ranking, reaction simulation, or molecular dynamics.
"""

from __future__ import annotations

import hashlib
import json
import math
import platform

from .. import __version__
from ..chemistry.contracts import UnsupportedChemistryError, _rdkit
from ..geometry import GeometryAtom, GeometryBond, MolecularGeometry
from .linear import LinearChainResult


_MAX_RDKIT_INT = 2**31 - 1
_EMBED_TIMEOUT_SECONDS = 30


class ConformerGenerationError(ValueError):
    """Embedding or its independent geometric checks did not succeed."""


def _integer(value, name, *, minimum=1):
    if type(value) is not int or not minimum <= value <= _MAX_RDKIT_INT:
        raise ValueError(f"{name} must be an integer in [{minimum}, {_MAX_RDKIT_INT}]")
    return value


def _canonical_json(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _distance(first, second):
    return math.sqrt(math.fsum((a - b)**2 for a, b in zip(first, second)))


def _verify_stereochemistry(source, embedded, conformer_id=0):
    """Reconstruct specified stereo from coordinates after clearing copied tags.

    Atom ordering remains fixed. Atom-map numbers are provenance rather than
    chemical substituents, and are cleared on both comparison copies so that
    equivalent explicit H atoms do not gain artificial CIP distinctions.
    """
    Chem = _rdkit()
    conformer = embedded.GetConformer(conformer_id)
    if not conformer.Is3D():
        raise ConformerGenerationError("embedded conformer must be marked as 3D")
    expected, actual = Chem.Mol(source), Chem.Mol(embedded)
    for molecule in (expected, actual):
        for atom in molecule.GetAtoms():
            atom.SetAtomMapNum(0)
    Chem.AssignStereochemistry(expected, cleanIt=True, force=True)
    specified_atoms = {
        atom.GetIdx(): atom.GetChiralTag()
        for atom in expected.GetAtoms()
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    }
    allowed = {Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW}
    if any(tag not in allowed for tag in specified_atoms.values()):
        raise UnsupportedChemistryError("coordinate verification supports only specified tetrahedral atom stereo")
    defined_bonds = {Chem.BondStereo.STEREOE, Chem.BondStereo.STEREOZ,
                     Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOTRANS}
    specified_bonds = {
        bond.GetIdx(): bond.GetStereo()
        for bond in expected.GetBonds() if bond.GetStereo() in defined_bonds
    }
    ignored_bond_states = {Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY}
    if any(bond.GetStereo() not in defined_bonds | ignored_bond_states for bond in expected.GetBonds()):
        raise UnsupportedChemistryError("coordinate verification supports only specified E/Z double-bond stereo")

    # Existing atom tags, bond stereo and directional bonds must not count as
    # evidence that generated coordinates have the requested handedness.
    Chem.RemoveStereochemistry(actual)
    Chem.AssignStereochemistryFrom3D(actual, confId=conformer_id, replaceExistingTags=True)
    Chem.AssignStereochemistry(actual, cleanIt=True, force=True)
    for index, tag in specified_atoms.items():
        if actual.GetAtomWithIdx(index).GetChiralTag() != tag:
            raise ConformerGenerationError(f"generated coordinates do not preserve specified tetrahedral stereo at atom {index}")
    for index, stereo in specified_bonds.items():
        if actual.GetBondWithIdx(index).GetStereo() != stereo:
            raise ConformerGenerationError(f"generated coordinates do not preserve specified double-bond stereo at bond {index}")
    return {
        "specified_tetrahedral_centers_checked": len(specified_atoms),
        "specified_double_bonds_checked": len(specified_bonds),
        "method": "clear_atom_maps_and_stereo_on_probe_then_assign_from_3d",
        "unspecified_stereochemistry": "not_assigned_to_source_graph",
    }


def _check_embedded_geometry(source, embedded, conformer_id=0):
    Chem = _rdkit()
    if embedded.GetNumConformers() != 1 or embedded.GetNumAtoms() != source.GetNumAtoms():
        raise ConformerGenerationError("embedding must return one conformer without changing atom count")
    conformer = embedded.GetConformer(conformer_id)
    positions = tuple(tuple(float(value) for value in conformer.GetAtomPosition(i))
                      for i in range(embedded.GetNumAtoms()))
    if any(not math.isfinite(value) for point in positions for value in point):
        raise ConformerGenerationError("embedding returned nonfinite coordinates")
    closest = math.inf
    for i, first in enumerate(positions):
        for second in positions[:i]:
            distance = _distance(first, second)
            if not math.isfinite(distance) or distance < 0.35:
                raise ConformerGenerationError("embedding returned coincident or implausibly close atom positions")
            closest = min(closest, distance)

    periodic_table = Chem.GetPeriodicTable()
    bond_lengths = []
    for bond in embedded.GetBonds():
        left, right = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        distance = _distance(positions[left], positions[right])
        covalent_sum = math.fsum(periodic_table.GetRcovalent(embedded.GetAtomWithIdx(i).GetAtomicNum())
                                 for i in (left, right))
        if not math.isfinite(covalent_sum) or covalent_sum <= 0:
            raise UnsupportedChemistryError("a positive covalent-radius reference is required for bonded geometry checks")
        # Broad rejection limits only, not an equilibrium bond-length model.
        if not 0.55 * covalent_sum <= distance <= 1.65 * covalent_sum:
            raise ConformerGenerationError(f"generated bond {bond.GetIdx()} fails broad covalent-distance sanity bounds")
        bond_lengths.append(distance)
    if not bond_lengths:
        raise ConformerGenerationError("a connected molecular graph with bonds is required")
    stereo = _verify_stereochemistry(source, embedded, conformer_id)
    return positions, {
        "coordinates_finite": True,
        "all_atom_separation_floor_angstrom": 0.35,
        "minimum_atom_separation_angstrom": closest,
        "bond_distance_check": "0.55_to_1.65_times_sum_of_RDKit_covalent_radii",
        "minimum_bond_length_angstrom": min(bond_lengths),
        "maximum_bond_length_angstrom": max(bond_lengths),
        "stereochemistry": stereo,
    }


def generate_conformer(
    chain: LinearChainResult,
    *,
    seed: int,
    max_iterations: int = 1000,
    max_atoms: int = 2000,
) -> MolecularGeometry:
    """Generate one checked ETKDGv3 realization, retaining all explicit atoms.

    ``max_iterations`` bounds RDKit embedding iterations; its backend timeout
    is 30 seconds. This is not a hard subprocess deadline or a conformational
    search. Identical coordinates are expected only for the same graph, seed,
    options and compatible RDKit/runtime build.
    """
    if not isinstance(chain, LinearChainResult):
        raise ValueError("chain must be a LinearChainResult")
    _integer(seed, "seed", minimum=0)
    _integer(max_iterations, "max_iterations")
    _integer(max_atoms, "max_atoms")
    Chem = _rdkit()
    from rdkit import __version__ as rdkit_version
    from rdkit.Chem import rdDistGeom

    source = chain.mol
    if source.GetNumAtoms() > max_atoms:
        raise ValueError("chain atom count exceeds max_atoms for conformer generation")
    graph = chain.to_dict()
    if len(graph["atoms"]) != source.GetNumAtoms() or len(chain.lineage) != source.GetNumAtoms():
        raise ValueError("chain graph and atomic lineage have incompatible sizes")
    for index, record in enumerate(graph["atoms"]):
        atom = source.GetAtomWithIdx(index)
        if (record["atom_id"] != chain.lineage[index].atom_id
                or record["atom_map"] != atom.GetAtomMapNum()
                or record["element"] != atom.GetSymbol()
                or record["isotope"] != atom.GetIsotope()
                or record["formal_charge"] != atom.GetFormalCharge()):
            raise ValueError("chain atomic identity and lineage do not agree")
    source.RemoveAllConformers()
    embedded = Chem.Mol(source)
    identity = Chem.MolToSmiles(source, canonical=True, isomericSmiles=True, allHsExplicit=True)
    parameters = rdDistGeom.ETKDGv3()
    parameters.randomSeed = seed
    parameters.numThreads = 1
    parameters.enforceChirality = True
    parameters.maxIterations = max_iterations
    parameters.timeout = _EMBED_TIMEOUT_SECONDS
    parameters.clearConfs = True
    parameters.trackFailures = True
    try:
        conformer_id = rdDistGeom.EmbedMolecule(embedded, parameters)
    except (ValueError, RuntimeError, OverflowError) as exc:
        raise ConformerGenerationError("RDKit ETKDGv3 coordinate generation failed") from exc
    if conformer_id < 0:
        raise ConformerGenerationError(
            "RDKit ETKDGv3 did not produce a conformer within its configured limits; "
            f"failure counts={tuple(parameters.GetFailureCounts())}"
        )
    if Chem.MolToSmiles(embedded, canonical=True, isomericSmiles=True, allHsExplicit=True) != identity:
        raise ConformerGenerationError("embedding unexpectedly changed the mapped chemical graph")
    positions, validation = _check_embedded_geometry(source, embedded, conformer_id)
    records = graph["atoms"]
    atoms = tuple(GeometryAtom(
        atom_id=record["atom_id"], atom_map=record["atom_map"],
        element=record["element"], isotope=record["isotope"],
        formal_charge=record["formal_charge"],
        mass_g_mol=source.GetAtomWithIdx(index).GetMass(),
        lineage_json=_canonical_json(record),
    ) for index, record in enumerate(records))
    bonds = tuple(GeometryBond(
        atom_ids=(records[b.GetBeginAtomIdx()]["atom_id"], records[b.GetEndAtomIdx()]["atom_id"]),
        order=b.GetBondTypeAsDouble(), aromatic=b.GetIsAromatic(), stereo=str(b.GetStereo()),
    ) for b in source.GetBonds())
    provenance = {
        "method": "RDKit_ETKDGv3",
        "source_kind": "generated_molecular_geometry",
        "seed": seed,
        "software": {"rdkit": rdkit_version, "polychemmd": __version__,
                     "python": platform.python_version()},
        "parameters": {
            "max_iterations": max_iterations,
            "backend_timeout_seconds": _EMBED_TIMEOUT_SECONDS,
            "num_threads": 1,
            "enforce_chirality": True,
            "ETversion": parameters.ETversion,
            "use_macrocycle_torsions": parameters.useMacrocycleTorsions,
            "use_small_ring_torsions": parameters.useSmallRingTorsions,
        },
        "validation": validation,
        "source_graph_identity": {
            "canonical_smiles": chain.canonical_smiles,
            "mapped_smiles": chain.mapped_smiles,
            "molecular_formula": chain.molecular_formula,
            "formal_charge": chain.formal_charge,
            "degree_of_polymerization": chain.degree_of_polymerization,
            "caps": {"head": chain.head_cap, "tail": chain.tail_cap},
        },
        "force_field_assigned": False,
        "energy_minimization_performed": False,
        "energy_minimization_scope": "no separate post-embedding force-field minimization; ETKDG internal objective optimization is part of embedding",
        "conformational_search_performed": False,
        "simulation_performed": False,
        "equilibration_status": "not_performed",
        "material_property_validation": "not_performed",
    }
    return MolecularGeometry(
        atoms=atoms, bonds=bonds, coordinates_angstrom=positions,
        source_graph_sha256=hashlib.sha256(_canonical_json(graph).encode("utf-8")).hexdigest(),
        provenance_json=_canonical_json(provenance),
    )
