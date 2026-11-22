"""Build starting Ti3C2Tx basal surfaces and explicitly mapped polymer templates.

Lengths are in ångströms. The substrate is an idealized geometric construction,
with neither a force field nor a claim that the mixed surfaces are equilibrated.
"""

from __future__ import annotations

import json
import math
import random
from hashlib import sha256
from pathlib import Path

from ..builders.linear import build_linear
from ..chemistry.contracts import RepeatUnitSpec
from .microstructure import (
    acetate_sequence,
    correlated_sites,
    stream_seeds,
    surface_stats,
    triangular_edges,
)

# Triangular in-plane spacing and approximate normal positions relative to
# central Ti. The first length is |a|/18 of the native 18-cell interface.
# For z: median position in the upper face (14 < z < 19.5 Å) of the
# upper sheet in Armstrong's deposited Stack3x3_H2Osys_OH_47.data minus
# the median central Ti position. Flatten the layers; recover the OH bond length with an in-plane component.
A = 54.8927 / 18
C_Z, TI_Z, O_Z, F_Z, H_Z = 1.18934, 2.41927, 3.46238, 3.48834, 4.36813
SURFACE_REFERENCE = {
    "doi": "10.1021/acsnano.5c20046",
    "member": "Stack3x3_H2Osys_OH_47.data",
    "member_sha256": "ec194edb8b6a5d1923269a1dee4bdf2ba047964ad930a5d52a887baa112fa873",
    "site_reference": "10.1021/acsomega.2c05913, section 3.1 and figure 1",
    "model": "Planar median layers; occupancy on CCP hollows over central Ti. "
    "A common basal lattice is imposed at all requested OH/F fractions. "
    "OH normal separation is 0.90575 Å and O-H length 0.94 Å (native bond type 9); "
    "its azimuth is seeded among six lattice directions.",
}
# Transformed two-ended units, independently assembled with PolyChemMD.
REPEATS = {
    "PEO": "[*:1]CCO[*:2]",
    "PVA": "[*:1]CC(O)[*:2]",
    "PVP": "[*:1]CC(N1CCCC1=O)[*:2]",
    "PEtOx": "[*:1]N(C(=O)CC)CC[*:2]",
    "PEI": "[*:1]NCC[*:2]",
    "PCE": "[*:1]CC(C(=O)O)[*:2]",
}


def _integer(value, name, low, high):
    if type(value) is not int or not low <= value <= high:
        raise ValueError(f"{name}: expected an integer from {low} to {high}")
    return value


def _real(value, name, low, high):
    if (
        type(value) not in (int, float)
        or not math.isfinite(value)
        or not low <= value <= high
    ):
        raise ValueError(f"{name}: expected a finite number from {low} to {high}")
    return float(value)


def _choice(value, name, choices):
    if value not in choices:
        raise ValueError(f"{name}: choose one of {', '.join(choices)}")
    return value


def _nearest_count(value):
    """Round nonnegative occupancy with half-up ties, independent of float quirks."""
    return math.floor(value + 0.5 + 1e-12)


def _select_sites(sites, n, size, mode, seed):
    """Either uniform sampling without replacement or a compact minimum-image patch."""
    rng = random.Random(seed)
    labels = list(range(len(sites)))
    if mode == "random":
        return set(rng.sample(labels, n))
    cx, cy = sites[rng.randrange(len(sites))]

    def squared_distance(i):
        dx, dy = sites[i][0] - cx, sites[i][1] - cy
        return (dx - size[0] * round(dx / size[0])) ** 2 + (
            dy - size[1] * round(dy / size[1])
        ) ** 2

    return set(sorted(labels, key=squared_distance)[:n])


def _sheet(nx, ny, fractions, mode, seed, sheet_index, targets):
    import ase
    import numpy as np
    from ase import Atoms
    from ase.build import make_supercell

    # Primitive (60-degree) hex basis. The lower half is the inversion partner
    # of the upper. The determinant-two transform makes a rectangular basal cell.
    ay = math.sqrt(3) * A
    primitive = Atoms(
        ["Ti", "C", "Ti", "C", "Ti"],
        positions=[
            (A / 2, ay / 6, -TI_Z),
            (A, ay / 3, -C_Z),
            (0, 0, 0),
            (A / 2, ay / 6, C_Z),
            (A, ay / 3, TI_Z),
        ],
        cell=[[A, 0, 0], [A / 2, ay / 2, 0], [0, 0, 20]],
        pbc=[True, True, False],
    )
    # compound primitive-to-orthogonal and in-plane tiling transforms
    cell = make_supercell(
        primitive, np.array([[nx, 0, 0], [-ny, 2 * ny, 0], [0, 0, 1]]), wrap=False
    )
    lx, ly = nx * A, ny * ay
    atoms = []
    sites = []

    def add(el, x, y, z, layer, site=None):
        atoms.append(
            dict(
                element=el,
                x=float(x % lx),
                y=float(y % ly),
                z=float(z),
                role="slab",
                layer=layer,
                site=site,
                sheet_index=sheet_index,
                molecule_id=f"slab-{sheet_index}",
            )
        )

    for symbol, (x, y, z) in zip(cell.get_chemical_symbols(), cell.positions):
        # Avoid representing a coordinate numerically close to zero as L.
        x, y = round(x, 10) % lx, round(y, 10) % ly
        if math.isclose(x, lx, abs_tol=1e-8):
            x = 0.0
        if math.isclose(y, ly, abs_tol=1e-8):
            y = 0.0
        layer = (
            "middle_Ti" if abs(z) < 1e-8 else ("top_" if z > 0 else "bottom_") + symbol
        )
        add(symbol, x, y, z, layer)
        if layer == "middle_Ti":
            sites.append((x, y))
    # Each surface has one termination site for every primitive Ti3C2 unit.
    edges = triangular_edges(sites, nx, ny, A)
    face_summary = {}
    for sign, face in ((1, "top"), (-1, "bottom")):
        fraction = fractions[face]
        n = _nearest_count(fraction * len(sites))
        site_seed = seed + 65537 * sheet_index + (1 if sign == 1 else 2)
        if mode == "correlated":
            oh, stats = correlated_sites(len(sites), n, edges, targets[face], site_seed)
        else:
            oh = _select_sites(sites, n, (lx, ly), mode, site_seed)
            stats = surface_stats(oh, edges, len(sites))
        oh_rng = random.Random(seed + 65537 * sheet_index + 1901 * sign)
        radial = math.sqrt(0.94**2 - (H_Z - O_Z) ** 2)
        for i, (x, y) in enumerate(sites):
            if i in oh:
                add("O", x, y, sign * O_Z, face + "_OH", i)
                phi = oh_rng.randrange(6) * math.pi / 3
                add(
                    "H",
                    x + radial * math.cos(phi),
                    y + radial * math.sin(phi),
                    sign * H_Z,
                    face + "_OH",
                    i,
                )
            else:
                add("F", x, y, sign * F_Z, face + "_F", i)
        face_summary[face] = dict(
            requested_OH=fraction,
            n_OH=n,
            n_F=len(sites) - n,
            realized_OH=n / len(sites),
            site_seed=site_seed,
            neighbor_statistics=stats,
        )
    return atoms, face_summary, ase.__version__


def _polymer_template(config, out, center):
    from rdkit import Chem, rdBase
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

    from ..builders.conformers import _check_embedded_geometry

    p, dp, seed = config["polymer"], config["dp"], config["seed_streams"]["conformer"]
    head = "methyl" if p in ("PEtOx", "PEI") else "H"
    tail = "methyl" if p in ("PVA", "PVP", "PCE") else "H"
    chain = build_linear(RepeatUnitSpec(p, REPEATS[p], 1, 2), dp, head, tail)
    mol = Chem.RWMol(chain.mol)
    lineages = [(r.atom_id, r.repeat_index) for r in chain.lineage]
    residues = []
    if p == "PEO":
        # Change the carbon-side H cap to OH, giving HO-(CH2CH2O)n-H.
        i = next(i for i, a in enumerate(chain.lineage) if a.atom_id == "cap:head:H")
        mol.ReplaceAtom(i, Chem.Atom("O"))
        lineages[i] = ("cap:head:OH:O", None)
    if p == "PVA":
        residues, sequence_statistics = acetate_sequence(config)
        for i, r in enumerate(chain.lineage):
            if r.repeat_index in residues and r.source_atom_index == 3:
                h = next(
                    a.GetIdx()
                    for a in mol.GetAtomWithIdx(i).GetNeighbors()
                    if a.GetAtomicNum() == 1
                )
                mol.ReplaceAtom(h, Chem.Atom("C"))
                ox = mol.AddAtom(Chem.Atom("O"))
                methyl = mol.AddAtom(Chem.Atom("C"))
                mol.AddBond(h, ox, Chem.BondType.DOUBLE)
                mol.AddBond(h, methyl, Chem.BondType.SINGLE)
                label = f"repeat:{r.repeat_index}:acetate"
                lineages[h] = (label + ":C", r.repeat_index)
                lineages.extend(
                    [(label + ":O", r.repeat_index), (label + ":CH3", r.repeat_index)]
                )
    grafts = []
    if p == "PCE":
        # Explicit acid-form AA backbone with a periodic series of neutral
        # methoxy-oligo(ethylene glycol) acrylate grafts. This is a specified
        # comb oligomer, not an assignment for an unspecified commercial PCE.
        grafts = list(range(0, dp, config["graft_every"]))
        for i, r in enumerate(chain.lineage):
            if r.repeat_index not in grafts or r.source_atom_index != 5:
                continue
            h = next(
                a.GetIdx()
                for a in mol.GetAtomWithIdx(i).GetNeighbors()
                if a.GetAtomicNum() == 1
            )
            label = f"repeat:{r.repeat_index}:graft"
            mol.ReplaceAtom(h, Chem.Atom("C"))
            lineages[h] = (label + ":0", r.repeat_index)
            end = h
            for k, element in enumerate(("CCO" * config["side_dp"] + "C")[1:], 1):
                idx = mol.AddAtom(Chem.Atom(element))
                mol.AddBond(end, idx, Chem.BondType.SINGLE)
                lineages.append((f"{label}:{k}", r.repeat_index))
                end = idx
    mol = mol.GetMol()
    for atom in mol.GetAtoms():
        atom.SetAtomMapNum(0)
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    Chem.SanitizeMol(mol)
    for i in range(len(lineages), mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(i)
        parent = atom.GetNeighbors()[0].GetIdx()
        lineages.append((f"{lineages[parent][0]}:newH:{i}", lineages[parent][1]))

    # Source SMILES leave tacticity unspecified. Select one finite atactic
    # realization with a distinct seed stream and explicitly report its graph.
    rng = random.Random(config["seed_streams"]["tacticity"])
    if p in ("PVA", "PVP", "PCE"):
        for i, cip in Chem.FindMolChiralCenters(
            mol, includeUnassigned=True, useLegacyImplementation=False
        ):
            if cip == "?" and mol.GetAtomWithIdx(i).GetSymbol() == "C":
                mol.GetAtomWithIdx(i).SetChiralTag(
                    rng.choice(
                        [
                            Chem.ChiralType.CHI_TETRAHEDRAL_CCW,
                            Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                        ]
                    )
                )
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(mol)) != 1 or Chem.GetFormalCharge(mol) != 0:
        raise ValueError(
            "polymer template must be a connected, charge-neutral molecule"
        )
    conf = Chem.Mol(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    params.timeout = 25
    if AllChem.EmbedMolecule(conf, params) < 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(conf, params) < 0:
            raise ValueError("ETKDGv3 could not generate the requested chain")
    ff_name = config["conformer_ff"]
    if ff_name == "MMFF94" and AllChem.MMFFHasAllMoleculeParams(conf):
        ff = AllChem.MMFFGetMoleculeForceField(
            conf, AllChem.MMFFGetMoleculeProperties(conf)
        )
    elif ff_name == "UFF" and AllChem.UFFHasAllMoleculeParams(conf):
        ff = AllChem.UFFGetMoleculeForceField(conf)
    else:
        raise ValueError(f"incomplete {ff_name} parameters for this isolated oligomer")
    if ff is None or ff.Minimize(maxIts=3000):
        raise ValueError(f"{ff_name} minimization did not converge for seed {seed}")
    coords, checks = _check_embedded_geometry(mol, conf)
    formula = rdMolDescriptors.CalcMolFormula(conf)
    smiles = Chem.MolToSmiles(Chem.RemoveHs(conf), isomericSmiles=True)
    for a in conf.GetAtoms():
        a.SetAtomMapNum(a.GetIdx() + 1)
    # Placement translates the optimized conformer; atom maps and internal
    # coordinates are retained both in the JSON preview and in the SDF.
    shift = [
        target - sum(pt[axis] for pt in coords) / len(coords)
        for axis, target in enumerate(center)
    ]
    coords = [[v + d for v, d in zip(pt, shift)] for pt in coords]
    for i, point in enumerate(coords):
        conf.GetConformer().SetAtomPosition(i, point)
    conf.SetProp("_Name", f"{p} DP{dp}")
    with Chem.SDWriter(str(out / "polymer.sdf")) as writer:
        writer.write(conf)
    atoms = [
        dict(
            element=a.GetSymbol(),
            x=float(x),
            y=float(y),
            z=float(z),
            role="polymer",
            atom_map=i + 1,
            atom_id=lineages[i][0],
            repeat_index=lineages[i][1],
            formal_charge=a.GetFormalCharge(),
            molecule_id="polymer-0",
        )
        for i, (a, (x, y, z)) in enumerate(zip(conf.GetAtoms(), coords))
    ]
    bonds = [
        [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
        for b in conf.GetBonds()
    ]
    info = dict(
        name=p,
        degree_of_polymerization=dp,
        formal_charge=0,
        molecular_formula=formula,
        molar_mass_g_mol=float(Descriptors.MolWt(conf)),
        isomeric_smiles=smiles,
        head_cap="OH" if p == "PEO" else head,
        tail_cap=tail,
        stereochemistry="seeded atactic realization"
        if p in ("PVA", "PVP", "PCE")
        else "achiral",
        acetate_residues_1based=[r + 1 for r in residues],
        realized_hydrolysis=1 - len(residues) / dp if p == "PVA" else None,
        force_field_scope="single isolated oligomer",
        conformer_method=f"ETKDGv3 / {ff_name}",
        minimized_energy_kJ_mol=float(ff.CalcEnergy()) * 4.184,
        geometry_checks=checks,
        rdkit_version=rdBase.rdkitVersion,
    )
    if p == "PVA":
        info["sequence_statistics"] = sequence_statistics
    if p == "PCE":
        info.update(
            architecture="neutral acrylic acid / methoxy-PEG acrylate model comb",
            ionization_state="free acid; no pH or commercial-grade assignment",
            side_chain_EO_units=config["side_dp"],
            graft_residues_1based=[r + 1 for r in grafts],
        )
    if p == "PEI":
        info.update(
            architecture="linear polyethyleneimine",
            ionization_state="neutral free base",
        )
    return atoms, bonds, info


def build_geometry(config: dict, out: Path) -> dict:
    """Write an idealized ASE substrate and a mapped neutral RDKit chain template.

    The initial gap is an outer-nucleus geometric thickness. Accessible volume
    is the geometric envelope, without van der Waals or polymer-volume correction.
    The preview displays one template. The packing stage places requested copies.
    """
    if not isinstance(config, dict):
        raise ValueError("geometry requires a config object")
    c = dict(config)
    shape = _choice(c.get("geometry", "single"), "geometry", ("single", "slit"))
    nx = _integer(c.get("nx", 6), "nx", 1, 30)
    ny = _integer(c.get("ny", 4), "ny", 1, 30)
    c["seed"] = _integer(c.get("seed", 17), "seed", 0, 2**31 - 1)
    c["dp"] = _integer(c.get("dp", 8), "dp", 1, 80)
    c["polymer"] = _choice(c.get("polymer", "PEO"), "polymer", (*REPEATS, "bare"))
    c["hydrolysis"] = _real(c.get("hydrolysis", 1.0), "hydrolysis", 0, 1)
    c["acetate_pattern"] = _choice(
        c.get("acetate_pattern", "random"),
        "acetate_pattern",
        ("random", "clustered", "explicit", "dyads"),
    )
    if c["acetate_pattern"] != "explicit" and "acetate_positions" in c:
        raise ValueError("acetate_positions requires acetate_pattern=explicit")
    if c["acetate_pattern"] != "dyads" and (
        {"acetate_dyads", "acetate_dyad_fraction"} & c.keys()
    ):
        raise ValueError("acetate dyad targets require acetate_pattern=dyads")
    c["conformer_ff"] = _choice(
        c.get("conformer_ff", "MMFF94"), "conformer_ff", ("MMFF94", "UFF")
    )
    c["side_dp"] = _integer(c.get("side_dp", 2), "side_dp", 1, 12)
    c["graft_every"] = _integer(c.get("graft_every", 4), "graft_every", 1, 48)
    if c["polymer"] != "PCE" and (c["side_dp"] != 2 or c["graft_every"] != 4):
        raise ValueError("side_dp and graft_every apply only to PCE")
    if (
        c["polymer"] == "PCE"
        and c["dp"] + (1 + (c["dp"] - 1) // c["graft_every"]) * c["side_dp"] > 80
    ):
        raise ValueError(
            "PCE backbone DP + total side-chain EO units exceeds the 80 unit construction limit"
        )
    if c["polymer"] != "PVA" and (
        c["hydrolysis"] != 1 or c["acetate_pattern"] != "random"
    ):
        raise ValueError("hydrolysis and acetate pattern apply only to PVA")
    n_chains = _integer(
        c.get("chains", 0 if c["polymer"] == "bare" else 1), "chains", 0, 96
    )
    if (c["polymer"] == "bare") != (n_chains == 0):
        raise ValueError(
            "choose bare with zero chains or a polymer with at least one chain"
        )
    fractions = {
        face: _real(c.get("oh_" + face, 0.5), "oh_" + face, 0, 1)
        for face in ("top", "bottom")
    }
    pattern = _choice(
        c.get("pattern", "random"), "pattern", ("random", "patch", "correlated")
    )
    targets = {}
    for face in fractions:
        key = "surface_correlation_" + face
        if pattern == "correlated":
            targets[face] = _real(c.get(key), key, -1, 1)
        elif key in c:
            raise ValueError(f"{key} requires pattern=correlated")
    c["seed_streams"] = stream_seeds(c)
    slit_gap = _real(c.get("gap_A", 22), "gap_A", 2, 200)
    water_height = _real(c.get("water_height_A", 22), "water_height_A", 2, 200)
    vacuum = _real(c.get("vacuum_padding_A", 0), "vacuum_padding_A", 0, 100)
    sheets = []
    centers = []
    summaries = []
    try:
        for i in range(2 if shape == "slit" else 1):
            atoms, summary, ase_version = _sheet(
                nx, ny, fractions, pattern, c["seed_streams"]["surface"], i, targets
            )
            target_bottom = 2.0 if i == 0 else max(a["z"] for a in sheets) + slit_gap
            center = target_bottom - min(a["z"] for a in atoms)
            for a in atoms:
                a["z"] += center
            sheets.extend(atoms)
            centers.append(center)
            summaries.append(summary)
        out = Path(out)
        out.mkdir(parents=True, exist_ok=True)
    except ImportError as exc:
        raise RuntimeError(
            f"Missing geometry dependency ({exc.name}); install the Studio environment"
        ) from exc
    surface_top = max(a["z"] for a in sheets if a["sheet_index"] == 0)
    lo = surface_top
    hi = (
        lo + water_height
        if shape == "single"
        else min(a["z"] for a in sheets if a["sheet_index"] == 1)
    )
    box = [
        nx * A,
        ny * math.sqrt(3) * A,
        (hi if shape == "single" else max(a["z"] for a in sheets) + 2.0) + vacuum,
    ]
    try:
        chain, bonds, info = (
            _polymer_template(c, out, [box[0] / 2, box[1] / 2, (lo + hi) / 2])
            if n_chains
            else ([], [], {"name": "bare", "molar_mass_g_mol": 0.0})
        )
    except ImportError as exc:
        raise RuntimeError(
            f"Missing geometry dependency ({exc.name}); install the Studio environment"
        ) from exc
    all_atoms = sheets + chain
    manifest = dict(
        schema_version=1,
        units="angstrom",
        geometry=shape,
        box_A=box,
        pbc=[True, True, False],
        atoms=all_atoms,
        bonds=[[i + len(sheets), j + len(sheets), k] for i, j, k in bonds],
        templates={"slab": sheets, "polymer": chain},
        polymer_bonds=bonds,
        polymer=info,
        requested_chains=n_chains,
        sheet_centers_A=centers,
        accessible_z_A=[lo, hi],
        vacuum_padding_A=vacuum,
        accessible_volume_A3=(hi - lo) * box[0] * box[1],
        volume_convention="geometric outer-nuclei gap; excludes substrate only",
        center_spacing_A=centers[1] - centers[0] if len(centers) == 2 else None,
        face_occupancy=summaries,
        surface_reference=SURFACE_REFERENCE,
        ase_version=ase_version,
        coordinates_status="initial substrate and one isolated template; simulation not performed",
        fixed_geometry=True,
        site_pattern=pattern,
        topology_scope="polymer covalent bonds only; substrate is a geometry template",
        seed=c["seed"],
        seed_streams=c["seed_streams"],
    )
    xyz = f"{len(all_atoms)}\nASE Ti3C2Tx and optional isolated oligomer preview | angstrom\n"
    xyz += "".join(
        f"{a['element']} {a['x']:.9f} {a['y']:.9f} {a['z']:.9f}\n" for a in all_atoms
    )
    (out / "structure.xyz").write_text(xyz, encoding="utf-8")
    manifest["xyz_sha256"] = sha256(xyz.encode()).hexdigest()
    (out / "structure.json").write_text(
        json.dumps(manifest, allow_nan=False), encoding="utf-8"
    )
    return manifest
