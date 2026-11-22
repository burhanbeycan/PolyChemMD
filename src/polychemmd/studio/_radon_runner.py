"""Reference calculation executed only inside the external RadonPy interpreter."""

import hashlib
import importlib.metadata
import importlib.util
import json
import math
from pathlib import Path
import sys
from datetime import datetime, timezone


def capabilities():
    from radonpy.core import poly  # noqa: F401 — probe the actual graph builder
    from radonpy.ff.gaff2 import GAFF2  # noqa: F401 — probe the parameterizer
    from rdkit import rdBase

    return dict(
        engine="RadonPy",
        version=importlib.metadata.version("radonpy-pypi"),
        rdkit=rdBase.rdkitVersion,
        python=sys.version.split()[0],
        psi4_detected=importlib.util.find_spec("psi4") is not None,
        workflow="Linear oligomer, seeded ETKDG and MMFF94 optimization, GAFF2 typing audit",
        citation="10.1038/s41524-022-00906-4",
    )


def write(name, data):
    Path(name).write_text(
        json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )


def run(request):
    from radonpy.core import poly
    from radonpy.ff.gaff2 import GAFF2
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolDescriptors

    versions = capabilities()
    monomers = {
        "PEO": ("*CCO*", 3),
        "PVA": ("*CC(O)*", 3),
        "PVP": ("*CC(N1CCCC1=O)*", 8),
        "PEtOx": ("*CCN(C(=O)CC)*", 7),
    }
    s, heavy_per_repeat = monomers[request["polymer"]]
    # These are finite alkyl-capped model molecules. PEO here is explicitly not a PE glycol diol.
    mol = poly.polymerize_MolFromSmiles(s, n=request["dp"], terminal="C")
    if mol is None:
        raise RuntimeError("RadonPy did not construct the molecular graph")
    Chem.SanitizeMol(mol)
    if (
        mol.GetNumHeavyAtoms() != request["dp"] * heavy_per_repeat + 2
        or Chem.GetFormalCharge(mol) != 0
    ):
        raise ValueError("Unexpected capped graph size or charge")
    params = AllChem.ETKDGv3()
    params.randomSeed = request["seed"]
    params.numThreads = 1
    params.maxIterations = 150
    using_random_coords = False
    if AllChem.EmbedMolecule(mol, params) < 0:
        using_random_coords = True
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) < 0:
            raise RuntimeError("Could not embed the oligomer")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("MMFF94 missing atom parameters")
    code = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94", maxIters=1200)
    if code != 0:
        raise RuntimeError(
            "MMFF94 geometry minimization did not converge within 1200 iterations"
        )
    if not GAFF2().ff_assign(mol, charge=None):
        raise ValueError("RadonPy GAFF2 assignment failed")
    conf = mol.GetConformer()
    atoms = []
    for atom in mol.GetAtoms():
        index = atom.GetIdx()
        pos = conf.GetAtomPosition(index)
        if not all(math.isfinite(v) for v in (pos.x, pos.y, pos.z)):
            raise ValueError("Nonfinite atom coordinates")
        atom.SetAtomMapNum(index + 1)
        atoms.append(
            dict(
                element=atom.GetSymbol(),
                x=pos.x,
                y=pos.y,
                z=pos.z,
                role="polymer",
                atom_map=index + 1,
                gaff2=atom.GetProp("ff_type"),
            )
        )
    bonds = [
        [b.GetBeginAtomIdx(), b.GetEndAtomIdx(), b.GetBondTypeAsDouble()]
        for b in mol.GetBonds()
    ]
    structure = {
        "atoms": atoms,
        "bonds": bonds,
        "units": "angstrom",
        "pbc": [False, False, False],
    }
    write("structure.json", structure)
    Chem.MolToMolFile(mol, "polymer.sdf")
    write("radon_input.json", request)
    result = dict(
        engine=versions["engine"],
        software=versions,
        request=request,
        repeat_smiles=s,
        formula=rdMolDescriptors.CalcMolFormula(mol),
        n_atoms=len(atoms),
        heavy_atoms=mol.GetNumHeavyAtoms(),
        gaff2_typed=len(atoms),
        stereochemistry="unassigned repeat stereocenters, one embedded realization",
        caps="methyl substituent at each repeat-unit attachment",
        charge_model="no numerical partial charges assigned",
        random_coordinate_retry=using_random_coords,
        evidence="isolated neutral model oligomer; geometry generation and atom typing",
        source_doi=versions["citation"],
        utc=datetime.now(timezone.utc).isoformat(),
        input_sha256=hashlib.sha256(Path("radon_input.json").read_bytes()).hexdigest(),
        sdf_sha256=hashlib.sha256(Path("polymer.sdf").read_bytes()).hexdigest(),
    )
    write("radon_result.json", result)


if __name__ == "__main__":
    if sys.argv[1:] == ["version"]:
        print(json.dumps(capabilities()))
    elif sys.argv[1:] == ["oligomer"]:
        run(json.load(sys.stdin))
    else:
        raise SystemExit("Select version or oligomer")
