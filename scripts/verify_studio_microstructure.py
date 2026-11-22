"""Build matched PVA–Ti3C2Tx geometries and check the prescribed sequence and sites."""

import argparse
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from polychemmd.studio.geometry import build_geometry


def verify():
    root = Path(__file__).resolve().parents[1]
    config = json.loads(
        (root / "examples/studio/controlled_interface.json").read_text()
    )
    rows = []
    with TemporaryDirectory() as temp:
        out = Path(temp)
        for target in (-1 / 3, 0.25):
            for dyads in (0, 1, 2):
                request = dict(
                    config, surface_correlation_top=target, acetate_dyads=dyads
                )
                m = build_geometry(request, out)
                p = m["polymer"]
                f = [s["top"]["neighbor_statistics"] for s in m["face_occupancy"]]
                if p["sequence_statistics"]["acetate_dyads"] != dyads:
                    raise RuntimeError("dyad mismatch")
                for face, s in zip(m["face_occupancy"], f, strict=True):
                    if (
                        face["top"]["n_OH"] != 24
                        or abs(s["oh_correlation"] - target) > 1e-8
                    ):
                        raise RuntimeError("surface composition/correlation mismatch")
                rows.append(
                    {
                        "requested_acetate_dyads": dyads,
                        "requested_top_face_correlation": target,
                        "per_sheet_top_face_correlation": [
                            s["oh_correlation"] for s in f
                        ],
                        "acetate_positions_1based": p["acetate_residues_1based"],
                        "sequence_count_at_dyads": p["sequence_statistics"][
                            "sequence_count_at_dyads"
                        ],
                        "top_face_OH_count_per_sheet": [
                            s["top"]["n_OH"] for s in m["face_occupancy"]
                        ],
                        "atom_count": len(m["atoms"]),
                        "molecular_formula": p["molecular_formula"],
                        "xyz_sha256": m["xyz_sha256"],
                    }
                )
    return {
        "method": "Matched finite microstructure controls, Å coordinates. No dynamics sampled.",
        "recipe": config,
        "ase_version": m["ase_version"],
        "rdkit_version": p["rdkit_version"],
        "cases": rows,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
