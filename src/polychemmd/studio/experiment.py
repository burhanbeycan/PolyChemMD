"""Paired optical and lamellar measurements, with source-resolved literature records."""

from __future__ import annotations

import csv
import io
import json
import math
import statistics
from collections import defaultdict
from functools import lru_cache
from importlib.resources import files
from pathlib import Path

PARK = "10.1126/sciadv.abl5299"
XING = "10.1002/eom2.12300"
SOURCES = {
    PARK: "Park et al. · Science Advances (2022), Fig. 2B and Fig. 1E",
    XING: "Xing et al. · EcoMat (2023), Table 2",
    "10.1039/D4SC04935G": "Di et al. · Chemical Science (2025), Fig. S26",
    "10.1557/s43578-025-01588-6": "Fotius et al. · Journal of Materials Research (2025), Fig. 2b",
    "10.1039/D5NR00450K": "Naqvi et al. · Nanoscale (2025), Fig. 1b",
    "10.1002/advs.202500953": "Poliukhova et al. · Advanced Science (2025), Fig. 2e",
}
XING_ROW = dict(
    batch="Xing Table 2",
    bare_dry_nm=1.346,
    coated_dry_nm=1.406,
    bare_wet_nm=1.527,
    coated_wet_nm=1.415,
)
DRY_WET = ("bare_dry_nm", "coated_dry_nm", "bare_wet_nm", "coated_wet_nm")


def _read(name: str) -> list[dict]:
    data = files("polychemmd.studio").joinpath("data", name).read_text(encoding="utf-8")
    return list(csv.DictReader(io.StringIO(data)))


@lru_cache(maxsize=1)
def atlas() -> dict:
    """Read audited, numerical subset of the Atlas (all records are secondary use)."""
    measured = []
    for row in _read("park_fig2b.csv"):
        medium = row["medium"]
        condition = "Published SW mixture" if medium == "SW" else f"1 M {medium}"
        measured.append(
            dict(
                id=f"park-{row['polymer']}-{medium}",
                polymer=row["polymer"],
                method="A750",
                condition=condition,
                value=float(row["absorbance_au"]),
                unit="AU",
                source_doi=PARK,
                panel="Fig. 2B",
                state="freeze-dried then redispersed",
                nominal_polymer_wt_pct=int(row["nominal_polymer_wt_pct"]),
                wavelength_nm=750,
                reading_bounds=[
                    float(row["digitization_lo_au"]),
                    float(row["digitization_hi_au"]),
                ],
                note="Absolute plotted absorbance. Reading interval reflects digitization. Measurement time and replicate uncertainty unreported.",
            )
        )
    for row in _read("xrd_pairs.csv"):
        measured.append(
            dict(
                id=row["id"],
                polymer=row["polymer"],
                method="XRD Δd002",
                condition=row["composition"],
                value=float(row["delta_nm"]),
                unit="nm",
                source_doi=row["doi"],
                panel=row["fig"],
                state=row["specimen"],
                bare_nm=float(row["baseline_nm"]),
                coated_nm=float(row["loaded_nm"]),
                note=row["qualified_note"],
            )
        )
    measured.append(
        dict(
            id="XING-WET",
            polymer="polyethyleneimine",
            method="XRD Δd002",
            condition="PEI:MXene 0.5:1 by mass in feed",
            value=-0.112,
            unit="nm",
            source_doi=XING,
            panel="Table 2",
            state="PAN-supported, water-immersed",
            bare_nm=XING_ROW["bare_wet_nm"],
            coated_nm=XING_ROW["coated_wet_nm"],
            note="Paired to dry measurements in this same study. One study-level observation.",
        )
    )
    for row in _read("di_loading.csv"):
        measured.append(
            dict(
                id=f"di-{row['sample']}",
                polymer=row["sample"],
                method="TGA",
                condition=row["protocol"],
                value=float(row["loading_wt_percent"]),
                unit="wt%",
                source_doi=row["doi"],
                panel="Fig. S26 and article text",
                state="washed composite",
                note="Same mass feed for PDDA and PAH, separate feed for PIL. "
                + row["uncertainty"],
            )
        )
    return {
        "records": measured,
        "sources": [
            dict(doi=k, title=v, url="https://doi.org/" + k) for k, v in SOURCES.items()
        ],
    }


def _text(row: dict, key: str) -> str:
    s = row.get(key)
    if not isinstance(s, str) or not 1 <= len(s.strip()) <= 120:
        raise ValueError(f"{key} must be a non-empty label (max 120 characters)")
    return s.strip()


def _number(row: dict, key: str, positive: bool = False) -> float:
    x = row.get(key)
    if isinstance(x, bool) or not isinstance(x, (int, float)) or not math.isfinite(x):
        raise ValueError(f"{key} must be a finite number")
    if positive and x <= 0:
        raise ValueError(f"{key} must be positive")
    return float(x)


def _summary(values):
    values = list(values)
    if not all(math.isfinite(v) for v in values):
        raise ValueError("Calculated value outside finite range")
    return dict(
        n_pairs=len(values),
        mean=statistics.mean(values) if values else None,
        sample_sd=statistics.stdev(values) if len(values) > 1 else None,
    )


def run_experiment(config: dict, out: Path) -> dict:
    """Apply matched-condition equations. A built-in reference is never labeled user data."""
    assay = config.get("assay")
    source = config.get("reference")
    if source is not None and (
        source != "xing-pei" or assay != "xrd" or "rows" in config
    ):
        raise ValueError(
            "The named XRD reference and user observations are mutually exclusive"
        )
    rows = [XING_ROW] if source else config.get("rows")
    if (
        not isinstance(rows, list)
        or not 1 <= len(rows) <= 1000
        or not all(isinstance(r, dict) for r in rows)
    ):
        raise ValueError("Supply 1–1000 observation rows")
    observed, seen = [], set()
    if assay == "xrd":
        for r in rows:
            batch = _text(r, "batch")
            if batch in seen:
                raise ValueError("Duplicate batch within paired XRD observations")
            seen.add(batch)
            a, b, c, d = (_number(r, k, True) for k in DRY_WET)
            observed.append(
                {
                    "batch": batch,
                    **dict(zip(DRY_WET, (a, b, c, d))),
                    "dry_effect_nm": b - a,
                    "wet_effect_nm": d - c,
                    "wet_minus_dry_nm": (d - c) - (b - a),
                }
            )
        summary = _summary(r["wet_minus_dry_nm"] for r in observed)
        summary.update(
            unit="nm", comparison="(coated − bare) wet minus (coated − bare) dry"
        )
    elif assay == "optical":
        pairs = defaultdict(dict)
        for r in rows:
            batch, medium, polymer = (
                _text(r, k) for k in ("batch", "medium", "polymer")
            )
            key = (batch, medium, polymer)
            if key in seen:
                raise ValueError(
                    "Combine technical readings first. Duplicate sample within batch and medium."
                )
            seen.add(key)
            a0, a24, b0, b24 = (_number(r, k) for k in ("A0", "A24", "B0", "B24"))
            baseline = a0 - b0
            if baseline <= 0:
                raise ValueError(
                    f"Non-positive corrected baseline for {batch}, {medium}, {polymer}"
                )
            retention = (a24 - b24) / baseline
            if not math.isfinite(retention):
                raise ValueError("Corrected R24 is outside the finite range")
            observed.append(
                dict(
                    batch=batch,
                    medium=medium,
                    polymer=polymer,
                    A0=a0,
                    A24=a24,
                    B0=b0,
                    B24=b24,
                    corrected_A0=baseline,
                    R24=retention,
                )
            )
            pairs[batch, medium][polymer] = retention
        by_medium = defaultdict(list)
        missing = []
        for (batch, medium), measures in pairs.items():
            if {"PCE", "PEO"} <= measures.keys():
                by_medium[medium].append(
                    dict(batch=batch, difference=measures["PCE"] - measures["PEO"])
                )
            else:
                missing.append(dict(batch=batch, medium=medium))
        summary = dict(
            unit="dimensionless",
            comparison="PCE − PEO, paired within batch and medium",
            by_medium={
                m: {**_summary(x["difference"] for x in group), "paired": group}
                for m, group in by_medium.items()
            },
            unpaired_batch_media=missing,
        )
    else:
        raise ValueError("assay must be xrd or optical")
    result = dict(
        assay=assay,
        evidence="published secondary analysis"
        if source
        else "user-supplied observations",
        source_doi=XING if source else None,
        rows=observed,
        summary=summary,
    )
    # Preserve raw, corrected and paired values alongside the human-readable view.
    out.mkdir(parents=True, exist_ok=True)
    (out / "experiment.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False),
        encoding="utf-8",
    )
    with (out / "experiment.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=observed[0].keys())
        writer.writeheader()
        for row in observed:
            # A user-supplied batch identifier can otherwise become a spreadsheet formula.
            writer.writerow(
                {
                    k: (
                        "'" + v
                        if isinstance(v, str) and v.startswith(("=", "+", "-", "@"))
                        else v
                    )
                    for k, v in row.items()
                }
            )
    return result
