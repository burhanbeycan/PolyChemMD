"""Fit a number distribution of finite PVA oligomers and draw explicit sequences."""

import csv
import math
from pathlib import Path
import random

from .geometry import _integer, _real
from .microstructure import acetate_sequence


def masses_for_lengths(lengths, hydrolysis):
    from rdkit import Chem

    weight = Chem.GetPeriodicTable().GetAtomicWeight
    # H and methyl caps used by the atomistic PVA builder, plus C2H2O per acetylation.
    cap = weight("C") + 4 * weight("H")
    alcohol = 2 * weight("C") + 4 * weight("H") + weight("O")
    acetyl = 2 * weight("C") + 2 * weight("H") + weight("O")
    counts = [math.floor(n * (1 - hydrolysis) + 0.5 + 1e-12) for n in lengths]
    return [cap + n * alcohol + k * acetyl for n, k in zip(lengths, counts)], counts


def number_weights(masses, mn, dispersity):
    """Maximum-entropy number fractions on the supplied discrete mass support."""
    import numpy as np
    from scipy.optimize import brentq
    from scipy.special import softmax

    m = np.asarray(masses, dtype=float)
    if not m[0] <= mn <= m[-1]:
        raise ValueError(
            f"Requested Mn is outside the supported oligomer mass range {m[0]:g}..{m[-1]:g} g/mol."
        )
    index = np.searchsorted(m, mn)
    low = max(0, index - 1)
    high = min(index, len(m) - 1)
    lower = (mn - m[low]) * (m[high] - mn)
    upper = (mn - m[0]) * (m[-1] - mn)
    variance = (dispersity - 1) * mn**2
    scale = float(m[-1] - m[0])
    tol = max(1e-10, 2e-12 * scale**2)
    if variance < lower - tol or variance > upper + tol:
        raise ValueError(
            f"At this Mn and DP range, attainable dispersity is {1 + lower / mn**2:.9g}..{1 + upper / mn**2:.9g}."
        )

    def two_points(a, b):
        p = np.zeros_like(m)
        if a == b:
            p[a] = 1
        else:
            p[b] = (mn - m[a]) / (m[b] - m[a])
            p[a] = 1 - p[b]
        return p

    if variance <= lower + tol:
        return two_points(low, high)
    if variance >= upper - tol:
        return two_points(0, len(m) - 1)
    x = (m - mn) / scale
    target = variance / scale**2

    def root(fn):
        bound = 1.0
        while fn(-bound) > 0 or fn(bound) < 0:
            bound *= 2
            if bound > 2**48:
                raise ValueError(
                    "Ill-conditioned moment constraints at the edge of the DP range."
                )
        return brentq(fn, -bound, bound, xtol=1e-12)

    def profile(b):
        a = root(lambda a: float(x @ softmax(a * x + b * x * x)))
        return softmax(a * x + b * x * x)

    b = root(lambda b: float((x * x) @ profile(b)) - target)
    return profile(b)


def run_ensemble(config, out):
    mn = _real(config.get("mn_g_mol"), "mn_g_mol", 1, 1e9)
    dispersity = _real(config.get("dispersity"), "dispersity", 1, 1e7)
    hydrolysis = _real(config.get("hydrolysis", 0.75), "hydrolysis", 0, 1)
    dyad = _real(
        config.get("acetate_dyad_fraction", 0.08), "acetate_dyad_fraction", 0, 1
    )
    dpmin = _integer(config.get("dp_min", 2), "dp_min", 2, 79)
    dpmax = _integer(config.get("dp_max", 40), "dp_max", dpmin + 1, 80)
    sample = _integer(config.get("n_chains", 12), "n_chains", 1, 128)
    seed = _integer(config.get("seed", 17), "seed", 0, 2147483646)
    lengths = list(range(dpmin, dpmax + 1))
    masses, counts = masses_for_lengths(lengths, hydrolysis)
    weights = number_weights(masses, mn, dispersity).tolist()
    fitted = sum(m * w for m, w in zip(masses, weights))
    d_fit = sum(m * m * w for m, w in zip(masses, weights)) / fitted**2
    if not math.isclose(fitted, mn, rel_tol=1e-8) or not math.isclose(
        d_fit, dispersity, rel_tol=1e-8
    ):
        raise ValueError(
            "The fitted distribution did not reproduce both molecular-weight moments."
        )
    rng = random.Random(seed)
    indices = rng.choices(range(len(lengths)), weights=weights, k=sample)
    members = []
    for i in indices:
        n = lengths[i]
        sequence_seed = rng.randrange(2**31)
        residues, stats = acetate_sequence(
            dict(
                dp=n,
                hydrolysis=hydrolysis,
                acetate_pattern="dyads",
                acetate_dyad_fraction=dyad,
                seed_streams={"sequence": sequence_seed},
            )
        )
        members.append(
            dict(
                dp=n,
                hydrolysis=hydrolysis,
                sequence_seed=sequence_seed,
                mass_g_mol=masses[i],
                acetate_positions=[r + 1 for r in residues],
                acetate_dyads=stats["acetate_dyads"],
            )
        )
    out = Path(out)
    out.mkdir(exist_ok=True, parents=True)
    with (out / "number-distribution.csv").open("w", newline="", encoding="utf-8") as f:
        table = csv.writer(f)
        table.writerow(["DP", "acetate_count", "model_mass_g_mol", "number_fraction"])
        table.writerows(zip(lengths, counts, masses, weights))
    msample = sum(masses[i] for i in indices) / sample
    n_units = sum(lengths[i] for i in indices)
    return dict(
        distribution="Maximum entropy with uniform prior over supported DP values; number-weighted finite oligomer model",
        fitted_mn_g_mol=fitted,
        fitted_dispersity=d_fit,
        fitted_mw_g_mol=fitted * d_fit,
        sample_mn_g_mol=msample,
        sample_dispersity=sum(masses[i] ** 2 for i in indices) / sample / msample**2,
        requested_hydrolysis=hydrolysis,
        sample_hydrolysis=1 - sum(counts[i] for i in indices) / n_units,
        target_acetate_dyad_fraction=dyad,
        sample_acetate_dyad_fraction=sum(c["acetate_dyads"] for c in members)
        / (n_units - sample),
        dp_range=[dpmin, dpmax],
        n_chains=sample,
        seed=seed,
        members=members,
        end_groups="H / methyl as in Studio PVA",
        fractions_basis="all adjacent repeat pairs; finite-chain rounding per chain",
    )
