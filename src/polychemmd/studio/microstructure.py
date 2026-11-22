"""Finite, composition-preserving chain and surface microstructures."""

from __future__ import annotations

import math
import random
from hashlib import sha256
from itertools import combinations


def stream_seeds(config):
    """Legacy recipes replay their old coordinates. New controls use independent streams."""
    base = config["seed"]
    names = ("surface", "sequence", "tacticity", "conformer")
    independent = (
        config.get("pattern") == "correlated"
        or config.get("acetate_pattern") in ("dyads", "explicit")
        or any(f"{n}_seed" in config for n in names)
    )
    result = {}
    for name in names:
        key = name + "_seed"
        if key in config:
            value = config[key]
            if type(value) is not int or not 0 <= value < 2**31:
                raise ValueError(f"{key}: expected an integer in [0, 2^31-1]")
        elif independent:
            digest = sha256(f"polychemmd/interface/v1/{name}/{base}".encode()).digest()
            value = int.from_bytes(digest[:4], "big") % (2**31)
        else:
            value = base + (1777 if name == "tacticity" else 0)
        result[name] = value
    return result


def dyad_bounds(dp, count):
    if not count:
        return 0, 0
    return max(0, count - (dp - count + 1)), count - 1


def conditional_sequence(dp, count, dyads, rng):
    """Uniform exact draw for distinct chain ends; returns zero-based acetate indices."""
    lo, hi = dyad_bounds(dp, count)
    if type(dyads) is not int or not lo <= dyads <= hi:
        raise ValueError(
            f"acetate_dyads: finite chain admits integer counts {lo}..{hi}"
        )
    if not count:
        return []
    runs = count - dyads
    # An ordered positive composition of acetates and occupied gaps among the
    # intervening OH sites identify every admissible binary chain exactly once.
    cuts = [0, *sorted(rng.sample(range(1, count), runs - 1)), count]
    slots = sorted(rng.sample(range(dp - count + 1), runs))
    residues = []
    for j, slot in enumerate(slots):
        residues.extend(range(slot + cuts[j], slot + cuts[j + 1]))
    return residues


def acetate_sequence(config):
    dp = config["dp"]
    mode = config["acetate_pattern"]
    count = math.floor(dp * (1 - config["hydrolysis"]) + 0.5 + 1e-12)
    rng = random.Random(config["seed_streams"]["sequence"])
    lo, hi = dyad_bounds(dp, count)
    requested = config.get("acetate_dyads")
    fraction = config.get("acetate_dyad_fraction")
    if mode == "explicit":
        positions = config.get("acetate_positions")
        if (
            not isinstance(positions, list)
            or len(positions) != count
            or any(type(x) is not int or not 1 <= x <= dp for x in positions)
            or len(set(positions)) != count
        ):
            raise ValueError(
                f"acetate_positions must contain {count} distinct 1-based integers in 1..{dp}, consistent with hydrolysis"
            )
        residues = sorted(p - 1 for p in positions)
    elif mode == "dyads":
        if (requested is None) == (fraction is None):
            raise ValueError("provide either acetate_dyads or acetate_dyad_fraction")
        if fraction is not None:
            if (
                dp == 1
                or type(fraction) not in (int, float)
                or not math.isfinite(fraction)
                or not 0 <= fraction <= 1
            ):
                raise ValueError(
                    "acetate_dyad_fraction requires DP > 1 and a finite fraction of all (DP-1) dyads in [0,1]"
                )
            dyads = min(hi, max(lo, math.floor(fraction * (dp - 1) + 0.5 + 1e-12)))
        else:
            dyads = requested
        residues = conditional_sequence(dp, count, dyads, rng)
    elif mode == "random":
        residues = sorted(rng.sample(range(dp), count))
    else:
        start = rng.randrange(dp - count + 1) if count else 0
        residues = list(range(start, start + count))
    occupied = set(residues)
    counts = [0, 0, 0]
    for i in range(dp - 1):
        counts[int(i in occupied) + int(i + 1 in occupied)] += 1
    dyads = counts[2]
    runs = count - dyads
    multiplicity = (
        math.comb(count - 1, runs - 1) * math.comb(dp - count + 1, runs) if count else 1
    )
    achieved = dyads / (dp - 1) if dp > 1 else None
    return residues, {
        "acetate_dyads": dyads,
        "acetate_dyad_fraction": achieved,
        "OH_OH_dyads": counts[0],
        "OH_OAc_dyads_both_orders": counts[1],
        "all_dyads": dp - 1,
        "acetate_runs": runs,
        "mean_acetate_run_length": count / runs if runs else None,
        "sequence_count_at_dyads": multiplicity,
        "feasible_dyads": [lo, hi],
        "requested_dyads": requested,
        "requested_dyad_fraction": fraction,
        "dyad_fraction_error": achieved - fraction if fraction is not None else None,
        "target_mode": mode,
    }


def triangular_edges(sites, nx, ny, spacing):
    """Three forward bonds per site on a periodic six-neighbour triangular lattice."""
    width, height = 2 * nx, 2 * ny
    uy = spacing * math.sqrt(3) / 2
    keys = [(round(2 * x / spacing) % width, round(y / uy) % height) for x, y in sites]
    lookup = {key: i for i, key in enumerate(keys)}
    if len(lookup) != 2 * nx * ny or len(keys) != len(lookup):
        raise ValueError("unexpected triangular surface-site registry")
    return [
        (i, lookup[((x + dx) % width, (y + dy) % height)])
        for i, (x, y) in enumerate(keys)
        for dx, dy in ((2, 0), (1, 1), (-1, 1))
    ]


def _oo_edges(vertices, edges):
    return sum(a in vertices and b in vertices for a, b in edges)


def _correlation(occupied_edges, edge_count, fraction):
    if fraction in (0, 1):
        return None
    return (occupied_edges / edge_count - fraction**2) / (fraction * (1 - fraction))


def surface_stats(oh, edges, sites, target=None, status="observed", attempts=0):
    n, e = len(oh), len(edges)
    fraction = n / sites
    oo = _oo_edges(oh, edges)
    corr = _correlation(oo, e, fraction)
    loops = sum(a == b for a, b in edges)
    # A self-image is always correlated with itself. Distinct sites sampled
    # uniformly without replacement have correlation -1/(sites-1).
    null = (loops - (e - loops) / (sites - 1)) / e if corr is not None else None
    return {
        "total_edges": e,
        "oh_oh_edges": oo,
        "oh_f_edges": 6 * n - 2 * oo,
        "f_f_edges": e - 6 * n + oo,
        "self_image_edges": loops,
        "oh_correlation": corr,
        "fixed_composition_random_mean": null,
        "requested_correlation": target,
        "status": status,
        "search_evaluations": attempts,
        "absolute_correlation_error": abs(corr - target)
        if target is not None
        else None,
    }


def correlated_sites(sites, n, edges, target, seed):
    """Exact enumeration for small cells; bounded fixed-composition swaps otherwise."""
    if n == 0 or n == sites:
        raise ValueError(
            "a correlation target needs both OH and F on the realized face"
        )
    rng = random.Random(seed)
    f = n / sites
    desired = len(edges) * (f * f + target * f * (1 - f))
    labels = list(range(sites))
    rng.shuffle(labels)
    selected = set(labels[:n])
    pairs = _oo_edges(selected, edges)
    best = set(selected)
    best_error = abs(pairs - desired)
    evaluations = 1
    # All edge counts are integers, allowing early proof of the optimal error.
    lattice_bound = abs(desired - round(desired))
    exact_enum = math.comb(sites, n) <= 24000
    proven = lambda v: v <= lattice_bound + 1e-10
    if not proven(best_error) and exact_enum:
        # Enumerate the minority species; use degree-six bond counting to
        # recover OH/OH counts if the enumerated vertices are F.
        enumerate_oh = n <= sites // 2
        population = min(n, sites - n)
        incident = [[] for _ in labels]
        for i, j in edges:
            incident[i].append(j)
        for combo in combinations(labels, population):
            chosen = set(combo)
            inner = sum(j in chosen for i in chosen for j in incident[i])
            current = inner if enumerate_oh else len(edges) - 6 * population + inner
            evaluations += 1
            error = abs(current - desired)
            if error < best_error:
                best, best_error = (
                    chosen if enumerate_oh else set(labels) - chosen,
                    error,
                )
            if proven(best_error):
                break
    elif not proven(best_error):
        incident = [set() for _ in labels]
        for idx, (i, j) in enumerate(edges):
            incident[i].add(idx)
            incident[j].add(idx)
        on, off = labels[:n], labels[n:]
        budget = min(40000, max(4000, 80 * sites))
        for iteration in range(budget):
            i, j = rng.randrange(n), rng.randrange(sites - n)
            a, b = on[i], off[j]
            affected = incident[a] | incident[b]
            delta = 0
            for idx in affected:
                u, v = edges[idx]
                old = u in selected and v in selected
                new_u = u != a and (u in selected or u == b)
                new_v = v != a and (v in selected or v == b)
                delta += (new_u and new_v) - old
            cost = abs(pairs + delta - desired) - abs(pairs - desired)
            temp = 0.05 + 1.6 * (1 - iteration / budget) ** 2
            evaluations += 1
            if cost <= 0 or rng.random() < math.exp(-cost / temp):
                selected.remove(a)
                selected.add(b)
                on[i], off[j] = b, a
                pairs += delta
                error = abs(pairs - desired)
                if error < best_error:
                    best, best_error = set(selected), error
                    if proven(best_error):
                        break
    if best_error < 1e-10:
        status = "exact"
    elif proven(best_error) or exact_enum:
        status = "nearest_feasible_proven"
    else:
        status = "best_found_bounded_search"
    return best, surface_stats(best, edges, sites, target, status, evaluations)
