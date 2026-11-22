"""Finite-chain sequence enumeration and periodic basal-lattice controls."""

import random
import tempfile
import unittest
from collections import Counter
from importlib.util import find_spec
from itertools import combinations, pairwise
from math import sqrt
from pathlib import Path

from polychemmd.studio.geometry import build_geometry
from polychemmd.studio.microstructure import (
    acetate_sequence,
    conditional_sequence,
    correlated_sites,
    stream_seeds,
    surface_stats,
    triangular_edges,
)


def sequence(dp=12, count=3, **overrides):
    config = {
        "dp": dp,
        "hydrolysis": 1 - count / dp,
        "acetate_pattern": "dyads",
        "seed_streams": {"sequence": 37},
    }
    return acetate_sequence({**config, **overrides})


class FiniteSequenceTests(unittest.TestCase):
    def test_dp12_exact_multiplicities_and_direct_sampler(self):
        groups = Counter(
            sum(b == a + 1 for a, b in pairwise(s)) for s in combinations(range(12), 3)
        )
        self.assertEqual(groups, {0: 120, 1: 90, 2: 10})
        for d, expected in groups.items():
            _, info = sequence(acetate_dyads=d)
            self.assertEqual(info["sequence_count_at_dyads"], expected)
            produced = {
                tuple(conditional_sequence(12, 3, d, random.Random(i)))
                for i in range(2000)
            }
            self.assertEqual(len(produced), expected)
            self.assertTrue(
                all(sum(b == a + 1 for a, b in pairwise(s)) == d for s in produced)
            )

    def test_all_binary_sequences_up_to_dp9(self):
        for n in range(1, 10):
            for k in range(n + 1):
                all_counts = Counter(
                    sum(b == a + 1 for a, b in pairwise(s))
                    for s in combinations(range(n), k)
                )
                for d, total in all_counts.items():
                    residues, stats = sequence(n, k, acetate_dyads=d)
                    self.assertEqual(len(set(residues)), k)
                    self.assertEqual(stats["sequence_count_at_dyads"], total)
                    self.assertEqual(stats["acetate_dyads"], d)
                    self.assertEqual(
                        n - 1,
                        d + stats["OH_OH_dyads"] + stats["OH_OAc_dyads_both_orders"],
                    )

    def test_nmr_fraction_exact_composition_and_explicit_recipe(self):
        for value, achieved in [(0.15, 2), (0.0, 0), (0.8, 2)]:
            _, stats = sequence(acetate_dyad_fraction=value)
            self.assertAlmostEqual(stats["acetate_dyad_fraction"], achieved / 11)
            self.assertEqual(stats["requested_dyad_fraction"], value)
        r, stats = sequence(acetate_pattern="explicit", acetate_positions=[7, 5, 6])
        self.assertEqual(r, [4, 5, 6])
        self.assertEqual(stats["acetate_dyads"], 2)
        for settings in [
            {"acetate_dyads": True},
            {"acetate_dyads": 3},
            {"acetate_dyads": 0, "acetate_dyad_fraction": 0},
            {"acetate_dyad_fraction": float("nan")},
            {"acetate_pattern": "explicit", "acetate_positions": [5, 5, 6]},
            {"acetate_pattern": "explicit", "acetate_positions": [1, True, 12]},
        ]:
            with self.subTest(settings=settings), self.assertRaises(ValueError):
                sequence(**settings)
        with self.assertRaises(ValueError):
            sequence(1, 1, acetate_dyad_fraction=0)
        _, stats = sequence(1, 1, acetate_dyads=0)
        self.assertIsNone(stats["acetate_dyad_fraction"])

    def test_independent_stream_overrides(self):
        config = {"seed": 7, "pattern": "correlated"}
        seeds = stream_seeds(config)
        self.assertEqual(len(set(seeds.values())), 4)
        self.assertEqual(
            {**seeds, "surface": 9}, stream_seeds(dict(config, surface_seed=9))
        )
        self.assertEqual(stream_seeds({"seed": 7})["conformer"], 7)
        with self.assertRaises(ValueError):
            stream_seeds(dict(config, conformer_seed=True))


def lattice(nx, ny):
    sites = [
        (x + 0.5 * (y % 2), sqrt(3) / 2 * y) for x in range(nx) for y in range(2 * ny)
    ]
    return triangular_edges(sites, nx, ny, 1.0)


class SurfaceStatisticsTests(unittest.TestCase):
    def test_torus_edges_small_periodic_images_and_finite_null(self):
        for nx, ny in [(1, 1), (2, 1), (2, 2), (3, 2)]:
            n = 2 * nx * ny
            edges = lattice(nx, ny)
            degrees = Counter(i for pair in edges for i in pair)
            self.assertEqual(set(degrees.values()), {6})
            self.assertEqual(len(edges), 3 * n)
            stats = [
                surface_stats(set(s), edges, n) for s in combinations(range(n), n // 2)
            ]
            self.assertAlmostEqual(
                sum(s["oh_correlation"] for s in stats) / len(stats),
                stats[0]["fixed_composition_random_mean"],
            )

    def test_best_attainable_checked_against_exhaustive_complement_rich_and_poor(self):
        edges = lattice(3, 2)
        for k in [2, 4, 9, 11]:
            for target in [-1, -0.2, 0.25, 1]:
                true_error = min(
                    abs(surface_stats(set(c), edges, 12)["oh_correlation"] - target)
                    for c in combinations(range(12), k)
                )
                sites, s = correlated_sites(12, k, edges, target, 55)
                self.assertEqual(len(sites), k)
                self.assertAlmostEqual(s["absolute_correlation_error"], true_error)
                self.assertIn(s["status"], ("exact", "nearest_feasible_proven"))
        self.assertIsNone(surface_stats(set(), edges, 12)["oh_correlation"])
        with self.assertRaises(ValueError):
            correlated_sites(12, 12, edges, 0, 4)

    def test_large_torus_swap_updater_recounts(self):
        edges = lattice(6, 4)
        for expected in [-0.3, 0, 0.25, 1]:
            sites, s = correlated_sites(48, 24, edges, expected, 7)
            self.assertEqual(len(sites), 24)
            self.assertEqual(
                s["oh_oh_edges"], sum(a in sites and b in sites for a, b in edges)
            )
            self.assertEqual(
                s["total_edges"], s["oh_oh_edges"] + s["f_f_edges"] + s["oh_f_edges"]
            )
            self.assertEqual((sites, s), correlated_sites(48, 24, edges, expected, 7))


@unittest.skipUnless(
    find_spec("ase") and find_spec("rdkit"), "requires Studio chemistry extras"
)
class InverseGeometryTests(unittest.TestCase):
    def test_real_dp12_builder_dyads_and_four_faces(self):
        cfg = {
            "polymer": "PVA",
            "dp": 12,
            "hydrolysis": 0.75,
            "acetate_pattern": "dyads",
            "acetate_dyads": 1,
            "pattern": "correlated",
            "geometry": "slit",
            "nx": 3,
            "ny": 2,
            "surface_correlation_top": 0.3,
            "surface_correlation_bottom": -1 / 3,
            "conformer_seed": 17,
        }
        with tempfile.TemporaryDirectory() as temp:
            m = build_geometry(cfg, Path(temp))
            other = build_geometry({**cfg, "surface_seed": 777}, Path(temp))
        self.assertEqual(m["templates"]["polymer"], other["templates"]["polymer"])
        self.assertEqual(m["polymer"]["sequence_statistics"]["acetate_dyads"], 1)
        residues = m["polymer"]["acetate_residues_1based"]
        marked = {
            a["repeat_index"] + 1
            for a in m["templates"]["polymer"]
            if ":acetate:" in a["atom_id"]
        }
        self.assertEqual(set(residues), marked)
        self.assertEqual(len(m["face_occupancy"]), 2)
        for sheet in m["face_occupancy"]:
            for face, target in [("top", 0.3), ("bottom", -1 / 3)]:
                self.assertEqual(sheet[face]["n_OH"], 6)
                self.assertEqual(
                    sheet[face]["neighbor_statistics"]["requested_correlation"], target
                )


@unittest.skipUnless(find_spec("ase"), "requires Studio surface extra")
class VacuumPaddingTests(unittest.TestCase):
    def test_padding_changes_upper_cell_only(self):
        for shape in ["single", "slit"]:
            config = {"polymer": "bare", "geometry": shape, "nx": 2, "ny": 1}
            with tempfile.TemporaryDirectory() as tmp:
                before = build_geometry(config, Path(tmp))
                after = build_geometry({**config, "vacuum_padding_A": 15}, Path(tmp))
            self.assertEqual(before["atoms"], after["atoms"])
            self.assertEqual(before["accessible_z_A"], after["accessible_z_A"])
            self.assertEqual(
                before["accessible_volume_A3"], after["accessible_volume_A3"]
            )
            self.assertEqual(before["box_A"][2] + 15, after["box_A"][2])
            self.assertEqual(after["vacuum_padding_A"], 15)
