from importlib.util import find_spec
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from polychemmd.studio.ensemble import masses_for_lengths, number_weights, run_ensemble
from polychemmd.studio.geometry import build_geometry


@unittest.skipUnless(
    find_spec("scipy") and find_spec("rdkit"), "requires the Studio chemistry stack"
)
class TestPVAEnsemble(unittest.TestCase):
    def test_discrete_extreme_distributions(self):
        masses = [100, 200, 300]
        for mn, d, expected in [
            (200, 1, [0, 1, 0]),
            (200, 1.25, [0.5, 0, 0.5]),
            (150, 1 + 2500 / 150**2, [0.5, 0.5, 0]),
        ]:
            self.assertEqual(number_weights(masses, mn, d).tolist(), expected)
        for mn, d in [(150, 1), (200, 2), (350, 1.5)]:
            with self.assertRaises(ValueError):
                number_weights(masses, mn, d)

    def test_interior_unique_flat_prior_and_molecular_weight_convention(self):
        m = [100, 150, 200, 250, 300]
        # A uniform number distribution realizes mean 200 and variance 5000.
        p = number_weights(m, 200, 1.125)
        for weight in p:
            self.assertAlmostEqual(weight, 0.2, places=10)
        expected, _ = masses_for_lengths([8], 0.75)
        self.assertAlmostEqual(expected[0], 452.541, places=3)

    @unittest.skipUnless(
        find_spec("ase"), "requires ASE to reconstruct the selected chain"
    )
    def test_queue_recipe_to_graph_mass_and_dyad(self):
        import json

        with TemporaryDirectory() as tmp:
            folder = Path(tmp)
            spec = dict(mn_g_mol=600, dispersity=1.12, hydrolysis=0.75, seed=17)
            result = run_ensemble(spec, folder)
            self.assertAlmostEqual(result["fitted_mn_g_mol"], 600.0, places=7)
            self.assertAlmostEqual(result["fitted_dispersity"], 1.12, places=8)
            self.assertEqual(result, run_ensemble(spec, folder))
            selected = result["members"][0]
            recipe = {
                k: selected[k]
                for k in ("dp", "acetate_positions", "sequence_seed", "hydrolysis")
            }
            model = build_geometry(
                {**recipe, "polymer": "PVA", "acetate_pattern": "explicit"}, folder
            )
            self.assertAlmostEqual(
                model["polymer"]["molar_mass_g_mol"], selected["mass_g_mol"], places=7
            )
            self.assertEqual(
                model["polymer"]["sequence_statistics"]["acetate_dyads"],
                selected["acetate_dyads"],
            )
            self.assertTrue((folder / "number-distribution.csv").is_file())
            self.assertTrue(json.dumps(result, allow_nan=False))

    def test_no_large_chain_clipping(self):
        with TemporaryDirectory() as tmp, self.assertRaisesRegex(ValueError, "outside"):
            run_ensemble(dict(mn_g_mol=50000, dispersity=1.2), Path(tmp))
