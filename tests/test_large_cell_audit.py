"""Independent synthetic audit checks, not simulated material data."""

from copy import deepcopy
import importlib.util
from pathlib import Path
from types import SimpleNamespace
import unittest


SCRIPT = Path(__file__).resolve().parents[1]/"scripts/verify_large_cell_v07.py"
SPEC = importlib.util.spec_from_file_location("polychemmd_large_auditor_test", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def singlepoint():
    energies = {key: 1. for key in AUDIT.COMPONENTS[:-1]}
    energies["total"] = 6.
    return {"atom_ids": ["synthetic:1", "synthetic:2"], "assignment_canonical_sha256": "a"*64,
            "energy_components_kj_mol": energies, "forces_kj_mol_angstrom": [[1., 2., 3.], [-1., -2., -3.]]}


class LargeAuditSensitivityTests(unittest.TestCase):
    def test_exact_synthetic_cutoff_difference_and_per_atom_accounting(self):
        original = singlepoint(); variant = deepcopy(original)
        variant["energy_components_kj_mol"]["lj"] += .02
        variant["energy_components_kj_mol"]["total"] += .02
        variant["forces_kj_mol_angstrom"][0][0] += .03
        result = AUDIT.compare_singlepoints(original, variant, "cutoff_12")
        self.assertAlmostEqual(result["total_energy_difference_per_atom_kj_mol"], .01)
        self.assertAlmostEqual(result["max_force_component_difference_kj_mol_angstrom"], .03)
        self.assertAlmostEqual(result["rms_force_component_difference_kj_mol_angstrom"], .03/(6**.5))
        self.assertEqual(result["review_status"], "diagnostics_passed")
        self.assertEqual(result["bulk_property_convergence"], "not_established")
        self.assertIn("change the truncated Hamiltonian", result["interpretation"])

    def test_review_failure_retained_without_adjusting_limits(self):
        original = singlepoint(); variant = deepcopy(original)
        variant["energy_components_kj_mol"]["lj"] += 1.
        variant["energy_components_kj_mol"]["total"] += 1.
        result = AUDIT.compare_singlepoints(original, variant, "cutoff_12")
        self.assertEqual(result["review_status"], "requires_review")
        self.assertEqual(result["predeclared_development_review_limits"], AUDIT.SENSITIVITY_REVIEW_LIMITS["cutoff_12"])

    def test_reciprocal_change_cannot_change_realspace_model(self):
        original = singlepoint(); variant = deepcopy(original)
        variant["energy_components_kj_mol"]["coulomb_kspace"] += 1e-6
        variant["energy_components_kj_mol"]["total"] += 1e-6
        self.assertEqual(AUDIT.compare_singlepoints(original, variant, "reciprocal_20")["review_status"], "diagnostics_passed")
        variant["energy_components_kj_mol"]["coulomb_realspace"] += .001
        variant["energy_components_kj_mol"]["total"] += .001
        with self.assertRaisesRegex(AUDIT.LargeCellAuditError, "fixed model component"):
            AUDIT.compare_singlepoints(original, variant, "reciprocal_20")

    def test_mismatched_assignment_or_atom_order_rejected(self):
        original = singlepoint()
        for key, value in (("atom_ids", ["synthetic:2", "synthetic:1"]),
                           ("assignment_canonical_sha256", "b"*64)):
            variant = deepcopy(original); variant[key] = value
            with self.subTest(field=key), self.assertRaises(AUDIT.LargeCellAuditError):
                AUDIT.compare_singlepoints(original, variant, "cutoff_12")

    def test_nonfinite_or_malformed_forces_cannot_pass_maximum(self):
        original = singlepoint()
        for value in (True, float("nan"), float("inf"), 10**10000):
            variant = deepcopy(original); variant["forces_kj_mol_angstrom"][0][0] = value
            with self.subTest(type=type(value).__name__), self.assertRaises(AUDIT.LargeCellAuditError):
                AUDIT.compare_singlepoints(original, variant, "cutoff_12")
        variant = deepcopy(original); variant["forces_kj_mol_angstrom"][0] = [1.]
        with self.assertRaises(AUDIT.LargeCellAuditError):
            AUDIT.compare_singlepoints(original, variant, "cutoff_12")

    def test_inconsistent_total_is_rejected(self):
        original = singlepoint(); variant = deepcopy(original)
        variant["energy_components_kj_mol"]["total"] += .001
        with self.assertRaisesRegex(AUDIT.LargeCellAuditError, "total energy"):
            AUDIT.compare_singlepoints(original, variant, "cutoff_12")


class MotifReferenceTests(unittest.TestCase):
    def test_hand_calculated_single_site_simple_lattice(self):
        geometry = SimpleNamespace(coordinates_angstrom=((0., 0., 0.),))
        recipe = {"rotation_matrix": None, "repeats": [2, 2, 2], "spacing_angstrom": [2., 3., 4.],
                  "min_distance_angstrom": 1.5}
        result = AUDIT.motif_contact_reference(geometry, recipe)
        self.assertEqual(result["minimum_periodic_intermolecular_distance_angstrom"], 2.)
        self.assertEqual(result["farther_offset_distance_lower_bound_angstrom"], 4.)
        self.assertEqual(result["atom_pair_distance_evaluations"], 26)

    def test_nearby_forbidden_contact_fails_independent_reference(self):
        geometry = SimpleNamespace(coordinates_angstrom=((0., 0., 0.),))
        recipe = {"rotation_matrix": None, "repeats": [2, 2, 2], "spacing_angstrom": [1., 3., 4.],
                  "min_distance_angstrom": 1.5}
        with self.assertRaisesRegex(AUDIT.LargeCellAuditError, "forbidden contact"):
            AUDIT.motif_contact_reference(geometry, recipe)

    def test_farther_lattice_offsets_require_a_valid_bound(self):
        geometry = SimpleNamespace(coordinates_angstrom=((0., 0., 0.), (5.9, 0., 0.)))
        recipe = {"rotation_matrix": None, "repeats": [2, 2, 2], "spacing_angstrom": [3., 3., 4.],
                  "min_distance_angstrom": .01}
        with self.assertRaisesRegex(AUDIT.LargeCellAuditError, "Farther lattice offsets"):
            AUDIT.motif_contact_reference(geometry, recipe)


class PublicEvidenceTests(unittest.TestCase):
    def test_summary_preserves_scientific_aggregates_and_hashes_omissions(self):
        raw = {"final_state": {"atoms": [{"id": 1}], "density_g_cm3": .7},
               "energy_components_kj_mol": {"total": 1.}, "input_script": "trusted input",
               "engine_run": {"executable": "/private/lmp", "run_directory": "/private/run", "input_sha256": "a"*64}}
        before = deepcopy(raw); result = AUDIT.public_summary(raw)
        self.assertEqual(raw, before)
        self.assertNotIn("atoms", result["final_state"])
        self.assertEqual(result["final_state"]["density_g_cm3"], .7)
        self.assertEqual(result["engine_run"]["input_sha256"], "a"*64)
        omissions = result["public_summary_policy"]["omitted_fields"]
        self.assertEqual(len(omissions), 2)
        self.assertEqual(omissions[0]["canonical_sha256"], AUDIT.digest(AUDIT.canonical_bytes(before["final_state"]["atoms"])))


if __name__ == "__main__":
    unittest.main()
