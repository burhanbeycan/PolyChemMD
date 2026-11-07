"""Synthetic audit failure cases; these records are not molecular trajectories."""

from copy import deepcopy
import importlib.util
from pathlib import Path
import unittest


SCRIPT = Path(__file__).resolve().parents[1]/"scripts"/"verify_restart_v06.py"
SPEC = importlib.util.spec_from_file_location("polychemmd_restart_auditor_test", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def endpoint():
    return {"step": 240, "atom_count": 1,
            "box": {"origin_angstrom": [-.1, -.1, -.1], "lengths_angstrom": [24., 24., 24.]},
            "potential_energy_kj_mol": -3., "total_energy_kj_mol": 2.,
            "atoms": [{"engine_id": 1, "cell_atom_id": "synthetic/molecule:0/atom:1",
                       "molecule_index": 0, "charge_e": 0., "image_flags": [1, 0, 0],
                       "coordinates_angstrom": [25., 2., 3.],
                       "wrapped_coordinates_angstrom": [1., 2., 3.],
                       "velocities_angstrom_fs": [.01, -.01, .001],
                       "forces_kj_mol_angstrom": [1., 2., 3.]}]}


class RestartAuditTests(unittest.TestCase):
    def test_identical_synthetic_endpoint_passes_without_sampling_claim(self):
        record = AUDIT.compare_endpoints(endpoint(), endpoint())
        self.assertEqual(record["status"], "passed")
        self.assertEqual(set(record["maximum_absolute_errors"].values()), {0.})
        self.assertIn("not independent sampling", record["interpretation"])

    def test_periodic_equivalence_cannot_hide_unwrapped_discontinuity(self):
        moved = endpoint()
        moved["atoms"][0]["coordinates_angstrom"][0] += 24.
        with self.assertRaisesRegex(AUDIT.RestartAuditError, "position_component"):
            AUDIT.compare_endpoints(endpoint(), moved)

    def test_velocity_is_checked_even_when_positions_match(self):
        changed = endpoint()
        changed["atoms"][0]["velocities_angstrom_fs"][0] += 1e-5
        with self.assertRaisesRegex(AUDIT.RestartAuditError, "velocity_component"):
            AUDIT.compare_endpoints(endpoint(), changed)

    def test_identity_images_and_physical_step_cannot_be_relabelled(self):
        for key, value in (("cell_atom_id", "different"), ("molecule_index", 1),
                           ("engine_id", 2), ("charge_e", .06), ("image_flags", [0, 0, 0])):
            changed = endpoint()
            changed["atoms"][0][key] = value
            with self.subTest(key=key), self.assertRaises(AUDIT.RestartAuditError):
                AUDIT.compare_endpoints(endpoint(), changed)
        changed = endpoint()
        changed["step"] = 40
        with self.assertRaises(AUDIT.RestartAuditError):
            AUDIT.compare_endpoints(endpoint(), changed)

    def test_energy_and_box_changes_cannot_pass_coordinate_only_check(self):
        changed = endpoint()
        changed["total_energy_kj_mol"] += .01
        with self.assertRaisesRegex(AUDIT.RestartAuditError, "energy_kj_mol"):
            AUDIT.compare_endpoints(endpoint(), changed)
        changed = endpoint()
        changed["box"]["origin_angstrom"][0] += .01
        with self.assertRaisesRegex(AUDIT.RestartAuditError, "box_component"):
            AUDIT.compare_endpoints(endpoint(), changed)

    def test_nonfinite_endpoint_cannot_silently_pass_maximum(self):
        for value in (float("nan"), float("inf"), True):
            changed = endpoint()
            changed["atoms"][0]["velocities_angstrom_fs"][0] = value
            with self.subTest(value=value), self.assertRaises(AUDIT.RestartAuditError):
                AUDIT.compare_endpoints(endpoint(), changed)

    def test_public_sanitization_preserves_scientific_identity_and_input(self):
        raw = {"engine_run": {"executable": "/synthetic/private/lmp", "run_directory": "/synthetic/run",
                              "executable_sha256": "a"*64, "command": ["/synthetic/private/lmp", "-in", "input.lammps"]},
               "checkpoint": {"sha256": "b"*64}, "density_g_cm3": .123}
        before = deepcopy(raw)
        public = AUDIT.public_record(raw)
        self.assertEqual(raw, before)
        self.assertEqual(public["engine_run"]["executable"], "lmp")
        self.assertEqual(public["checkpoint"], before["checkpoint"])
        self.assertEqual(public["density_g_cm3"], .123)
        self.assertEqual(public["engine_run"]["executable_sha256"], "a"*64)

    def test_actual_thermo_conversion_path_retains_short_sampling_failure(self):
        # Synthetic numbers exercise the unit boundary; they are not MD data.
        columns = ["Step", "Temp", "Press", "Volume"]
        result = {"engine_run": {"thermo_blocks": [
                  {"columns": columns, "rows": [[200, 300., 1., 1000.]]},
                  {"columns": columns, "rows": [[200+i*10, 300., 1., 1000.+i] for i in range(5)]}]}}
        assignment = {"atoms": [{"mass_g_mol": 1.008}, {"mass_g_mol": 1.008}]}
        reviewed = AUDIT.assess_continuation_sampling(result, assignment)
        row = reviewed["canonical_thermo_rows"][0]
        self.assertAlmostEqual(row["pressure_bar"], 1.01325)
        self.assertAlmostEqual(row["density_g_cm3"], 2.016/(6.02214076e23*1000.*1e-24))
        self.assertEqual(reviewed["analysis"]["sampling_status"], "insufficient_sampling")
        self.assertEqual(reviewed["analysis"]["n_retained_frames"], 5)
        self.assertAlmostEqual(reviewed["analysis"]["retained_duration_ps"], .01)
        self.assertFalse(reviewed["analysis"]["production_ready"])


if __name__ == "__main__":
    unittest.main()
