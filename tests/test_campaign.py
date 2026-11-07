"""Campaign declarations remain plans and expose infeasible bulk-state requests."""

import copy
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from polychemmd.protocols.campaign import DensityCampaignSpec, plan_density_campaign


ROOT = Path(__file__).resolve().parents[1]


def recipe():
    return json.loads((ROOT / "examples/density_campaign.json").read_text())


def large_recipe():
    return json.loads((ROOT / "examples/density_campaign_large.json").read_text())


def with_chain_counts(data, counts):
    data["carbon_counts"] = counts
    for index, rep in enumerate(data["replicates"]):
        base = 10000 + index * 10000
        rep["conformer_seeds"] = [base + i for i in range(len(counts))]
        rep["packing_seed"] = base + 5000
        rep["velocity_seed"] = base + 6000
    return data


class DensityCampaignTests(unittest.TestCase):
    def plan(self, data=None):
        return plan_density_campaign(DensityCampaignSpec.from_mapping(recipe() if data is None else data))

    def test_plan_distinguishes_three_cells_from_nine_paired_runs(self):
        result = self.plan()
        self.assertEqual(result["sampling_accounting"]["independent_cell_recipes_requested"], 3)
        self.assertEqual(len(result["work_requests"]), 9)
        self.assertEqual(result["sampling_accounting"]["independent_cells_executed"], 0)
        self.assertFalse(result["sampling_accounting"]["continued_segments_are_new_independent_cells"])
        self.assertFalse(result["sampling_accounting"]["seed_uniqueness_proves_statistical_independence"])
        self.assertEqual(result["resource_accounting"]["engine_runs_requested"], 45)
        self.assertEqual(result["resource_accounting"]["integration_steps_requested"], 450000)
        self.assertEqual(result["sampling_accounting"]["production_duration_ps_per_setting_per_cell"], 10)

    def test_current_small_cell_cannot_reach_declared_candidate_bulk_density(self):
        result = self.plan()
        gate = result["state_capacity_gate"]
        # Independent mass/volume calculation for four C10H22, not an observed density.
        mass_per_cell_g = 4 * (10 * 12.011 + 22 * 1.008) / 6.02214076e23
        min_volume_cm3 = (2 * 11) ** 3 * 1e-24
        expected = mass_per_cell_g / min_volume_cm3
        self.assertAlmostEqual(gate["strict_density_upper_bound_g_cm3"], expected, places=14)
        self.assertLess(expected, .1)
        self.assertEqual(gate["status"], "incompatible_target_range_with_current_cell_and_cutoff_bound")
        self.assertFalse(result["execution_supported"])
        self.assertFalse(gate["geometric_bound_is_sufficient_for_liquid_preparation"])

    def test_missing_target_density_blocks_state_qualification(self):
        data = recipe(); data["target_density_range_g_cm3"] = None
        self.assertEqual(self.plan(data)["state_capacity_gate"]["status"], "blocked_target_density_hypothesis_not_declared")

    def test_partial_density_range_is_flagged(self):
        data = recipe(); data["target_density_range_g_cm3"] = [.05, .2]
        self.assertEqual(self.plan(data)["state_capacity_gate"]["status"], "target_range_partly_incompatible_with_cell_and_cutoff_bound")

    def test_geometric_feasibility_cannot_qualify_liquid_or_sampling(self):
        data = recipe(); data["target_density_range_g_cm3"] = [.05, .08]
        result = self.plan(data)
        self.assertIn("box_cutoff_bound_only", result["state_capacity_gate"]["status"])
        self.assertEqual(result["target_state"]["phase_observed"], "not_assessed")
        self.assertEqual(result["qualification_gates"]["status"], "not_evaluated")
        self.assertFalse(result["production_ready"])
        self.assertEqual(result["density_prediction_validation"], "not_performed")
        self.assertFalse(result["simulation_performed"])
        self.assertFalse(result["external_comparison"]["comparison_performed"])

    def test_tiny_run_budget_is_explicitly_exceeded(self):
        data = recipe(); data["resources"]["max_engine_runs"] = 44
        result = self.plan(data)
        self.assertEqual(result["resource_accounting"]["run_count_budget_status"], "exceeds_declared_run_count")
        self.assertFalse(result["execution_supported"])
        self.assertIsNone(result["resource_accounting"]["projection"])

    def test_measured_calibration_is_partial_scheduling_projection_only(self):
        data = recipe()
        data["calibration"] = {"hardware_label": data["resources"]["hardware_label"],
                               "engine_executable_sha256": "ab" * 32, "atom_count": 128,
                               "integration_steps": 1000, "wall_seconds": 2.,
                               "cutoff_angstrom": 10., "reciprocal_nmax": 16}
        result = self.plan(data)
        projection = result["resource_accounting"]["projection"]
        self.assertEqual(projection["baseline_integration_seconds_estimated"], 300.)
        self.assertFalse(projection["measurement_authenticated_by_planner"])
        self.assertFalse(projection["whole_campaign_budget_verified"])
        self.assertIn("cutoff/Ewald variants", projection["excludes"])

    def test_calibration_must_match_hardware_size_and_baseline(self):
        data = recipe()
        data["calibration"] = {"hardware_label": data["resources"]["hardware_label"],
                               "engine_executable_sha256": "ab" * 32, "atom_count": 128,
                               "integration_steps": 1000, "wall_seconds": 2.,
                               "cutoff_angstrom": 10., "reciprocal_nmax": 16}
        for key, wrong in (("hardware_label", "different_host"), ("atom_count", 256),
                           ("cutoff_angstrom", 8.), ("reciprocal_nmax", 15),
                           ("engine_executable_sha256", "unidentified"), ("wall_seconds", 0.)):
            trial = copy.deepcopy(data); trial["calibration"][key] = wrong
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.plan(trial)

    def test_exact_roundtrip_and_immutable_snapshot(self):
        data = recipe(); spec = DensityCampaignSpec.from_mapping(data)
        self.assertEqual(spec.to_dict(), data)
        data["carbon_counts"][0] = 50
        self.assertEqual(spec.to_dict()["carbon_counts"][0], 10)
        copy_data = spec.to_dict(); copy_data["replicates"][0]["cell_id"] = "modified"
        self.assertEqual(spec.to_dict()["replicates"][0]["cell_id"], "decane_cell_01")
        self.assertEqual(spec, DensityCampaignSpec(spec.payload_json))
        self.assertEqual(self.plan()["campaign_canonical_sha256"], plan_density_campaign(spec)["campaign_canonical_sha256"])

    def test_planning_has_no_engine_or_file_side_effects(self):
        with tempfile.TemporaryDirectory() as temporary:
            before = list(Path(temporary).iterdir())
            with patch("subprocess.run", side_effect=AssertionError("must not launch")), patch("pathlib.Path.write_text", side_effect=AssertionError("must not write")):
                result = self.plan()
            self.assertFalse(result["simulation_performed"])
            self.assertEqual(list(Path(temporary).iterdir()), before)

    def test_distinct_primary_preparations_are_required(self):
        for edit in (lambda data: data["replicates"].pop(),
                     lambda data: data["replicates"][1].update(cell_id="decane_cell_01"),
                     lambda data: data["replicates"][1].update(packing_seed=data["replicates"][0]["packing_seed"]),
                     lambda data: data["replicates"][0].update(conformer_seeds=[1])):
            data = recipe(); edit(data)
            with self.assertRaises(ValueError):
                self.plan(data)

    def test_unknown_chemistry_phase_and_fields_fail(self):
        trials = []
        for key, value in (("chemistry", "polyurethane"), ("phase_target", "liquid_validated"),
                           ("schema_version", 1), ("unsupported", 0)):
            data = recipe(); data[key] = value; trials.append(data)
        data = recipe(); data["resources"]["gpu"] = True; trials.append(data)
        for data in trials:
            with self.assertRaises(ValueError):
                self.plan(data)

    def test_boolean_nonfinite_and_invalid_seed_rejected(self):
        for key in ("temperature_k", "pressure_bar"):
            for value in (True, float("nan"), float("inf"), 0, -1, "300", 10**1000):
                data = recipe(); data[key] = value
                with self.subTest(key=key, value=str(value)[:20]), self.assertRaises(ValueError):
                    self.plan(data)
        data = recipe(); data["replicates"][0]["velocity_seed"] = True
        with self.assertRaises(ValueError): self.plan(data)
        data = recipe(); data["carbon_counts"][0] = True
        with self.assertRaises(ValueError): self.plan(data)

    def test_over_capacity_or_initial_box_cutoff_conflict_rejected(self):
        data = recipe(); data["carbon_counts"] = [50] * 4
        with self.assertRaises(ValueError): self.plan(data)
        data = recipe(); data["initialization"]["density_g_cm3"] = .7
        with self.assertRaisesRegex(ValueError, "cubic box"):
            self.plan(data)

    def test_single_variable_numerical_sensitivity_required(self):
        for setting, field, value in (("cutoff_variant", "ewald_alpha_per_angstrom", .4),
                                      ("cutoff_variant", "cutoff_angstrom", 9.),
                                      ("ewald_variant", "cutoff_angstrom", 11.),
                                      ("ewald_variant", "reciprocal_nmax", 15)):
            data = recipe(); data["numerics"][setting][field] = value
            with self.subTest(setting=setting, field=field), self.assertRaises(ValueError):
                self.plan(data)

    def test_bounded_chunks_and_declared_duration_gate(self):
        for field, value in (("segment_steps", 10001), ("production_steps", 40001),
                             ("thermo_stride_steps", 3), ("timestep_fs", 2.), ("npt_steps", 6000)):
            data = recipe(); data["integration"][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                self.plan(data)
        data = recipe(); data["sampling_gates"]["min_production_duration_ps"] = 11
        with self.assertRaises(ValueError): self.plan(data)
        data = recipe(); data["sampling_gates"]["finite_size_assessment_required"] = False
        with self.assertRaises(ValueError): self.plan(data)

    def test_malformed_density_hypothesis_fails(self):
        for value in ([.8, .65], [.7, .7], [.5], [False, .8], "unknown"):
            data = recipe(); data["target_density_range_g_cm3"] = value
            with self.subTest(value=value), self.assertRaises(ValueError):
                self.plan(data)

    def test_extreme_finite_inputs_cannot_emit_infinity_or_zero_duration(self):
        data = recipe(); data["initialization"]["density_g_cm3"] = 1e-308
        with self.assertRaisesRegex(ValueError, "volume"):
            self.plan(data)
        data = recipe()
        for setting in data["numerics"].values():
            setting["cutoff_angstrom"] = 1e-308
        data["numerics"]["cutoff_variant"]["cutoff_angstrom"] = 2e-308
        with self.assertRaisesRegex(ValueError, "density upper bound"):
            self.plan(data)
        data = recipe(); data["integration"]["timestep_fs"] = 5e-324
        data["integration"]["nvt_steps"] = 1
        data["sampling_gates"]["min_production_duration_ps"] = 5e-324
        with self.assertRaisesRegex(ValueError, "duration"):
            self.plan(data)
        data = recipe()
        data["calibration"] = {"hardware_label": data["resources"]["hardware_label"],
                               "engine_executable_sha256": "ab" * 32, "atom_count": 128,
                               "integration_steps": 1000, "wall_seconds": 5e-324,
                               "cutoff_angstrom": 10., "reciprocal_nmax": 16}
        with self.assertRaisesRegex(ValueError, "seconds per step"):
            self.plan(data)

    def test_legacy_plan_matches_preserved_v06_evidence_exactly(self):
        archived = json.loads((ROOT / "evidence/density_campaign_plan_v06.json").read_text())["payload"]
        self.assertEqual(self.plan(), archived)
        self.assertNotIn("assignment_representation", DensityCampaignSpec.from_mapping(recipe()).to_dict())
        self.assertNotIn("assignment_capacity", self.plan())

    def test_large_compact_recipe_is_box_admissible_but_has_no_executed_cells(self):
        result = self.plan(large_recipe())
        self.assertEqual(result["target_state"]["atom_count_per_cell"], 2304)
        self.assertEqual(result["target_state"]["molecule_count_per_cell"], 72)
        self.assertEqual(result["assignment_capacity"]["atom_capacity"], 8192)
        self.assertEqual(result["assignment_capacity"]["max_atoms_per_chain"], 512)
        self.assertEqual(result["assignment_capacity"]["intramolecular_pairs_requested"], 35712)
        self.assertEqual(result["assignment_capacity"]["engine_nonbonded_rules_representation"], "sparse_topological_exceptions")
        mass_g = 72 * (10 * 12.011 + 22 * 1.008) / 6.02214076e23
        self.assertAlmostEqual(result["state_capacity_gate"]["strict_density_upper_bound_g_cm3"], mass_g / ((24 ** 3) * 1e-24), places=13)
        self.assertGreater(result["state_capacity_gate"]["strict_density_upper_bound_g_cm3"], .8)
        self.assertGreater(result["initial_state_request"]["cubic_box_length_angstrom"], 24.)
        self.assertEqual(result["initial_state_request"]["density_role"], "construction_input_not_prediction")
        self.assertIn("box_cutoff_bound_only", result["state_capacity_gate"]["status"])
        self.assertEqual(result["sampling_accounting"]["independent_cell_recipes_requested"], 3)
        self.assertEqual(result["sampling_accounting"]["independent_cells_executed"], 0)
        self.assertFalse(result["sampling_accounting"]["tiling_one_cell_creates_independent_preparations"])
        self.assertFalse(result["execution_supported"])
        self.assertFalse(result["production_ready"])
        self.assertFalse(result["simulation_performed"])
        self.assertFalse(result["state_capacity_gate"]["geometric_bound_is_sufficient_for_liquid_preparation"])

    def test_compact_opt_in_is_explicit_and_rejects_unknown_representation(self):
        data = large_recipe(); data.pop("assignment_representation")
        with self.assertRaises(ValueError): self.plan(data)
        for representation in (None, True, [], "dense_v3", "sparse_topological_exceptions"):
            data = large_recipe(); data["assignment_representation"] = representation
            with self.subTest(representation=representation), self.assertRaises(ValueError):
                self.plan(data)
        data = recipe(); data["assignment_representation"] = "explicit_all_pairs_v1"
        self.assertEqual(self.plan(data)["assignment_capacity"]["atom_capacity"], 512)

    def test_compact_global_atom_and_per_chain_limits(self):
        data = with_chain_counts(large_recipe(), [2] * 1024)
        result = self.plan(data)
        self.assertEqual(result["target_state"]["atom_count_per_cell"], 8192)
        data = with_chain_counts(large_recipe(), [10] * 257)
        with self.assertRaisesRegex(ValueError, "8192-atom"):
            self.plan(data)
        data = with_chain_counts(large_recipe(), [171])
        with self.assertRaises(ValueError): self.plan(data)

    def test_compact_intramolecular_audit_budget_is_not_global_pair_count(self):
        # 15 C170H342 chains fit the 512-atom-per-chain and aggregate audit budget;
        # 16 fit 8192 atoms but exceed two million intramolecular audit pairs.
        accepted = self.plan(with_chain_counts(large_recipe(), [170] * 15))
        self.assertEqual(accepted["target_state"]["atom_count_per_cell"], 7680)
        self.assertEqual(accepted["assignment_capacity"]["intramolecular_pairs_requested"], 1962240)
        with self.assertRaisesRegex(ValueError, "intramolecular-pair audit budget"):
            self.plan(with_chain_counts(large_recipe(), [170] * 16))

    def test_compact_calibration_uses_selected_size_without_certifying_budget(self):
        data = large_recipe()
        data["calibration"] = {"hardware_label": data["resources"]["hardware_label"],
                               "engine_executable_sha256": "ab" * 32, "atom_count": 2304,
                               "integration_steps": 1000, "wall_seconds": 2.,
                               "cutoff_angstrom": 10., "reciprocal_nmax": 16}
        result = self.plan(data)
        self.assertEqual(result["resource_accounting"]["projection"]["baseline_integration_seconds_estimated"], 300.)
        self.assertFalse(result["resource_accounting"]["projection"]["whole_campaign_budget_verified"])
        data["calibration"]["atom_count"] = 128
        with self.assertRaisesRegex(ValueError, "atom count"):
            self.plan(data)
        data["calibration"]["atom_count"] = 8193
        with self.assertRaises(ValueError): self.plan(data)

    def test_dense_construction_strategy_is_named_without_phase_claim(self):
        data = large_recipe()
        self.assertEqual(data["initialization"]["strategy"], "seeded_rigid_then_npt")
        self.assertEqual(self.plan(data)["target_state"]["phase_observed"], "not_assessed")
        legacy = recipe(); legacy["initialization"]["strategy"] = "seeded_rigid_then_npt"
        with self.assertRaises(ValueError): self.plan(legacy)
        data["initialization"]["strategy"] = []
        with self.assertRaises(ValueError): self.plan(data)

    def test_wrong_api_types_fail(self):
        for data in (None, [], "{}", {}):
            with self.assertRaises(ValueError): DensityCampaignSpec.from_mapping(data)
        with self.assertRaises(ValueError): DensityCampaignSpec("not JSON")
        with self.assertRaises(ValueError): plan_density_campaign(recipe())


if __name__ == "__main__":
    unittest.main()
