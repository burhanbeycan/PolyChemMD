"""Synthetic auditor invariants; these fixtures are not simulated polymers."""
from copy import deepcopy
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1]/"scripts/verify_qualification_v08.py"
SPEC = importlib.util.spec_from_file_location("polychemmd_qualification_auditor_test", SCRIPT)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def synthetic_accounting_cell():
    atoms = []
    for molecule in range(147):
        for element, count, mass in (("C", 10, 12.011), ("H", 22, 1.008)):
            atoms.extend({"molecule_index": molecule, "element": element, "mass_g_mol": mass}
                         for _ in range(count))
    box = [39., 35.7, 35.7]
    return {"atoms": atoms, "molecules": [{} for _ in range(147)],
            "box": {"lengths_angstrom": box},
            "construction_density_g_cm3": math.fsum(atom["mass_g_mol"] for atom in atoms)/(.602214076*math.prod(box))}


class QualificationAccountingTests(unittest.TestCase):
    def test_predeclared_resource_accounting_contains_all_actual_requests(self):
        self.assertEqual(len(AUDIT.NUMERICAL_SPEC["cutoffs_angstrom"])+1+2, AUDIT.MAX_ENGINE_INVOCATIONS)
        self.assertEqual(AUDIT.MAX_ENGINE_INVOCATIONS, 8)
        self.assertEqual(AUDIT.NUMERICAL_SPEC["total_wall_seconds"], 120.)
        self.assertLess(AUDIT.NUMERICAL_SPEC["total_wall_seconds"], AUDIT.TOTAL_WALL_SECONDS)
        self.assertEqual(AUDIT.PREPARATION["cutoff_angstrom"], max(AUDIT.NUMERICAL_SPEC["cutoffs_angstrom"]))
        self.assertEqual(AUDIT.PREPARATION["reciprocal_nmax"], AUDIT.NUMERICAL_SPEC["refined_reciprocal_nmax"])
        self.assertEqual(AUDIT.PREPARATION["nvt_steps"]+AUDIT.PREPARATION["npt_steps"]+AUDIT.CONTINUATION_STEPS, 30)

    def test_geometric_accounting_remains_one_cell_without_phase_claim(self):
        result = AUDIT.audit_cell(synthetic_accounting_cell())
        self.assertEqual(result["atoms"], 4704)
        self.assertEqual(result["molecules"], 147)
        self.assertEqual(result["independently_prepared_cells"], 1)
        self.assertFalse(result["molecules_and_variants_are_independent_cells"])
        self.assertEqual(result["phase"], "not_assigned")
        self.assertAlmostEqual(result["construction_density_g_cm3"], .6987592396767455)

    def test_wrong_molecular_composition_and_density_are_rejected(self):
        original = synthetic_accounting_cell()
        for mutation in ("element", "density", "molecule"):
            changed = deepcopy(original)
            if mutation == "element": changed["atoms"][0]["element"] = "O"
            elif mutation == "density": changed["construction_density_g_cm3"] *= 1.01
            else: changed["atoms"][0]["molecule_index"] = 147
            with self.subTest(mutation=mutation), self.assertRaises(AUDIT.QualificationAuditError):
                AUDIT.audit_cell(changed)

    def test_cutoff_incompatible_or_wrong_box_is_not_accepted(self):
        changed = synthetic_accounting_cell(); changed["box"]["lengths_angstrom"][1] = 33.
        with self.assertRaises(AUDIT.QualificationAuditError): AUDIT.audit_cell(changed)
        with patch.dict(AUDIT.NUMERICAL_SPEC, {"cutoffs_angstrom": [10., 18.]}):
            with self.assertRaisesRegex(AUDIT.QualificationAuditError, "cutoff"):
                AUDIT.audit_cell(synthetic_accounting_cell())

    def test_published_recipe_matches_predeclared_construction(self):
        path = SCRIPT.parents[1]/AUDIT.RECIPE_PATH
        self.assertEqual(json.loads(path.read_text()), AUDIT.RECIPE)
        self.assertEqual(AUDIT.RECIPE["repeats"], [3, 7, 7])
        self.assertLess(2*max(AUDIT.NUMERICAL_SPEC["cutoffs_angstrom"]), min(
            a*b for a,b in zip(AUDIT.RECIPE["repeats"], AUDIT.RECIPE["spacing_angstrom"])))

    def test_child_phase_cannot_launch_without_its_whole_outer_allowance(self):
        self.assertEqual(AUDIT.require_numerical_phase_budget(110.), 130.)
        with self.assertRaisesRegex(AUDIT.QualificationAuditError, "Remaining outer"):
            AUDIT.require_numerical_phase_budget(110.001)
        for value in (True, -1., float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(AUDIT.QualificationAuditError):
                AUDIT.require_numerical_phase_budget(value)


class EvidenceBoundaryTests(unittest.TestCase):
    def test_short_pilot_review_cannot_be_promoted(self):
        original = {"status": "requires_preparation", "production_ready": False, "equilibrium_status": "not_assessed"}
        self.assertEqual(AUDIT.checked_preparation_review(original), original)
        for key, value in (("production_ready", True), ("equilibrium_status", "equilibrated"), ("status", "validated")):
            changed = dict(original, **{key: value})
            with self.subTest(field=key), self.assertRaises(AUDIT.QualificationAuditError):
                AUDIT.checked_preparation_review(changed)

    def test_source_change_is_detected_without_running_engine(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory); path = root/"module.py"; path.write_text("original source\n")
            record = {"module.py": {"sha256": hashlib.sha256(path.read_bytes()).hexdigest()}}
            with patch.object(AUDIT, "REPOSITORY", root):
                AUDIT.assert_sources_unchanged(record)
                path.write_text("changed source\n")
                with self.assertRaisesRegex(AUDIT.QualificationAuditError, "source changed"):
                    AUDIT.assert_sources_unchanged(record)

    def test_partial_failure_retains_actual_text_but_no_native_binary(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory); private = output/"private_engine_runs"
            run = private/"unfinished"; run.mkdir(parents=True)
            (run/"input.lammps").write_text("run 100\n")
            (run/"log.lammps").write_text("incomplete actual log\n")
            (run/"checkpoint.restart").write_bytes(b"native binary")
            result = AUDIT.retain_failed_inputs(private, output)
            self.assertEqual(set(result), {"unfinished/input.lammps", "unfinished/log.lammps"})
            self.assertEqual((output/result["unfinished/log.lammps"]["public_path"]).read_text(), "incomplete actual log\n")
            self.assertFalse(list((output/"failed_inputs").rglob("*.restart")))
            self.assertIn("not_success", result["unfinished/log.lammps"]["evidence"])

    def test_non_utf8_partial_log_does_not_hide_original_failure(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory); private = output/"private_engine_runs"
            run = private/"unfinished"; run.mkdir(parents=True)
            (run/"log.lammps").write_bytes(b"\xff\xfeinvalid log")
            result = AUDIT.retain_failed_inputs(private, output)
            item = result["unfinished/log.lammps"]
            self.assertIsNone(item["public_path"])
            self.assertEqual(item["sha256"], hashlib.sha256((run/"log.lammps").read_bytes()).hexdigest())
            self.assertIn("hash_only", item["evidence"])

    def test_public_summary_retains_static_pressure_and_hashes_array_omissions(self):
        raw = {"pressure_bar": 12., "pressure_tensor_bar": [[1., 0., 0.], [0., 2., 0.], [0., 0., 33.]],
               "forces_kj_mol_angstrom": [[1., 2., 3.]], "final_state": {"atoms": [{"id": 1}]},
               "engine_run": {"input_sha256": "a"*64, "run_directory": "/private/host"}}
        before = deepcopy(raw); result = AUDIT.public_summary(raw)
        self.assertEqual(raw, before)
        self.assertEqual(result["pressure_bar"], 12.)
        self.assertEqual(result["pressure_tensor_bar"], before["pressure_tensor_bar"])
        self.assertNotIn("forces_kj_mol_angstrom", result)
        self.assertNotIn("atoms", result["final_state"])
        omissions = result["public_summary_policy"]["omitted_fields"]
        self.assertTrue(any(item["canonical_sha256"] == AUDIT.digest(AUDIT.canonical_bytes(before["forces_kj_mol_angstrom"])) for item in omissions))

    def test_repeated_static_chain_rows_have_explicit_hash_and_keep_aggregates(self):
        raw = {"structure_snapshot": {"per_chain": [{"molecule_index": 0, "radius": 2.}],
            "nematic_order_parameter": 1., "mean_end_to_end_distance_angstrom": 8.},
            "configuration_statistics_are_equilibrium_averages": False}
        result = AUDIT.public_summary(raw)
        self.assertIn("per_chain", raw["structure_snapshot"])
        self.assertNotIn("per_chain", result["structure_snapshot"])
        self.assertEqual(result["structure_snapshot"]["nematic_order_parameter"], 1.)
        self.assertFalse(result["configuration_statistics_are_equilibrium_averages"])
        omission = next(item for item in result["public_summary_policy"]["omitted_fields"] if item["field"].endswith("per_chain"))
        self.assertEqual(omission["row_count"], 1)
        self.assertEqual(omission["canonical_sha256"], AUDIT.digest(AUDIT.canonical_bytes(raw["structure_snapshot"]["per_chain"])))


if __name__ == "__main__": unittest.main()
