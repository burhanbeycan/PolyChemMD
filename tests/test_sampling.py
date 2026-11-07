"""Independent scalar examples and labeled synthetic sampling review cases."""

import copy
import math
import random
import unittest

from polychemmd.analysis.sampling import (
    DensitySamplingSpec, _correlation, _inefficiency,
    analyze_density_sampling, summarize_density_replicates,
)


def synthetic_rows(n=2048, seed=901):
    # Controlled independent random numbers for tests, never material MD data.
    rng = random.Random(seed)
    return [dict(step=i*10,temperature_k=300+rng.uniform(-.5,.5),
                 pressure_bar=1+rng.uniform(-.5,.5),
                 density_g_cm3=.8+rng.uniform(-.005,.005)) for i in range(n)]


def test_spec(**kwargs):
    return DensitySamplingSpec(timestep_fs=1,block_size=64,min_duration_ps=10,**kwargs)


def replicate(identity, analysis, status="completed"):
    return dict(independent_cell_id=identity,state_id="n-decane-300K-1bar",
                protocol_id="declared-protocol-A",status=status,sampling_analysis=analysis)


class SamplingSpecTests(unittest.TestCase):
    def test_defaults_and_round_trip(self):
        spec = DensitySamplingSpec()
        self.assertEqual(DensitySamplingSpec.from_mapping(spec.to_dict()),spec)
        self.assertEqual(spec.min_duration_ps,100)

    def test_unknown_fields(self):
        with self.assertRaises(ValueError):
            DensitySamplingSpec.from_mapping({"min_duration":100})

    def test_boolean_nonfinite_and_bounds(self):
        for name,value in (("timestep_fs",True),("discard_steps",True),
                           ("min_duration_ps",float("nan")),("target_pressure_bar",float("inf")),
                           ("block_size",0),("min_blocks",15),("max_autocorrelation_lag",3),
                           ("max_abs_block_lag_one",1),("max_coarsened_se_ratio",1)):
            with self.subTest(name=name,value=value),self.assertRaises(ValueError):
                DensitySamplingSpec.from_mapping({name:value})


class SamplingReviewTests(unittest.TestCase):
    def test_hand_computed_correlations_and_monotone_pairs(self):
        # x=(0,1,2,3); mean=1.5; sum centered squares=5;
        # rho1=1.25/5, rho2=-1.5/5, rho3=-2.25/5.
        values=[0,1,2,3]
        self.assertAlmostEqual(_correlation(values,1),.25)
        self.assertAlmostEqual(_correlation(values,2),-.3)
        self.assertAlmostEqual(_correlation(values,3),-.45)
        result=_inefficiency(values,4)
        self.assertTrue(result["window_closed"])
        self.assertAlmostEqual(result["statistical_inefficiency"],1.5)
        self.assertAlmostEqual(result["effective_samples"],8/3)

    def test_independent_synthetic_pass_is_diagnostics_only(self):
        report=analyze_density_sampling(synthetic_rows(),test_spec())
        self.assertEqual(report["sampling_status"],"diagnostics_passed")
        self.assertEqual(report["flags"],[])
        self.assertEqual(report["equilibrium_status"],"not_assessed")
        self.assertFalse(report["production_ready"])
        self.assertEqual(report["independent_cell_count"],1)
        self.assertFalse(report["simulation_executed_by_analysis"])
        self.assertGreater(report["conditional_autocorrelation_standard_error_g_cm3"],0)

    def test_short_actual_style_series_cannot_pass(self):
        report=analyze_density_sampling(synthetic_rows(5),DensitySamplingSpec())
        self.assertEqual(report["sampling_status"],"insufficient_sampling")
        self.assertFalse(report["gates"]["physical_duration"])
        self.assertNotIn("conditional_block_standard_error_g_cm3",report)

    def test_predeclared_duration_cannot_be_replaced_by_stationarity(self):
        report=analyze_density_sampling(synthetic_rows(),DensitySamplingSpec(timestep_fs=1,block_size=64))
        self.assertEqual(report["sampling_status"],"insufficient_sampling")
        self.assertFalse(report["gates"]["physical_duration"])
        self.assertTrue(report["gates"]["nonconstant_density_and_blocks"])

    def test_duration_gate_uses_analyzed_frames_not_omitted_tail(self):
        rows=synthetic_rows(2100)
        spec=DensitySamplingSpec(timestep_fs=1,block_size=64,min_duration_ps=20.8)
        report=analyze_density_sampling(rows,spec)
        self.assertAlmostEqual(report["retained_duration_ps"],20.99)
        self.assertAlmostEqual(report["used_duration_ps"],20.47)
        self.assertFalse(report["gates"]["physical_duration"])
        self.assertEqual(report["n_omitted_trailing_frames"],52)

    def test_constant_density_is_inconclusive(self):
        rows=synthetic_rows()
        for row in rows:
            row["density_g_cm3"] = .8
        report=analyze_density_sampling(rows,test_spec())
        self.assertEqual(report["sampling_status"],"inconclusive")
        self.assertIsNone(report["conditional_normal_half_width_g_cm3"])
        self.assertIsNone(report["autocorrelation"]["effective_samples"])

    def test_identical_block_means_do_not_establish_sampling(self):
        rows=synthetic_rows()
        for i,row in enumerate(rows):
            row["density_g_cm3"] = .79 if i%2 else .81
        report=analyze_density_sampling(rows,test_spec())
        self.assertEqual(report["sampling_status"],"inconclusive")
        self.assertFalse(report["gates"]["nonconstant_density_and_blocks"])

    def test_correlated_synthetic_ar1_requires_longer_blocks(self):
        rows=synthetic_rows(4096)
        rng=random.Random(75)
        value=0
        for row in rows:
            value=.95*value+rng.gauss(0,.001)
            row["density_g_cm3"] = .8+value
        report=analyze_density_sampling(rows,test_spec())
        self.assertGreater(report["autocorrelation"]["statistical_inefficiency"],10)
        self.assertLess(report["autocorrelation"]["effective_samples"],len(rows)/10)
        self.assertFalse(report["gates"]["block_correlation_time_separation"])
        self.assertEqual(report["sampling_status"],"requires_review")

    def test_unclosed_correlation_window_is_failure(self):
        report=analyze_density_sampling(synthetic_rows(),test_spec(max_autocorrelation_lag=4))
        # Exact hand monotone trend guarantees positive low-lag pairs.
        trend=_inefficiency(list(range(200)),4)
        self.assertFalse(trend["window_closed"])
        self.assertEqual(trend["status"],"positive_window_not_closed")

    def test_hand_computed_block_standard_error(self):
        rows=[dict(step=i,temperature_k=300,pressure_bar=1,density_g_cm3=.8+.01*(i//2)) for i in range(32)]
        report=analyze_density_sampling(rows,DensitySamplingSpec(timestep_fs=1,block_size=2,min_duration_ps=.001))
        # Sixteen block means .80,.81,...,.95; sample variance=.0001*16*17/12.
        self.assertAlmostEqual(report["conditional_block_standard_error_g_cm3"],math.sqrt(.0001*16*17/12/16))
        self.assertFalse(report["gates"]["density_split_half_drift"])
        self.assertFalse(report["gates"]["density_four_window_drift"])

    def test_pressure_checks_mean_instead_of_individual_excursions(self):
        rows=synthetic_rows()
        for i,row in enumerate(rows):
            row["pressure_bar"] = 1001 if i%2 else -999
        report=analyze_density_sampling(rows,test_spec(pressure_mean_tolerance_bar=1))
        self.assertEqual(report["mean_pressure_bar"],1)
        self.assertTrue(report["gates"]["pressure_sample_mean"])

    def test_target_mean_state_mismatch_requires_review(self):
        rows=synthetic_rows()
        for row in rows:
            row["temperature_k"]+=50
            row["pressure_bar"]+=500
        report=analyze_density_sampling(rows,test_spec())
        self.assertFalse(report["gates"]["temperature_sample_mean"])
        self.assertFalse(report["gates"]["pressure_sample_mean"])
        self.assertEqual(report["sampling_status"],"requires_review")

    def test_discard_relative_to_first_step(self):
        rows=synthetic_rows()
        for row in rows:
            row["step"]+=2000
        report=analyze_density_sampling(rows,test_spec(discard_steps=100))
        self.assertEqual(report["n_discarded_frames"],10)
        self.assertEqual(report["n_retained_frames"],2038)

    def test_all_discarded_is_insufficient_not_synthetic_success(self):
        report=analyze_density_sampling(synthetic_rows(5),test_spec(discard_steps=10000))
        self.assertEqual(report["status"],"insufficient_sampling")
        self.assertEqual(report["n_used_frames"],0)

    def test_invalid_rows_rejected(self):
        valid=synthetic_rows(5)
        cases=[[],"text",[dict(valid[0],step=True)],[dict(valid[0],density_g_cm3=0)],
               [dict(valid[0],temperature_k=float("nan"))],[dict(valid[0],pressure_bar=True)],
               [dict(valid[0],extra=1)],valid[:2]+[dict(valid[2],step=21)],
               valid[:2]+[valid[1]]]
        for rows in cases:
            with self.subTest(rows=rows),self.assertRaises(ValueError):
                analyze_density_sampling(rows,test_spec())

    def test_autocorrelation_work_budget_is_explicit(self):
        with self.assertRaisesRegex(ValueError,"work budget"):
            analyze_density_sampling(synthetic_rows(4096),test_spec(max_autocorrelation_lag=4096))


class ReplicateReviewTests(unittest.TestCase):
    def setUp(self):
        self.analysis=analyze_density_sampling(synthetic_rows(),test_spec())
        self.records=[]
        for identity,mean in (("cell-A",.79),("cell-B",.8),("cell-C",.81)):
            analysis=copy.deepcopy(self.analysis)
            analysis["mean_density_g_cm3"]=mean
            self.records.append(replicate(identity,analysis))
        self.ids=["cell-A","cell-B","cell-C"]

    def test_independent_cell_standard_error_is_not_frame_uncertainty(self):
        report=summarize_density_replicates(self.records,expected_cell_ids=self.ids)
        self.assertEqual(report["status"],"replicate_diagnostics_passed")
        self.assertEqual(report["n_declared_independent_cells"],3)
        self.assertAlmostEqual(report["mean_density_g_cm3"],.8)
        self.assertAlmostEqual(report["between_cell_standard_deviation_g_cm3"],.01)
        self.assertAlmostEqual(report["conditional_between_cell_standard_error_g_cm3"],.01/math.sqrt(3))
        self.assertFalse(report["frames_are_independent_cells"])
        self.assertEqual(report["independence_status"],"declared_not_verified")

    def test_failed_cell_retained_prevents_reported_aggregate(self):
        self.records[1]=replicate("cell-B",None,status="failed")
        report=summarize_density_replicates(self.records,expected_cell_ids=self.ids)
        self.assertEqual(report["n_diagnostics_passed_cells"],2)
        self.assertIsNone(report["mean_density_g_cm3"])
        self.assertEqual(report["cell_outcomes"][1]["execution_status"],"failed")

    def test_missing_failed_cell_is_rejected(self):
        with self.assertRaisesRegex(ValueError,"every declared"):
            summarize_density_replicates(self.records[:2],expected_cell_ids=self.ids)

    def test_duplicate_continuations_cannot_be_replicates(self):
        self.records[1]["independent_cell_id"]="cell-A"
        with self.assertRaisesRegex(ValueError,"repeated frames"):
            summarize_density_replicates(self.records,expected_cell_ids=self.ids)

    def test_mixed_states_protocols_or_specs_are_rejected(self):
        for name,value in (("state_id","other-state"),("protocol_id","other-protocol")):
            records=copy.deepcopy(self.records)
            records[1][name]=value
            with self.subTest(name=name),self.assertRaises(ValueError):
                summarize_density_replicates(records,expected_cell_ids=self.ids)
        self.records[1]["sampling_analysis"]["spec"]["target_temperature_k"]=301
        with self.assertRaises(ValueError):
            summarize_density_replicates(self.records,expected_cell_ids=self.ids)

    def test_identical_cell_means_do_not_imply_zero_uncertainty(self):
        for record in self.records:
            record["sampling_analysis"]["mean_density_g_cm3"] = .8
        report=summarize_density_replicates(self.records,expected_cell_ids=self.ids)
        self.assertEqual(report["status"],"incomplete_or_inconclusive")
        self.assertIsNone(report["conditional_between_cell_standard_error_g_cm3"])

    def test_inconsistent_passed_summaries_rejected(self):
        for field,value in (("production_ready",True),("equilibrium_status","equilibrated"),
                            ("status","requires_review"),("flags",["physical_duration"]),
                            ("used_duration_ps",.001),("n_used_frames",1),
                            ("mean_temperature_k",900)):
            records=copy.deepcopy(self.records)
            records[0]["sampling_analysis"][field]=value
            with self.subTest(field=field),self.assertRaises(ValueError):
                summarize_density_replicates(records,expected_cell_ids=self.ids)
        self.records[0]["sampling_analysis"]["gates"]["effective_samples"]=False
        with self.assertRaisesRegex(ValueError,"every predeclared"):
            summarize_density_replicates(self.records,expected_cell_ids=self.ids)

    def test_wrong_expected_identity_rejected(self):
        with self.assertRaisesRegex(ValueError,"identities"):
            summarize_density_replicates(self.records,expected_cell_ids=["cell-A","cell-B","cell-D"])

    def test_only_two_cells_cannot_close_replicate_review(self):
        report=summarize_density_replicates(self.records[:2],expected_cell_ids=self.ids[:2])
        self.assertEqual(report["status"],"incomplete_or_inconclusive")
        self.assertIsNone(report["mean_density_g_cm3"])


if __name__=="__main__":
    unittest.main()
