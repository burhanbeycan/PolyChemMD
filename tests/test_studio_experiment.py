import csv
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from polychemmd.studio.experiment import atlas, run_experiment


class MeasurementEvidenceTests(unittest.TestCase):
    def test_literature_observables_are_traceable(self):
        data = atlas()
        ids = {r["id"] for r in data["records"]}
        self.assertEqual(len(ids), 39)
        cited = {s["doi"] for s in data["sources"]}
        self.assertTrue(all(r["source_doi"] in cited for r in data["records"]))
        pce = next(r for r in data["records"] if r["id"] == "park-PCE-NaCl")
        self.assertAlmostEqual(pce["value"], 1.3322)
        self.assertTrue(
            pce["reading_bounds"][0] < pce["value"] < pce["reading_bounds"][1]
        )

    def test_published_wet_minus_dry_effect(self):
        with TemporaryDirectory() as directory:
            out = Path(directory)
            r = run_experiment({"assay": "xrd", "reference": "xing-pei"}, out)
            self.assertEqual(r["summary"]["n_pairs"], 1)
            self.assertIsNone(r["summary"]["sample_sd"])
            self.assertAlmostEqual(r["rows"][0]["dry_effect_nm"], 0.060)
            self.assertAlmostEqual(r["rows"][0]["wet_effect_nm"], -0.112)
            self.assertAlmostEqual(r["summary"]["mean"], -0.172)
            self.assertEqual(json.loads((out / "experiment.json").read_text()), r)
            with (out / "experiment.csv").open() as handle:
                self.assertEqual(len(list(csv.DictReader(handle))), 1)
            self.assertEqual(r["evidence"], "published secondary analysis")

    def test_optical_with_matched_and_unmatched_batches(self):
        rows = [
            dict(
                batch="A",
                medium="salt",
                polymer="PCE",
                A0=1.2,
                A24=1.1,
                B0=0.2,
                B24=0.3,
            ),
            dict(
                batch="A",
                medium="salt",
                polymer="PEO",
                A0=0.9,
                A24=0.4,
                B0=0.1,
                B24=0.2,
            ),
            dict(
                batch="B",
                medium="water",
                polymer="PCE",
                A0=0.8,
                A24=0.2,
                B0=0.2,
                B24=0.3,
            ),
        ]
        with TemporaryDirectory() as directory:
            r = run_experiment({"assay": "optical", "rows": rows}, Path(directory))
        self.assertAlmostEqual(r["rows"][0]["R24"], 0.8)
        self.assertAlmostEqual(r["rows"][1]["R24"], 0.25)
        self.assertLess(r["rows"][2]["R24"], 0)
        self.assertAlmostEqual(r["summary"]["by_medium"]["salt"]["mean"], 0.55)
        self.assertEqual(
            r["summary"]["unpaired_batch_media"], [{"batch": "B", "medium": "water"}]
        )

    def test_invalid_measurements_do_not_generate_results(self):
        valid = dict(
            batch="A", medium="salt", polymer="PCE", A0=1, A24=0.8, B0=0.1, B24=0.1
        )
        with TemporaryDirectory() as directory:
            out = Path(directory)
            cases = [
                [valid, valid],
                [dict(valid, A0=True)],
                [dict(valid, B24=float("inf"))],
                [dict(valid, A0=0.05)],
                [dict(valid, batch="")],
            ]
            for rows in cases:
                with self.subTest(rows=rows), self.assertRaises(ValueError):
                    run_experiment(dict(assay="optical", rows=rows), out)
            with self.assertRaises(ValueError):
                run_experiment(dict(assay="xrd", reference="xing-pei", rows=[]), out)
            self.assertEqual(list(out.iterdir()), [])
