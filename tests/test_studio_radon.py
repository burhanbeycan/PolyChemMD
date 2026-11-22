import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

from polychemmd.studio.radon import run_radon, status


class ExternalReferenceTests(unittest.TestCase):
    def test_unavailable_optional_interpreter_is_explicit(self):
        with TemporaryDirectory() as directory:
            bad = str(Path(directory) / "missing")
            with patch.dict(os.environ, {"POLYCHEMMD_RADON_PYTHON": bad}):
                self.assertFalse(status()["available"])
                with self.assertRaisesRegex(RuntimeError, "unavailable"):
                    run_radon({}, Path(directory))
                for config in (
                    dict(polymer="invented"),
                    dict(dp=False),
                    dict(dp=17),
                    dict(seed=-1),
                ):
                    with self.subTest(config=config), self.assertRaises(ValueError):
                        run_radon(config, Path(directory))

    @unittest.skipUnless(
        os.getenv("POLYCHEMMD_RADON_PYTHON"),
        "isolated reference environment not selected",
    )
    def test_external_peo_stereo_and_mapping(self):
        with TemporaryDirectory() as directory:
            out = Path(directory)
            cfg = dict(polymer="PEO", dp=3, seed=17)
            result = run_radon(cfg, out / "a")
            repeated = run_radon(cfg, out / "b")
            self.assertEqual(result["formula"], "C8H18O3")
            self.assertEqual(result["sdf_sha256"], repeated["sdf_sha256"])
            self.assertEqual(result["gaff2_typed"], result["n_atoms"])
            structure = json.loads((out / "a" / "structure.json").read_text())
            self.assertEqual(len(structure["atoms"]), 29)
            self.assertEqual(len(structure["bonds"]), 28)
            self.assertEqual(
                [v["atom_map"] for v in structure["atoms"]], list(range(1, 30))
            )
