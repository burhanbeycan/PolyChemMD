"""Attributed numerical subset; no engine or upstream workflow is imported."""

from importlib.resources import files
import json

from .contracts import ParameterSet


MODEL_ID = "OPLS_AA_GROMACS_2025_3_N_ALKANE_SUBSET"


def load_alkane_parameters():
    """Return a fresh validated copy of the frozen CT/HC numerical subset.

    Per-type sources identify the exact upstream commit and source lines.
    Numerical transcription/conversion does not establish physical accuracy.
    """
    data = json.loads(files(__package__).joinpath("alkane_parameters.json").read_text(encoding="utf-8"))
    if data["schema_version"] != "1" or data["model_id"] != MODEL_ID:
        raise ValueError("unsupported alkane parameter schema/model")
    ParameterSet.from_mapping(data["parameter_set"])
    return data
