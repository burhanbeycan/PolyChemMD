"""Traceable result envelopes; no fabricated engine or material validation."""

from datetime import datetime, timezone
import hashlib
import platform

from . import __version__


def envelope(payload, input_bytes, record_kind):
    return {
        "result_schema_version": "0.1",
        "record_kind": record_kind,
        "software": {"name": "PolyChemMD", "version": __version__, "python": platform.python_version()},
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "input_sha256": hashlib.sha256(input_bytes).hexdigest(),
        "input_hash_basis": "exact input-file bytes",
        "payload": payload,
    }
