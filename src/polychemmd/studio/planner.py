"""A local planner can choose only predeclared, paired physical comparisons."""

from copy import deepcopy
import http.client
import json
import re

PAIRS = {
    "termination": [
        {
            "kind": "geometry",
            "config": {"polymer": "PEO", "dp": 8, "pattern": pattern, "seed": 17},
        }
        for pattern in ("random", "patch")
    ],
    "hydrolysis": [
        {
            "kind": "geometry",
            "config": {
                "polymer": "PVA",
                "dp": 12,
                "hydrolysis": 0.75,
                "acetate_pattern": p,
                "seed": 17,
            },
        }
        for p in ("random", "clustered")
    ],
    "model": [
        {
            "kind": "md",
            "config": {
                "md_profile": p,
                "md_sequence": "clustered",
                "md_steps": 20,
                "seed": 17,
            },
        }
        for p in ("iff-opls", "iff-charmm")
    ],
}


def suggest(mode="rules", model=None, objective="termination"):
    if (
        not isinstance(objective, str)
        or objective not in PAIRS
        or mode not in ("rules", "ollama")
    ):
        raise ValueError("Unknown planner mode or objective.")
    selected = objective
    if mode == "ollama":
        if not isinstance(model, str) or not re.fullmatch(r"[\w./:-]{1,100}", model):
            raise ValueError("Choose an installed local Ollama model.")
        schema = {
            "type": "object",
            "required": ["pair"],
            "additionalProperties": False,
            "properties": {"pair": {"type": "string", "enum": list(PAIRS)}},
        }
        body = {
            "model": model,
            "stream": False,
            "format": schema,
            "options": {"temperature": 0, "num_predict": 64},
            "prompt": f"Select one matched comparison. Objective {objective}. Allowed pairs: {', '.join(PAIRS)}. Return JSON.",
        }
        conn = http.client.HTTPConnection("127.0.0.1", 11434, timeout=30)
        try:
            conn.request(
                "POST",
                "/api/generate",
                json.dumps(body),
                {"Content-Type": "application/json"},
            )
            resp = conn.getresponse()
            raw = resp.read(65537)
            if resp.status != 200 or len(raw) > 65536:
                raise ValueError(f"Local model is unavailable (HTTP {resp.status}).")
            envelope = json.loads(raw)
            if not isinstance(envelope, dict) or envelope.get("done") is not True:
                raise ValueError("The local generation did not complete.")
            if not isinstance(envelope.get("response"), str):
                raise ValueError("The local generation returned no text response.")
            reply = json.loads(envelope["response"])
            if (
                not isinstance(reply, dict)
                or set(reply) != {"pair"}
                or not isinstance(reply["pair"], str)
                or reply["pair"] not in PAIRS
            ):
                raise ValueError(
                    "The returned proposal is outside the registered comparisons."
                )
            selected = reply["pair"]
        finally:
            conn.close()
    return {
        "mode": mode,
        "selected": selected,
        "jobs": deepcopy(PAIRS[selected]),
        "submitted": False,
    }
