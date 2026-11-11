"""Fingerprint the scientific source to prevent reuse across code changes."""

import hashlib
from functools import cache
from pathlib import Path


@cache
def scientific_fingerprint() -> str:
    root = Path(__file__).parent
    digest = hashlib.sha256()
    for name in [
        "parameters",
        "ancestry",
        "catalogue",
        "fitness",
        "fixation",
        "structure",
        "metrics",
        "simulation",
    ]:
        digest.update(name.encode())
        digest.update((root / f"{name}.py").read_bytes())
    return digest.hexdigest()
