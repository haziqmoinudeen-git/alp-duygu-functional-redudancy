"""One-time provenance migration for outputs generated during this refactor.

Requires complete saved-reference regression and reruns scientific tests. Only
metadata changes; numerical array checksums are verified before and after writing.
"""

import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

from enzyme_specialisation.io import scientific_fingerprint

if __name__ == "__main__":
    validation = json.loads(Path("reports/output_validation.json").read_text())
    if validation["complete"] != 200 or validation["regression_errors"]:
        raise RuntimeError("Complete original-grid regression is required before attestation")
    subprocess.run([sys.executable, "-m", "pytest", "-q"], check=True)
    report = []
    for p in sorted(Path("results/paper/trajectories").glob("*.npz")):
        with np.load(p, allow_pickle=False) as z:
            data = {k: z[k] for k in z.files}
        meta = json.loads(str(data["metadata"]))
        if meta.get("scientific_fingerprint") == scientific_fingerprint():
            continue
        before = hashlib.sha256(
            b"".join(data[k].tobytes() for k in sorted(data) if k != "metadata")
        ).hexdigest()
        meta.update(
            scientific_fingerprint=scientific_fingerprint(),
            fingerprint_attestation="Post-refactor scientific tests and full saved-250-endpoint regression; metadata-only migration",
        )
        data["metadata"] = np.array(json.dumps(meta))
        tmp = p.with_suffix(".tmp")
        with tmp.open("wb") as f:
            np.savez_compressed(f, **data)
        with np.load(tmp, allow_pickle=False) as z:
            after = hashlib.sha256(
                b"".join(z[k].tobytes() for k in sorted(z.files) if k != "metadata")
            ).hexdigest()
        assert before == after
        os.replace(tmp, p)
        report.append({"file": p.name, "numerical_arrays_sha256": before})
    Path("reports/fingerprint_attestation.json").write_text(
        json.dumps(
            {"scientific_fingerprint": scientific_fingerprint(), "outputs": report}, indent=2
        )
    )
