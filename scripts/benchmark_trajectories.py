"""Time two complete 250-substitution trajectories against the scalar reference implementation."""

import importlib.util
import json
import resource
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from enzyme_specialisation import model

spec = importlib.util.spec_from_file_location(
    "original_full", Path("data/provenance/original_sweep.py")
)
old = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = old
spec.loader.exec_module(old)
rows = []
for n, k, chi in [(40, "loss_biased", 0.0), (120, "symmetric", 1.0)]:
    c = model.RunConfig(
        n, k, chi, 0, n * 1000000 + 7, 47, 8.0, 50000, 1000000, 1e-8, 250, 3e10, 1, True
    )
    times = []
    results = []
    for mod in [old, model]:
        start = time.perf_counter()
        results.append(mod.run_one_condition(mod.RunConfig(**asdict(c))))
        times.append(time.perf_counter() - start)
    np.testing.assert_array_equal(results[0]["A_final"], results[1]["A_final"])
    for key in results[0]["history"]:
        np.testing.assert_allclose(
            results[0]["history"][key], results[1]["history"][key], rtol=1e-13, atol=1e-14
        )
    rows.append(
        {
            "N": n,
            "kernel": k,
            "chi": chi,
            "substitutions": 250,
            "baseline_seconds": times[0],
            "optimised_seconds": times[1],
            "baseline_seconds_per_fixation": times[0] / 250,
            "optimised_seconds_per_fixation": times[1] / 250,
            "speedup": times[0] / times[1],
            "peak_process_MiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            "identical_final_activity": True,
        }
    )
    Path("reports/full_trajectory_benchmark.json").write_text(json.dumps(rows, indent=2))
    print(rows[-1], flush=True)
