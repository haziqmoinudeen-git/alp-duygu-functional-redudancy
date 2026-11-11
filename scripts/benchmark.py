"""Profile scalar and vectorised expansion evaluation on matched workloads."""

import cProfile
import importlib.util
import json
import pstats
import resource
import sys
import time
from dataclasses import asdict
from pathlib import Path

import numpy as np

from enzyme_specialisation import model

s = importlib.util.spec_from_file_location("original", Path("data/provenance/original_sweep.py"))
old = importlib.util.module_from_spec(s)
sys.modules["original"] = old
s.loader.exec_module(old)
rows = []
for n, k, chi in [
    (40, "loss_biased", 0.0),
    (120, "symmetric", 1.0),
    (40, "symmetric", 0.1),
    (120, "loss_biased", 0.325),
]:
    c = model.RunConfig(
        n, k, chi, 0, 1000000 * n + 7, 47, 8.0, 50000, 1000000, 1e-8, 3, 3e10, 1, True
    )
    outputs = []
    times = []
    for mod in [old, model]:
        cfg = mod.RunConfig(**asdict(c))
        start = time.perf_counter()
        profile = cProfile.Profile()
        profile.enable()
        r = mod.run_one_condition(cfg)
        profile.disable()
        times.append(time.perf_counter() - start)
        outputs.append(r)
        if mod is old and n == 40 and k == "loss_biased":
            with open("reports/profile_baseline.txt", "w") as f:
                pstats.Stats(profile, stream=f).sort_stats("cumtime").print_stats(25)
    for key in outputs[0]["history"]:
        np.testing.assert_allclose(
            outputs[0]["history"][key], outputs[1]["history"][key], rtol=1e-13, atol=1e-14
        )
    np.testing.assert_array_equal(outputs[0]["A_final"], outputs[1]["A_final"])
    rows.append(
        dict(
            N=n,
            kernel=k,
            chi=chi,
            substitutions=3,
            baseline_seconds=times[0],
            optimised_seconds=times[1],
            speedup=times[0] / times[1],
            peak_process_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024,
            identical_final_activity=True,
        )
    )
    Path("reports/benchmark.json").write_text(json.dumps(rows, indent=2))
    print(rows[-1], flush=True)
Path("reports/performance_benchmark.md").write_text(
    "# Measured benchmark\n\nThree-substitution workloads including ancestry and catalogue setup. Profile overhead applies to both versions. Peak RSS is cumulative process peak, not isolated per case. These timings cover the measured three-substitution workloads.\n\n"
    + "\n".join(
        f"- N={r['N']}, {r['kernel']}, chi={r['chi']}: {r['baseline_seconds']:.3f}s → {r['optimised_seconds']:.3f}s; {r['speedup']:.2f}×; optimised {r['optimised_seconds'] / 3:.3f}s/substitution including setup; peak {r['peak_process_MiB']:.1f} MiB. Final activity identical; all recorded metrics agree to 1e-13 relative tolerance."
        for r in rows
    )
    + "\n"
)
