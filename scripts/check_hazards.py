"""Compare cached loop hazards with a full recomputation."""

import inspect
import json
import sys
from pathlib import Path

import numpy as np

from enzyme_specialisation import model, simulation

source, start = inspect.getsourcelines(simulation.run_one_condition)
target = start + next(i for i, line in enumerate(source) if "if fixed % cfg.record_every" in line)
report = {}


def trace(frame, event, arg):
    if (
        event == "line"
        and frame.f_code.co_filename == simulation.__file__
        and frame.f_lineno == target
        and frame.f_locals["fixed"] == 1
    ):
        local = frame.f_locals
        cached = local["loop_hazards"].copy()
        fresh = np.array(
            [local["evaluate_loop_events"](i)["total_hazard"] for i in range(len(cached))]
        )
        diff = np.abs(cached - fresh)
        report.update(
            fixed=1,
            loops=len(cached),
            loops_differing=int(np.count_nonzero(diff)),
            max_absolute_hazard_difference=float(diff.max()),
            total_cached=float(cached.sum()),
            total_fresh=float(fresh.sum()),
            max_relative_hazard_difference=float(np.max(diff / np.maximum(fresh, 1e-300))),
        )
    return trace


cfg = model.RunConfig(
    40, "loss_biased", 0.325, 0, 40000007, 47, 8.0, 100, 1000000, 1e-8, 1, 3e10, 1, True
)
sys.settrace(trace)
try:
    model.run_one_condition(cfg)
finally:
    sys.settrace(None)
Path("reports/hazard_cache_comparison.json").write_text(json.dumps(report, indent=2))
print(report)
