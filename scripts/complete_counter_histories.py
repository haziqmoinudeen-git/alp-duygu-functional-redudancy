"""Fill event counters for the cN endpoint preceding 250, preserving states exactly."""

import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import yaml

from enzyme_specialisation import model
from enzyme_specialisation.experiments import identifier, registry


def fill(cfg):
    p = Path("results/paper/trajectories") / f"{identifier(cfg)}.npz"
    with np.load(p) as z:
        if "fixed_beneficial" in z:
            return
        before = {k: z[k].copy() for k in z.files}
        meta = json.loads(str(z["metadata"]))
    start = time.perf_counter()
    r = model.run_one_condition(cfg)
    for key in r["history"]:
        if key in before:
            np.testing.assert_array_equal(r["history"][key], before[key])
    np.testing.assert_array_equal(r["A_final"], before["A_final"])
    meta["counter_history_rerun_seconds"] = time.perf_counter() - start
    tmp = p.with_suffix(".tmp")
    with tmp.open("wb") as f:
        np.savez_compressed(f, **r["history"], A_final=r["A_final"], metadata=json.dumps(meta))
    os.replace(tmp, p)
    return identifier(cfg)


if __name__ == "__main__":
    c = yaml.safe_load(Path("configs/paper.yaml").read_text())
    jobs = [j for j in registry(c) if j.n_proteins == 40]
    with ProcessPoolExecutor(max_workers=4) as pool:
        for result in pool.map(fill, jobs):
            print(result, flush=True)
