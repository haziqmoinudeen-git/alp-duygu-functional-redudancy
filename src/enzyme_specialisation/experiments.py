"""Stable seed registry and atomic, resumable trajectory execution."""

import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import pandas as pd

from . import model
from .io import scientific_fingerprint

METRICS = {
    "mean_breadth": "mean_breadth",
    "strongest_share": "mean_strongest_share",
    "weakest_P": "weakest_P",
    "mean_redundancy": "mean_redundancy",
    "mean_total_activity": "mean_total_activity",
}


def identifier(c):
    return f"N{c.n_proteins}_{c.mutation_kernel}_chi{c.chi:g}_rep{c.replicate:03d}"


def registry(c):
    original = model.build_jobs(
        c["redundancies"],
        c["kernels"],
        c["chis"],
        5,
        c["base_seed"],
        8.0,
        c["network_swaps"],
        1000000,
        1e-8,
        250,
        c["max_generations"],
        1,
    )
    jobs = []
    for rep in range(c["replicates"]):
        for cfg in original:
            if cfg.replicate != rep % 5:
                continue
            if rep >= 5:
                d = hashlib.sha256(f"extension:{identifier(cfg)}:{rep}".encode()).digest()
                cfg = replace(
                    cfg,
                    replicate=rep,
                    ancestry_seed=int.from_bytes(
                        hashlib.sha256(
                            f"ancestry:{c['base_seed']}:{cfg.n_proteins}:{rep}".encode()
                        ).digest()[:4],
                        "little",
                    ),
                    evolution_seed=int.from_bytes(d[4:8], "little"),
                )
            jobs.append(
                replace(
                    cfg,
                    max_fixed=c.get("smoke_fixed") or max(250, int(c["exposure"] * cfg.n_proteins)),
                    save_history=True,
                )
            )
    return jobs


def atomic_json(p, v):
    t = Path(str(p) + ".tmp")
    t.write_text(json.dumps(v, indent=2))
    os.replace(t, p)


def valid(p, c):
    try:
        with np.load(p, allow_pickle=False) as z:
            return (
                json.loads(str(z["metadata"]))["config"] == asdict(c)
                and json.loads(str(z["metadata"])).get("scientific_fingerprint")
                == scientific_fingerprint()
                and int(z["fixed"][-1]) == c.max_fixed
                and np.isfinite(z["mean_breadth"]).all()
            )
    except (OSError, ValueError, KeyError):
        return False


def worker(cfg, output):
    path = Path(output) / "trajectories" / f"{identifier(cfg)}.npz"
    if valid(path, cfg):
        return identifier(cfg), 0.0, "cached"
    start = time.perf_counter()
    checkpoint = path.with_suffix(".checkpoint")
    r = model.run_one_condition(cfg, checkpoint)
    seconds = time.perf_counter() - start
    status = "complete" if r["summary"]["fixed_substitutions"] == cfg.max_fixed else "incomplete"
    meta = {
        "config": asdict(cfg),
        "scientific_fingerprint": scientific_fingerprint(),
        "summary": r["summary"],
        "runtime_seconds": seconds,
        "status": status,
    }
    temp = path.with_suffix(".tmp")
    with temp.open("wb") as f:
        np.savez_compressed(f, **r["history"], A_final=r["A_final"], metadata=json.dumps(meta))
    os.replace(temp, path)
    if status == "complete":
        checkpoint.unlink(missing_ok=True)
    return identifier(cfg), seconds, status


def run(config, workers=None, limit=None):
    output = Path(config["output"])
    (output / "trajectories").mkdir(parents=True, exist_ok=True)
    jobs = registry(config)
    if limit:
        jobs = jobs[:limit]
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL, text=True
        ).strip()
    except subprocess.CalledProcessError:
        commit = None
    atomic_json(
        output / "manifest.json",
        {
            "configuration": config,
            "python": platform.python_version(),
            "platform": platform.platform(),
            "architecture": platform.machine(),
            "cpu_count": os.cpu_count(),
            "git_commit": commit,
            "scientific_fingerprint": scientific_fingerprint(),
            "requested_workers": workers or config["workers"],
            "effective_workers": min(
                workers or config["workers"], max(1, (os.cpu_count() or 2) // 2), len(jobs)
            ),
            "rng": "numpy.PCG64",
            "dependencies": {
                n: importlib.metadata.version(n)
                for n in ["numpy", "pandas", "matplotlib", "pyyaml"]
            },
            "jobs": [asdict(c) for c in jobs],
        },
    )
    start = time.perf_counter()
    failures = []
    cached = 0
    with ProcessPoolExecutor(
        max_workers=min(workers or config["workers"], max(1, (os.cpu_count() or 2) // 2), len(jobs))
    ) as pool:
        pending = {pool.submit(worker, c, str(output)): c for c in jobs}
        for done, f in enumerate(as_completed(pending), 1):
            cfg = pending[f]
            try:
                name, seconds, status = f.result()
                cached += int(status == "cached")
                if status == "incomplete":
                    failures.append({"id": name, "error": "Stopping target not reached"})
                elapsed = time.perf_counter() - start
                print(
                    f"{done}/{len(jobs)} {name} fixed={cfg.max_fixed} {status} {seconds:.1f}s elapsed={elapsed:.1f}s ETA={elapsed * (len(jobs) - done) / done:.1f}s failures={len(failures)}",
                    flush=True,
                )
            except Exception as exc:
                failures.append({"id": identifier(cfg), "error": repr(exc)})
                print(f"FAILED {failures[-1]}", flush=True)
    atomic_json(
        output / "execution.json",
        {
            "elapsed_seconds": time.perf_counter() - start,
            "failures": failures,
            "cached_trajectories": cached,
            "requested_workers": workers or config["workers"],
        },
    )
    aggregate(output, config["exposure"])
    if failures:
        raise RuntimeError(
            f"{len(failures)} trajectories failed or stopped early; see execution.json"
        )


def aggregate(output, exposure=5):
    output = Path(output)
    records = []
    reference = pd.read_csv("data/reference/replicate_results.csv")
    counter_names = (
        "fixed_support_expansions",
        "fixed_beneficial",
        "fixed_neutralish",
        "fixed_deleterious",
    )
    for path in sorted((output / "trajectories").glob("*.npz")):
        with np.load(path, allow_pickle=False) as h:
            meta = json.loads(str(h["metadata"]))
            cfg = model.RunConfig(**meta["config"])
            for experiment, end in [
                ("equal_total_250", 250),
                ("equal_exposure_c5", int(exposure * cfg.n_proteins)),
            ]:
                where = np.flatnonzero(h["fixed"] == end)
                if not len(where):
                    continue
                idx = where[0]
                row = {
                    **meta["config"],
                    "experiment_id": experiment,
                    "fixed_substitutions": end,
                    "substitutions_per_protein": end / cfg.n_proteins,
                    "redundancy_r0": model.redundancy_from_n(cfg.n_proteins),
                    "M": model.M_FUNCTIONS,
                    "eta": model.eta_from_chi(cfg.chi),
                    "stop_reason": "recorded fixed-substitution endpoint",
                    "runtime_seconds": meta["runtime_seconds"],
                    "completion_status": "complete",
                }
                for col, key in METRICS.items():
                    row["initial_" + col] = float(h[key][0])
                    row["final_" + col] = float(h[key][idx])
                    row["delta_" + col] = float(h[key][idx] - h[key][0])
                for key in [
                    "mean_nonzero_functions",
                    "mean_active_functions",
                    "mean_latent_loop_support",
                ]:
                    row["final_" + key] = float(h[key][idx])
                row["effective_generations"] = float(h["time"][idx])
                matched = reference[
                    (reference.ancestry_seed == cfg.ancestry_seed)
                    & (reference.evolution_seed == cfg.evolution_seed)
                    & (reference.n_proteins == cfg.n_proteins)
                    & (reference.mutation_kernel == cfg.mutation_kernel)
                    & (reference.chi == cfg.chi)
                ]
                if all(key in h for key in counter_names):
                    row.update({key: int(h[key][idx]) for key in counter_names})
                    row["counter_source"] = "recorded_history"
                elif end == h["fixed"][-1]:
                    row.update({key: int(meta["summary"][key]) for key in counter_names})
                    row["counter_source"] = "recorded_final_summary"
                elif end == 250 and len(matched) == 1:
                    row.update({key: int(matched.iloc[0][key]) for key in counter_names})
                    row["counter_source"] = "matched_original_250_csv"
                else:
                    raise ValueError(f"Missing intermediate event counters: {path} at {end}")
                for key in [
                    "mean_nonzero_functions",
                    "mean_active_functions",
                    "mean_latent_loop_support",
                ]:
                    row["initial_" + key] = float(h[key][0])
                    row["delta_" + key] = float(h[key][idx] - h[key][0])
                records.append(row)
    if not records:
        return
    df = pd.DataFrame(records)
    for experiment, sub in df.groupby("experiment_id"):
        folder = output / experiment
        folder.mkdir(exist_ok=True)
        temp = folder / "replicate_results.csv.tmp"
        sub.to_csv(temp, index=False)
        os.replace(temp, folder / "replicate_results.csv")
        metrics = [c for c in sub if c.startswith(("initial_", "final_", "delta_"))]
        g = sub.groupby(["mutation_kernel", "chi", "redundancy_r0"])[metrics].agg(
            ["count", "mean", "std"]
        )
        g.columns = ["_".join(c) for c in g.columns]
        temp = folder / "sweep_summary.csv.tmp"
        g.to_csv(temp)
        os.replace(temp, folder / "sweep_summary.csv")
