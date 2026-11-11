"""Command line entry points."""

import argparse

import yaml

from .experiments import aggregate, run


def main():
    p = argparse.ArgumentParser()
    p.add_argument("command", choices=["simulate", "plot", "validate", "aggregate"])
    p.add_argument("--config", default="configs/paper.yaml")
    p.add_argument("--workers", type=int)
    p.add_argument("--limit", type=int)
    p.add_argument("--experiment", default="all")
    a = p.parse_args()
    with open(a.config) as f:
        c = yaml.safe_load(f)
    if a.command == "simulate":
        run(c, a.workers, a.limit)
    elif a.command == "aggregate":
        aggregate(c["output"], c["exposure"])
    elif a.command == "plot":
        from .plotting import plot_all

        plot_all(c["output"], a.experiment)
    else:
        from pathlib import Path

        from .experiments import identifier, registry, valid

        jobs = registry(c)
        good = sum(
            valid(Path(c["output"]) / "trajectories" / f"{identifier(j)}.npz", j) for j in jobs
        )
        print(f"Validated {good}/{len(jobs)} complete trajectories")
        if good != len(jobs):
            raise SystemExit(1)
