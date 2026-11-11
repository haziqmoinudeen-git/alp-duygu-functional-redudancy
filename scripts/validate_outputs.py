"""Validate completeness, paired endpoint provenance, and historical regression."""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from enzyme_specialisation.experiments import aggregate, identifier, registry, valid

c = yaml.safe_load(Path("configs/paper.yaml").read_text())
root = Path(c["output"])
jobs = registry(c)
complete = [j for j in jobs if valid(root / "trajectories" / f"{identifier(j)}.npz", j)]
aggregate(root, c["exposure"])
ref = pd.read_csv("data/reference/replicate_results.csv")
errors = []
maxdiff = 0.0
for cfg in complete:
    path = root / "trajectories" / f"{identifier(cfg)}.npz"
    with np.load(path, allow_pickle=False) as h:
        np.testing.assert_array_equal(h["fixed"], np.arange(cfg.max_fixed + 1))
        assert np.all(np.diff(h["time"]) >= 0)
        assert np.all((h["A_final"] >= 0) & (h["A_final"] <= 1 + 1e-15))
        r = ref[
            (ref.n_proteins == cfg.n_proteins)
            & (ref.mutation_kernel == cfg.mutation_kernel)
            & (ref.chi == cfg.chi)
            & (ref.replicate == cfg.replicate)
        ].iloc[0]
        for key, col in [
            ("mean_breadth", "final_mean_breadth"),
            ("mean_strongest_share", "final_strongest_share"),
            ("weakest_P", "final_weakest_P"),
            ("mean_redundancy", "final_mean_redundancy"),
            ("mean_total_activity", "final_mean_total_activity"),
        ]:
            d = abs(h[key][250] - r[col])
            maxdiff = max(maxdiff, d)
            if not np.isclose(h[key][250], r[col], rtol=1e-11, atol=1e-12):
                errors.append({"id": identifier(cfg), "metric": key, "difference": d})
result = {
    "complete": len(complete),
    "expected": len(jobs),
    "regression_errors": errors,
    "max_absolute_endpoint_difference": maxdiff,
}
Path("reports/output_validation.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
if errors or len(complete) != len(jobs):
    raise SystemExit(1)
