# Numerical validation

The complete dataset contains 200 trajectories and five replicates in each of 40 parameter cells. All five principal 250-substitution endpoint metrics match `data/reference/replicate_results.csv` with maximum absolute difference 0.0.

Additional numerical comparisons use the following tolerances:

- Shared endpoint columns, including event counters and effective generations: rtol=1e-11, atol=1e-12.
- Reference history at recorded substitution counts: rtol=1e-12, atol=1e-12.
- History fields in matched complete-trajectory benchmarks: rtol=1e-13, atol=1e-14. Final activity matrices are identical.

Fourteen automated tests cover ancestral degrees and calibration, activity bounds, structural feasibility, mutation catalogue probabilities and magnitudes, finite-effect fitness, fixation probabilities, summary metrics, latent-support restoration, seed assignment, deterministic event/hazard comparisons, reproducibility, and exact checkpoint continuation. Histories contain every substitution count through their target, have nondecreasing times, and satisfy activity bounds. Endpoint event-sign counters sum to the fixed-substitution count.

`make test` runs the tests. `make validate` checks the complete recorded dataset and writes [machine-readable validation results](output_validation.json). The [numerical methods](../docs/numerics.md) describe hazard caching and other algorithmic conventions; the [hazard comparison](hazard_cache_comparison.json) quantifies cached versus full-refresh values on a diagnostic workload.
