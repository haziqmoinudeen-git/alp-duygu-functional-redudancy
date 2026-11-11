# Data dictionary

## Reference data

`data/reference/replicate_results.csv` and `sweep_summary.csv` contain the 250-substitution reference results. The parameter manifest and reference history accompany them. Reference seed: ancestry 100000007, evolution 1220000107, N=100, loss-biased, chi=0.325, replicate 0. The trajectory is used in the catalytic-output and dynamics figures.

## Paired endpoint tables

`results/paper/{equal_total_250,equal_exposure_c5}/replicate_results.csv` has one row per seed and stopping regime. `experiment_id` prevents mixing the two regimes. `n_proteins`, `M`, `redundancy_r0`, `mutation_kernel`, `chi`, `replicate`, `ancestry_seed`, and `evolution_seed` identify conditions. `fixed_substitutions` is the recorded endpoint; `substitutions_per_protein` is this count divided by N. The `max_fixed` configuration field describes the shared trajectory target, not the row's endpoint.

Prefixes `initial_`, `final_`, `delta_` denote ancestor, endpoint, and endpoint minus ancestor. Metrics:

| Stem | Definition / units |
|---|---|
| mean_breadth | Mean inverse Simpson activity breadth; effective number of functions |
| strongest_share | Mean largest within-protein activity fraction; dimensionless |
| weakest_P | Minimum reaction performance; dimensionless |
| mean_redundancy | Mean inverse Simpson reaction contributor diversity; effective proteins |
| mean_total_activity | Mean sum of catalytic activity per protein; model activity units |
| mean_nonzero_functions | Mean count of positive protein activities |
| mean_active_functions | Mean count above 0.005 |
| mean_latent_loop_support | Mean latent reactions per loop |

`runtime_seconds` is wall time for the shared trajectory invocation, not the separate endpoint. On checkpoint resume it excludes earlier invocations. `completion_status` describes availability of that endpoint. Other configuration fields record Gamma shape, swaps, population, mutation rate, history cadence, and generation horizon. Aggregated columns end in `_count`, `_mean`, `_std` (sample SD, ddof=1).

## Trajectory NPZ

Each file contains `fixed` (integer substitution count), `time` (effective generations), `mean_breadth`, `mean_strongest_share`, `weakest_P`, `mean_total_activity`, `mean_redundancy`, `mean_nonzero_functions`, `mean_active_functions`, and `mean_latent_loop_support`. All histories begin at zero. `A_final` is the final N×120 activity matrix. `metadata` is JSON with full configuration, final counters, runtime, and completion status. NPZs can be read with `numpy.load(..., allow_pickle=False)`.

Manifest fields include environment, dependency versions, configuration, RNG, and seed registry. The repository includes the endpoint tables, complete histories, and publication figures.

Event counters `fixed_support_expansions`, `fixed_beneficial`, `fixed_neutralish`, and `fixed_deleterious` count fixed events to the row endpoint using a ±1e-12 classification. `counter_source` identifies direct history, final run summary, or the matched original 250-substitution CSV. `effective_generations` is the actual recorded time at the endpoint.
