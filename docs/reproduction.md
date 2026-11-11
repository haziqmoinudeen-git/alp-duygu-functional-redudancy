# Reproduction

## Environment and commands

Run commands from the repository root with uv. Python 3.12.7 and `uv.lock` specify the environment.

```bash
make install
make test
make smoke
make reproduce
make validate
make figures
```

`make reproduce` executes the paired experiment, validates the data, and generates figures. Completed trajectory files are checked and reused. Conventional installation is also supported with `python -m pip install -e '.[dev]'`; use uv for the locked environment.

## Experiment configurations

```bash
uv run python -m enzyme_specialisation simulate --config configs/paper.yaml --workers 4
uv run python -m enzyme_specialisation plot --experiment equal_total_250
uv run python -m enzyme_specialisation plot --experiment equal_exposure_c5
```

Both endpoint configurations run shared histories to `max(250,5N)` and write separate analysis directories. `--limit 1` selects a single trajectory for a small validation run. `configs/smoke_test.yaml` uses 100 network swaps and two substitutions; `configs/paper.yaml` uses the study parameters. `configs/extended.yaml` increases the replicate count to 20 per cell.

## Seeds and parallel execution

Replicates 0–4 use the five-replicate seed registry; replicates 5–19 use SHA-256-derived seeds. Each redundancy/replicate ancestor is shared across mutation kernels and incompatibility values. Evolution seeds identify individual trajectories. Scheduling does not determine seeds. The paper grid order determines the reference condition indices, so use that grid to reproduce the recorded seed assignments.

Worker count is configurable and capped at half the available CPUs. Memory limits reserve 512 MiB for the parent process and allow approximately 512 MiB per worker. The manifest records available resources. Progress messages give completed trajectories, condition, substitution target, elapsed time, estimated remaining time, and failures.

## Checkpointing

Checkpoints every 25 substitutions contain activity, latent support, cached hazards, RNG state, counters, and history. Restart with the same command and configuration to resume. Atomic replacement protects both checkpoints and completed outputs. Source fingerprints prevent cache reuse after scientific code changes.

An interrupted interval repeats at most 24 substitutions. Resuming performs ancestor setup before restoring the checkpoint state. Checkpoints use pickle and should be loaded only from trusted local runs. A trajectory reaching the $3\times10^{10}$-generation horizon before its substitution target is marked incomplete. Numerical failures are recorded for diagnosis.

## Figure data and exposure alignment

Endpoint tables and complete event histories provide the figure data. The 250-substitution reference table and sampled reference history are numerical comparison fixtures.

Exposure curves use a common 201-point grid from 0 to 5. `searchsorted(..., side="right")-1` selects the latest recorded state at or before each exposure. Ribbons show sample standard deviation across aligned trajectories. The [data dictionary](data_dictionary.md) defines columns, units, and runtime fields; [numerical methods](numerics.md) describe event rates and hazard caching.
