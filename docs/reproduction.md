# Reproduction

Run commands from the repository root with Python 3.12 and uv:

```sh
make install
make test
make smoke
make reproduce
make validate
make figures
```

`uv.lock` pins the environment. For a conventional environment: `python -m pip install -e '.[dev]'` (this does not enforce the lock).

`python -m enzyme_specialisation simulate --config configs/paper.yaml --workers 4` runs 200 histories to `max(250,5N)`. Both stopping-regime config files deliberately run this shared pipeline. Their results have distinct experiment IDs and directories. `--limit 1` supplies a single-condition validation run; it does not represent the full experiment. `configs/extended.yaml` extends to 20 replicates. Replicates 0–4 retain the original full-grid seed registry; 5–19 use SHA-256-derived seeds. Subsetting or reordering the original grid changes its historical condition index, so preserve the paper grid when reproducing original seeds.

Completed NPZs are validated before skipping. Checkpoints every 25 substitutions contain the full activity and support states, cached hazards, RNG state, counters, and history. Writes use atomic replacement. Restart with the same command and configuration. Loading checkpoints uses pickle: only load trusted local checkpoints. An interrupted interval repeats at most 24 substitutions. Resume currently repeats ancestor setup before restoring the checkpoint; this costs setup time but preserves continuation exactly.

The generation horizon is retained. A run that reaches it before the substitution target is marked incomplete and fails complete-grid validation. Worker count is configurable and capped at half the available CPUs, with a conservative 512 MiB memory allowance per worker and 512 MiB reserve. Available memory is recorded in the manifest. The driver has finite stopping criteria and bounded ancestral construction. Completed-job progress includes elapsed time, ETA, and failures; intermediate substitutions are recoverable from checkpoints. The driver does not automatically retry numerical failures.

Histories retain every fixed event. Exposure curves use `searchsorted(..., side="right")-1` on a common 201-point grid from 0 to 5: the last observed state at or before each exposure is carried forward. There is no interpolation or invented event. Ribbons show sample standard deviation across trajectories, not confidence intervals or standard errors. Only trajectories reaching exposure 5 contribute.

Canonical figures use the supplied 250-event CSV and the regenerated reference trajectory when available; the new trajectory-derived 250 endpoints are independently compared to them. The saved every-ten-substitution reference history is independently checked. New histories retain every event.

Completed results and checkpoints carry a fingerprint of the scientific source files. A changed scientific source invalidates cache reuse. The original results generated before fingerprinting were migrated only after scientific tests and full endpoint regression; numerical-array hashes were verified unchanged. `scripts/attest_legacy_outputs.py` documents that one-time migration and is not part of ordinary reproduction.
