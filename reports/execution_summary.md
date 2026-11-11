# Execution summary

## Delivered computation

- 200/200 primary trajectories complete, five seeds per each of 40 cells.
- 82,000 fixed substitutions executed in shared histories to max(250,5N).
- Two 200-row endpoint tables and two 40-cell aggregate tables under `results/paper/`.
- Complete event histories and final activity matrices in 200 compressed NPZs (approximately 6.8 MiB including tables/manifests).
- Original first-five ancestry/evolution seeds preserved; see the manifests and endpoint tables for the complete registry.
- All 250-substitution principal metrics match saved outputs exactly; additional shared numerical output columns and the reference history pass regression.
- No numerical simulation failure reported. One execution process was interrupted after 199 completions; the final trajectory resumed from its checkpoint. The reproduction command skipped those 199 completed runs and finished the remaining trajectory in 7.74 seconds, including repeated setup.

## Runtime and resources

The first invocation's last completed-job log records 199 completions at 693.8 seconds. Its execution process ended before writing a final batch runtime. Together with the resumed invocation this is approximately twelve minutes of primary execution; exact end-to-end orchestration time is not claimed. Per-trajectory recorded wall times sum to 2071.89 seconds and overlap across workers; the resumed trajectory's metadata covers its resumed invocation only. Six CPUs were visible; the initial request for four workers was capped at three by the configured half-core policy. Later cached reproduction used two workers.

Forty N=40 trajectories were rerun solely to retain intermediate event counters at the earlier 200-event endpoint. Their numerical histories and activity matrices were compared exactly to the already saved results. Their additional per-worker wall times sum to 277.79 seconds, overlapping with primary execution. These runs are disclosed rather than counted as new independent replicates.

Four short profiled benchmark cases and two complete 250-event benchmark trajectories ran separately. The full workloads show measured 10.30× and 8.48× acceleration, with identical final activity arrays; peak benchmark-process RSS was approximately 107 MiB. See `performance_benchmark.md` for limitations.

## Reproduction and outputs

`make install`, `make test`, `make smoke`, `make reproduce`, `make validate`, and plotting commands were executed successfully. A subsequent full reproduction invocation skipped all 200 validated trajectories. Atomic files and exact checkpoint continuation are tested. Source fingerprints invalidate cached trajectories/checkpoints after scientific source changes.

Six main presentations (including absolute and change exposure curves) and three alternatives have PDF/SVG/600-dpi PNG exports, with lower-resolution README previews. Main Figures 1 and 2 retain their canonical filenames. Figure/data checksums appear in `artifact_manifest.json`.

The initial simulation invocation began before a Git commit existed; its original manifest is preserved as `manifest.initial.json`. A later committed-source validation/reproduction manifest records the Git commit and scientific source fingerprint. Numerical outputs were attested after refactoring through tests and full saved-endpoint regression; their array checksums were unchanged. This provenance distinction is explicit.

## Conditional work and release status

The optional 20-seed expansion was not executed. Its additional 600 trajectories would require roughly another 35 minutes at the observed primary throughput (a projection, not a benchmark). The fully implemented `configs/extended.yaml` retains original seeds, adds deterministic SHA-256-derived seeds, and skips completed trajectories. Five seeds per cell are the completed analysis; no 20-seed result is claimed.

The local repository is prepared for GitHub review, with original files untouched, tests and CI configured, locked dependencies, citation metadata, documentation, and committed compact data. No remote was published. Public release still requires author licensing decisions and review of the inherited scientific discrepancies in `scientific_fidelity.md`. A cloud GitHub Actions run has not been observed; the equivalent local checks passed.
