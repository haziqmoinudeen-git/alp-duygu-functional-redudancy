# Regression validation

- **200/200** complete seeded histories; 40 parameter cells × five replicates.
- All five principal 250-substitution endpoint metrics match the supplied replicate CSV with **maximum absolute difference 0.0** across all 200 runs.
- All shared numerical columns in the regenerated 250-endpoint table, including event counters and effective generations, pass comparison against the source CSV at rtol=1e-11, atol=1e-12.
- The regenerated reference trajectory matches every supplied sampled history field at the recorded substitution counts (rtol=1e-12, atol=1e-12).
- Both endpoint tables contain 200 rows, with exactly five replicates per cell; the event-sign counters sum to the substitution endpoint.
- Histories contain every integer substitution count from zero to their target. Activity matrices satisfy bounds, and recorded times are nondecreasing.
- **14 tests pass**, including ancestor degree/activity/calibration/structure checks, mutation probability and magnitude checks, finite-effect fitness, fixation stability and preserved boundary behavior, scalar/vector structural equivalence, summary metrics, latent-support candidate restoration, deterministic reruns, source seed registry, event/hazard-array equivalence, and exact interrupted-checkpoint continuation.
- Selected deterministic event evaluations match the original candidate, selection, and hazard arrays byte for byte. Four short benchmark workloads and two complete 250-substitution benchmark trajectories retain identical final activity matrices.

These are implementation-regression results. They do not establish that the inherited hazard-cache scheme equals full recomputation; the independent audit finds a small discrepancy documented in `scientific_fidelity.md` and `hazard_cache_audit.json`.

`output_validation.json` is the machine-readable complete-grid result. `scripts/validate_outputs.py` reruns these output checks; `make test` reruns the automated test suite. The metadata-only fingerprint migration verifies that numerical array checksums remain unchanged and records them in `fingerprint_attestation.json`.

CITATION.cff passes cffconvert 2.0.0 validation against CFF schema 1.2.0. Local lint and formatting checks pass. The equivalent CI workflow is configured; no remote CI status is claimed.
