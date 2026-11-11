# Scientific fidelity audit

## Authority and provenance

The canonical publication implementation is `data/provenance/original_sweep.py`, an unchanged copy of `integrated_reaction_theory_v2_fast_sweep.py`. The manuscript PDF *Integrated_Reaction_Theory_V2 (3).pdf* was extracted and read, including Methods II A–F, Results III A–C, and both captions. Its parameters, metric definitions, quoted reference trajectory, and cell-mean ranges agree with the source and supplied results. `source_inventory.json` records relevant originals and checksums. The initial preliminary audit's statement that the manuscript was unavailable is superseded: the PDF was present and read during continuation.

## Preserved science

Fitness, aggregate reaction performance, equal weights, loop activity decomposition, ancestral degree randomisation, Gamma initialisation, K calibration, structural feasibility and chi mapping, activity bounds, mutation sign/magnitude distributions, catalogue normalisation, support expansion, latent support, exact finite-effect selection, origin–fixation sampling, numerical fixation handling, generation horizon, and summary metrics retain the original operations. No duplication, structural-investment model, neutral evolutionary control, penalty, alternative distribution, or additional biological assumption was introduced. The symmetric kernel remains a mutation condition under selection.

Code is separated into parameters, ancestry, catalogue, fitness, fixation, structure, metrics, and simulation modules. The original source remains available for direct comparison. Engineering changes batch scalar expansion evaluations, add atomic checkpoint/resume, retain every-event histories, derive paired endpoints, and provide deterministic experiment execution and plotting. RNG draw order is retained. Deterministic candidate and hazard comparison tests compare the returned arrays byte for byte against the original evaluation on a selected workload. Four profiled short workloads retain identical final activity matrices.

The historical five-seed registry is preserved. Exact continuation from an interrupted checkpoint is tested. Full-grid numerical results are reported separately in `output_validation.json` and `regression_validation.md`; no stronger equivalence claim should be inferred from a short benchmark alone.

## Authorised changes

Changes.docx requests retention of 250 substitutions, addition of cN with c=5, optional replicate expansion to 10–20, and breadth against x=n_fix/N. The implementation executes each history to max(250,5N), records genuine states at both endpoints, uses the same metric definitions, and aligns exposure curves by last-observation-carried-forward without interpolation. No new model condition was added.

## Inherited discrepancies — not silently corrected

1. **Hazard dependency refresh and expansion events.** The baseline updates cached total hazards for the selected protein and loops whose current latent support overlaps changed reactions. A support-expandable loop can target *any* reaction outside its current support. Its expansion selection coefficients can therefore change when a reaction outside that support changes. Those cached total hazards need not be refreshed by the baseline dependency rule. The chosen loop is evaluated afresh before selecting its event, while loop selection used the cache. This can produce an inconsistency between loop selection and event selection. The refactor preserves this behavior to reproduce the existing scientific baseline. It is an inherited implementation issue, not an optimisation introduced here. The pipeline must not be described as verified against complete hazard recomputation; the author should decide whether a separately versioned scientific correction and rerun are needed.

2. **Catalogue versus coordinate independence.** The manuscript describes independent coordinate signs. The original catalogue conditions on at least one nonzero sign and draws one shared negative magnitude and one shared positive magnitude for all coordinates with each sign. This exact catalogue is retained; it must not be replaced by independent per-coordinate magnitudes. The conditioning and magnitude coupling deserve explicit Methods wording.

3. **Origin-rate naming.** The source uses Ne × U_loop divided by the total number of loops for each loop's origin rate. This preserves total origin intensity across N, rather than taking U_loop literally as a per-loop rate. No rescaling has been made. Endpoint substitution sequences and exposure comparisons use the original dynamics; generation-time interpretation requires clarification.

4. **Numerical boundaries.** The source treats |s|<1e-14 as neutral but its positive/negative masks use strict inequalities. Exactly ±1e-14 falls into neither mask and returns zero. This rare boundary behavior is preserved and documented rather than altered under the guise of stability.

These audit items prevent an unconditional claim that every mathematical transition probability matches a complete recomputation of the written model. They do not prevent reproducing the supplied implementation and performing the authorised endpoint/exposure extensions. Author review is a release gate.

The independent diagnostic `scripts/audit_hazards.py` quantifies item 1 after the first substitution on a small validation workload. Of 160 loops, 107 cached total hazards differ from full recomputation; maximum absolute difference is 1.04e-14 and maximum relative difference 3.70e-6. Total cached/fresh hazards and exact workload details are recorded in `hazard_cache_audit.json` and the script. This is a demonstrable small discrepancy, not a hypothetical warning.
