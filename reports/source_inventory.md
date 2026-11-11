# Source inventory and provenance map

Original research files remain untouched. Relevant files were copied into this repository; no large datasets were duplicated. `source_inventory.json` contains source sizes and SHA-256 digests. Unrelated workspace filenames and private dataset inventories are excluded from this release directory.

| Original material | Classification | Decision |
|---|---|---|
| integrated_reaction_theory_v2_fast_sweep.py | 1: canonical implementation | Unchanged copy in data/provenance; extracted scientific modules preserve its operations |
| irt_v2_sweep_results/replicate_results.csv | 2: canonical publication data | Copied; 200 unique conditions at 250 substitutions |
| irt_v2_sweep_results/sweep_summary.csv | 2: canonical publication data | Copied as reference aggregate |
| irt_v2_sweep_results/reference_history.csv | 2: canonical publication data | Copied; reference replicate 0, every 10 substitutions |
| irt_v2_sweep_results/run_manifest.json | 2: canonical publication provenance | Copied; paper grid, alpha=8, 50,000 swaps |
| irt_v2_sweep_results/*.png | 2: previous individual plots | Retained in source workspace; regenerated assembled figures from numerical data |
| Integrated_Reaction_Theory_V2 (3).pdf | Manuscript authority | Methods, Results, captions read; extracted text retained in reports/manuscript.txt |
| Changes.docx | Extension specification | Copied and transcribed in full |
| original_model_neutral_control.py, original_model_neutral_control_fast.py | 4: alternative scientific experiments | Fixed earlier architecture/true function loss and neutral comparisons; excluded from pipeline |
| loop_local_random_overlap_very_strong_loss_30B_with_nulls.py | 4: alternative scientific experiments | Explicit neutral/null modes and previous activity architecture; excluded |
| loss_biased_structural_budget_break_test.py; my_break_test/ | 4: alternative scientific experiments | Inert reservoir and nonlinear structure-to-activity mapping; excluded |
| irt_neutral_seed7.npz; selected_vs_neutral_plots/ | 4: alternative experiment outputs | Neutral comparison provenance only; excluded |
| analysis.ipynb | 5: unrelated material | Cell inspection identifies COBRA microbial interaction analysis; excluded |
| saved_results/ | 5: unrelated material | Schemas/names identify microbial resource competition and formulation sweeps; excluded, not copied |
| Explicitly excluded PHLORA and SCN1A materials | 5: unrelated material | Left untouched; not part of publication pipeline |

The canonical source identifies the same paper, activity equations, chi mapping, ancestral architecture, mutation parameters, and factorial design as the manuscript. The original 200-row data have five replicates for each of 40 cells, all at 250 substitutions. The reference condition matches manuscript quoted values (breadth approximately 13.51→11.28; strongest share approximately 0.113→0.188). Ancestry seed 100000007 and evolution seed 1220000107 are taken from the original registry and CSV, not selected by visual appearance.

No archival alternatives were copied into the scientific package. The only historical executable included is the actual canonical source needed for regression and benchmark comparisons.
