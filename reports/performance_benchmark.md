# Measured benchmark

Three-substitution workloads including ancestry and catalogue setup. Profile overhead applies to both versions. Peak RSS is cumulative process peak, not isolated per case. Complete-trajectory timings are reported separately below.

- N=40, loss_biased, chi=0.0: 6.993s → 2.717s; 2.57×; optimised 0.906s/substitution including setup; peak 105.5 MiB. Final activity identical; all recorded metrics agree to 1e-13 relative tolerance.
- N=120, symmetric, chi=1.0: 14.053s → 3.122s; 4.50×; optimised 1.041s/substitution including setup; peak 106.3 MiB. Final activity identical; all recorded metrics agree to 1e-13 relative tolerance.
- N=40, symmetric, chi=0.1: 7.370s → 2.667s; 2.76×; optimised 0.889s/substitution including setup; peak 106.6 MiB. Final activity identical; all recorded metrics agree to 1e-13 relative tolerance.
- N=120, loss_biased, chi=0.325: 16.103s → 3.053s; 5.27×; optimised 1.018s/substitution including setup; peak 106.6 MiB. Final activity identical; all recorded metrics agree to 1e-13 relative tolerance.

## Complete 250-substitution workloads

Timed without profiling, including ancestry and setup. Other experiment workers were active; these are measured wall times on this machine, specific to this hardware and concurrent workload.

| N | Kernel | chi | Baseline (s) | Optimised (s) | Speedup | Optimised s/fixation |
|---|---|---|---|---|---|---|
| 40 | loss_biased | 0.0 | 47.797 | 4.640 | 10.30× | 0.01856 |
| 120 | symmetric | 1.0 | 118.191 | 13.932 | 8.48× | 0.05573 |

Both final activity matrices are identical to baseline; all recorded baseline history fields agree within 1e-13 relative and 1e-14 absolute tolerance. Cumulative process peak RSS was approximately 107 MiB. The short-case profile demonstrates the main cost in repeated scalar expansion fixation calls; the vectorised code batches those calls.
