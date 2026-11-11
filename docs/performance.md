# Performance

The baseline profile identifies scalar support-expansion feasibility and fixation calls as major repeated costs. The implementation batches expansion candidates in their original reaction/magnitude order, evaluates structural costs once per magnitude per loop, and vectorises exact selection and fixation calculations. Ordinary catalogues, event sampling, dependency refresh, and RNG consumption are unchanged. No Numba/GPU or new approximation is used.

See [measured benchmarks](../reports/performance_benchmark.md), [raw measurements](../reports/benchmark.json), and [baseline profile](../reports/profile_baseline.txt). Short benchmarks include initialisation, cached catalogue setup, and profiler overhead. They are not estimates of whole-batch acceleration. Trajectory reuse separately avoids duplicate runs for the two endpoints.
