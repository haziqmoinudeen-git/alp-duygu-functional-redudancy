# Computational performance

Expansion candidates are evaluated in contiguous arrays, with structural costs calculated once per magnitude per loop. Exact finite-effect selection and fixation calculations are vectorised. Mutation catalogues are cached by support size and mutation condition; loop dependencies determine hazard refreshes. Independent trajectories run in separate processes, subject to CPU and memory limits.

[Benchmark measurements](../reports/performance_benchmark.md) compare the vectorised and scalar reference implementations on matched workloads. The [raw timings](../reports/benchmark.json) and [profile](../reports/profile_baseline.txt) provide setup and memory details. Short profiled runs include initialisation and catalogue construction. Complete 250-substitution workloads were also timed separately.

Sharing trajectories between the two endpoints requires 82,000 fixed events, compared with 130,000 for separate five-replicate endpoint batches. Checkpointing and validated result reuse support interrupted and extended experiments.
