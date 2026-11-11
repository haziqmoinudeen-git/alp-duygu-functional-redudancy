# Simulation dataset

| Quantity | Count |
|---|---:|
| Factorial parameter cells | 40 |
| Replicates per cell | 5 |
| Trajectories | 200 |
| Fixed substitutions across shared histories | 82,000 |
| Endpoint records | 400 |

`results/paper/equal_total_250/` contains outcomes at 250 fixed substitutions. `results/paper/equal_exposure_c5/` contains outcomes at five substitutions per protein. Each directory contains a replicate-level CSV and a 40-cell summary table. The 200 compressed histories in `results/paper/trajectories/` include every event and final activity matrices.

Ancestry and evolution seeds are listed in each endpoint table and batch manifest. The batch manifests also record the RNG, environment, model parameters, and software versions. Runtime fields describe individual trajectory invocations; resumed-run values exclude earlier invocations. They are not whole-batch wall times.

The dataset uses five replicates per cell. The 20-replicate configuration is available for extended experiments. Publication figures include the two endpoint grids, catalytic-output and reference-dynamics panels, absolute and change-in-breadth exposure curves, and three alternative comparisons.
