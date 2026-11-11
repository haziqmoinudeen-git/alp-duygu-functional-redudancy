# Experiments

The factorial design crosses five initial redundancy conditions, two mutation kernels, and four structural incompatibility values. Five replicates in each of 40 parameter cells give 200 trajectories.

| Proteins, N | Initial redundancy | Equal-total endpoint | Equal-exposure endpoint |
|---|---|---|---|
| 40 | 5 | 250 | 200 |
| 60 | 7.5 | 250 | 300 |
| 80 | 10 | 250 | 400 |
| 100 | 12.5 | 250 | 500 |
| 120 | 15 | 250 | 600 |

Each history runs to the larger endpoint and records both states. Both comparisons use effective protein breadth, strongest-function share, weakest-reaction performance, effective reaction redundancy, and mean total catalytic activity. Ancestor-to-endpoint changes are computed within each replicate.

The exposure coordinate $x=n_{fix}/N$ compares proteome sizes at matched substitutions per protein. Absolute and ancestor-subtracted breadth curves extend to x=5 for every condition.

`configs/paper.yaml` specifies the five-replicate dataset; `configs/extended.yaml` specifies 20 replicates per cell. Seed assignment is deterministic and independent of parallel scheduling. The same redundancy/replicate ancestor is shared across kernel and incompatibility conditions.
