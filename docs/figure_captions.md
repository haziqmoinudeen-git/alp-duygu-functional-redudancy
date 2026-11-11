# Figure captions

1. Specialisation at 250 fixed substitutions. Panels show changes in mean effective breadth, strongest-function share, weakest-reaction performance, and mean effective reaction redundancy. Colour indicates chi; solid circles and dashed squares indicate loss-biased and symmetric kernels. Points are individual replicates, means and error bars are mean ± sample SD, five replicates per cell.
2. Catalytic output and evolutionary dynamics at 250 substitutions. A: endpoint mean total activity; the horizontal line is the ancestral 0.30. B: breadth, strongest share, and weakest performance divided by their ancestral values. C: absolute breadth. The reference is N=100, r0=12.5, loss-biased, chi=0.325, replicate 0, ancestry seed 100000007, evolution seed 1220000107. The history records every fixed substitution.
3. Same metrics and encodings as Figure 1 at 5N substitutions, five replicates per cell.
4. Same panels as Figure 2 for the 5N endpoint, with reference dynamics against substitutions per protein; N=100 reaches 500 substitutions at exposure 5.
5. Mean effective breadth versus fixed substitutions divided by N. Rows indicate mutation kernel, columns chi, line colours initial redundancy. Bands are ± sample SD, five seeds per complete cell. A common 201-point exposure grid uses last-observation-carried-forward alignment; no interpolation is used. The change version subtracts each trajectory's ancestor before averaging.

Alternative A shows cell-mean breadth changes as heatmaps, sharing one colour scale. B shows cell-mean specialisation and weakest-performance change at 250 substitutions; it implies no causal estimate. C pairs the two endpoints by the same condition and seeds; the diagonal denotes equality.

Export-specific replicate counts are generated in [figures/captions.md](../figures/captions.md). The primary experiment uses five per cell; extended exports use the actual paired replicate count.
