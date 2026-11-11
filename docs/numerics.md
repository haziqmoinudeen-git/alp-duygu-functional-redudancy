# Numerical methods

## Mutation catalogue

An ordinary mutation draws independent coordinate signs on a loop's latent support, conditional on at least one nonzero sign. One negative magnitude is shared by the negative coordinates and one positive magnitude by the positive coordinates. The catalogue probabilities include this conditioning and magnitude coupling. Loss to zero leaves latent support intact.

## Selection and fixation

Finite-effect selection coefficients are calculated from the change in log fitness and converted with `expm1`. Fixation uses exponential/logarithmic evaluation to avoid cancellation and overflow. The numerical neutral interval is $|s|<10^{-14}$. At exactly $s=\pm10^{-14}$, the strict sign masks return zero; this boundary convention is covered by a regression test. Fixed events are classified as beneficial or deleterious using a separate threshold of $\pm10^{-12}$.

## Event rates and generation time

For $L_{tot}$ loops, the origin rate assigned to each loop is $N_eU_{loop}/L_{tot}$, split between ordinary mutation and support expansion. Thus $N_eU_{loop}$ specifies the total origin intensity across the proteome. Generation-time interpretation uses this normalisation. The endpoint analyses compare fixed substitution counts and counts per protein.

## Cached hazards

After a substitution, loop hazards are refreshed for the selected protein and for loops whose latent supports overlap changed reactions. An expandable loop can also target reactions outside its current support. Its cached expansion contribution can therefore differ from a full recomputation after those reaction capacities change. Loop selection uses the cached total; candidate events on the selected loop are evaluated afresh.

A one-substitution diagnostic with $N=40$, loss-biased mutation, and $\chi=0.325$ finds a maximum absolute total-loop-hazard difference of $1.04\times10^{-14}$ and a maximum relative difference of $3.70\times10^{-6}$. The comparison involves 160 loops, of which 107 differ. These measurements characterise the cached procedure, which is not identical to a full-refresh hazard calculation. See [comparison data](../reports/hazard_cache_comparison.json); `uv run python scripts/check_hazards.py` reproduces this diagnostic.

## Trajectories and exposure

Histories contain each fixed substitution. Exposure curves align trajectories by taking the latest state with $n_{fix}/N\le x$ on a 201-point grid from 0 to 5. Means and sample-standard-deviation bands are calculated across replicates after alignment. The procedure does not interpolate between events.
