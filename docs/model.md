# Scientific model

The canonical source is `data/provenance/original_sweep.py`, checked against manuscript Methods II A–F. The manuscript is titled *A Minimal Model of Functional Redundancy as a Route to Enzyme Specialisation*, dated September 2026.

For equal protein abundance,

$$R_j=\sum_i a_{ij},\qquad P_j=\frac{R_j}{K_j+R_j},\qquad W=\prod_j P_j^{w_j},\quad w_j=1/M.$$

Activities sum over four loops. Ordinary mutations affect latent loop support; loss to zero never removes support. Each coordinate receives a sign, while a common negative magnitude and a common positive magnitude are drawn for the respective signed coordinates in the enumerated catalogue. The all-zero sign vector is excluded and the catalogue renormalised, exactly as in the supplied implementation.

Loop feasibility is $C=\sum_j L_j+\eta\sum_{j<k}L_jL_k\le1$, with $0\le a_{ij}\le1$. This constraint does not enter fitness directly. $\chi=\eta(k_{ref}-1)a_0/2$, $k_{ref}=4.25$, $a_0=0.02$.

The exact finite-effect selection coefficient is $s=W'/W-1$. Fixation is $1/N_e$ at neutrality and $(1-e^{-2s})/(1-e^{-2N_es})$ otherwise. The original numerical neutral threshold $|s|<10^{-14}$, exponent cutoffs, and underflow handling are retained. No alternate fixation approximation was introduced.

| Quantity | Preserved value |
|---|---|
| Reactions; functions per ancestor | 120; 15 |
| Loops; support sizes | 4; (5,4,4,4) |
| Total ancestral activity; Gamma shape | 0.30; 8 |
| Structural budget; individual activity cap | 1; 1 |
| $r_0$ | 5, 7.5, 10, 12.5, 15 |
| Proteins | 40, 60, 80, 100, 120 |
| $\chi$ | 0, 0.1, 0.325, 1 |
| Loss-biased signs (−,0,+) | (0.24,0.65,0.11) |
| Symmetric signs | (0.175,0.65,0.175) |
| Magnitudes | (0.002,0.010,0.050,0.200,0.500) |
| Positive probabilities | (0.9600,0.0350,0.0040,0.0008,0.0002) |
| Negative probabilities | (0.10,0.20,0.40,0.20,0.10) |
| Symmetric magnitude probabilities | elementwise midpoint of positive and negative |
| Expansion probability | 0.001 |
| Expansion magnitude probabilities | (0.9650,0.0300,0.0040,0.0008,0.0002) |
| Effective population; mutation rate | $10^6$; $10^{-8}$ |
| Degree-preserving swaps | 50,000 |
| Generation safety horizon | $3\times10^{10}$ |
| Active-function threshold | 0.005 |

The ancestor has fixed protein degree and balanced reaction degrees, with normalised Gamma weights. Loop decomposition splits overlap contributions; $K_j=R_j(0)$ gives $P_j(0)=0.5$. Proteins and the organismal reaction universe do not duplicate or grow.

Effective breadth is inverse Simpson concentration of within-protein activity shares. Strongest share is the maximum share; effective reaction redundancy is inverse Simpson concentration of between-protein contributions to each reaction. The plotted breadth, share, redundancy, and activity are arithmetic means over proteins or reactions as appropriate; weakest performance is a minimum over reactions.

The symmetric kernel is a mutation condition, **not a neutral evolutionary control**. See the [scientific audit](../reports/scientific_fidelity.md) for implementation discrepancies preserved for fidelity.
