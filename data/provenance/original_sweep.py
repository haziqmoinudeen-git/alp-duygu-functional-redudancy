#!/usr/bin/env python3
"""
Integrated Reaction Theory V2
Fast origin-fixation simulation + parameter sweep

This implementation follows the model specified in:
    "A Minimal Model of Functional Redundancy As A Route For Enzyme
    Specialisation" (September 2026).

Implemented model features
--------------------------
Activity and fitness
    A_ij >= 0
    R_j = sum_i A_ij
    P_j = R_j / (K_j + R_j)
    W   = product_j P_j**w_j

Loop-local mutation
    A_ij = sum_l L_i,l,j
    Each loop has a LATENT mutational support T_i,l.
    Ordinary loss to zero does NOT remove a reaction from T_i,l.
    A later ordinary positive mutation can therefore restore it.

Structural feasibility
    C_i,l = sum_j L_i,l,j
            + eta * sum_{j<k} L_i,l,j L_i,l,k
    A proposed state is allowed only if:
        C_i,l <= B_loop
        0 <= A_ij <= a_max

Dimensionless incompatibility
    chi = eta * (k_ref - 1) * a0 / 2

Mutation kernels
    1. loss_biased:
         p(-) > p(+), with the original separate positive and negative
         magnitude distributions.
    2. symmetric:
         p(-) = p(+), with a single matched magnitude distribution.

Rare support expansion
    With probability rho_exp among mutation origins, a loop can add one
    new reaction to its latent support, subject to the same structural rule.

Heterogeneous ancestral generalists
    Each protein initially supports k_p reactions, but its total starting
    activity A_i,0 is distributed heterogeneously across those reactions.
    Positive Gamma weights are normalized to A_i,0. This is equivalent to
    a symmetric Dirichlet distribution after normalization.

Redundancy sweep
    r0 ~= N * k_p / M.
    The sweep changes N while holding M and k_p fixed.
    For every condition:
        K_j = R_j(0)
    so every reaction starts at P_j(0) = 1/2.

Experimental design
    redundancy x mutation kernel x structural incompatibility

Numerical evaluation
--------------------
- Finite-effect selection coefficients are calculated from fitness ratios.
- Event hazards combine mutation origin rates and fixation probabilities.
- Mutation catalogues are vectorized and cached.
- Only hazards that can change after a substitution are recomputed.
- Independent sweep conditions can run in parallel.
- The default "quick" preset stops after a modest equal number of fixed
  substitutions for speed.
- The "paper" preset uses a larger grid and more replicates.
- Set --max-fixed 0 to use only the generation horizon.

Dependencies
------------
numpy
pandas
matplotlib

Examples
--------
Quick sweep:
    python integrated_reaction_theory_v2_fast_sweep.py

Larger sweep:
    python integrated_reaction_theory_v2_fast_sweep.py --preset paper --jobs 6

Full 30-billion-generation horizon instead of fixed-substitution cap:
    python integrated_reaction_theory_v2_fast_sweep.py \
        --preset quick --max-fixed 0 --max-generations 3e10

Single reference condition:
    python integrated_reaction_theory_v2_fast_sweep.py \
        --single --redundancy 12.5 --kernel loss_biased --chi 0.325 \
        --max-fixed 300
"""

from __future__ import annotations

import argparse
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict, dataclass
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ============================================================
# MODEL PARAMETERS
# ============================================================

M_FUNCTIONS = 120
PROTEIN_FUNCTIONS = 15

LOOPS_PER_PROTEIN = 4
LOOP_SIZES = (5, 4, 4, 4)
MAX_LOOP_SUPPORT = 5

INITIAL_TOTAL_ACTIVITY = 0.30
MAX_ACTIVITY = 1.0

LOOP_STRUCTURAL_BUDGET = 1.0

# Loss-biased reference kernel.
LOSS_SIGN_PROBS = np.array([0.24, 0.65, 0.11], dtype=float)

MAGNITUDES = np.array(
    [0.002, 0.010, 0.050, 0.200, 0.500],
    dtype=float,
)

POSITIVE_MAGNITUDE_PROBS = np.array(
    [0.9600, 0.0350, 0.0040, 0.0008, 0.0002],
    dtype=float,
)

NEGATIVE_MAGNITUDE_PROBS = np.array(
    [0.10, 0.20, 0.40, 0.20, 0.10],
    dtype=float,
)

# Symmetric mutation condition:
# preserve p0 and the total probability of a non-zero sign, but split
# gain/loss equally. Use the midpoint of the two reference magnitude
# distributions so Q_- = Q_+ while preserving the overall effect-size scale.
SYMMETRIC_SIGN_PROBS = np.array(
    [
        (LOSS_SIGN_PROBS[0] + LOSS_SIGN_PROBS[2]) / 2.0,
        LOSS_SIGN_PROBS[1],
        (LOSS_SIGN_PROBS[0] + LOSS_SIGN_PROBS[2]) / 2.0,
    ],
    dtype=float,
)

SYMMETRIC_MAGNITUDE_PROBS = (
    POSITIVE_MAGNITUDE_PROBS + NEGATIVE_MAGNITUDE_PROBS
) / 2.0

SUPPORT_EXPANSION_PROB = 1e-3
SUPPORT_EXPANSION_MAGNITUDES = MAGNITUDES.copy()
SUPPORT_EXPANSION_MAGNITUDE_PROBS = np.array(
    [0.9650, 0.0300, 0.0040, 0.0008, 0.0002],
    dtype=float,
)

EFFECTIVE_POP_SIZE = 1_000_000
U_LOOP = 1e-8

MAX_EFFECTIVE_GENERATIONS = 30_000_000_000.0

ACTIVE_THRESHOLD = 0.005

K_REF = float(np.mean(LOOP_SIZES))
A0_MEAN = INITIAL_TOTAL_ACTIVITY / PROTEIN_FUNCTIONS

W_FUNCTION = np.ones(M_FUNCTIONS, dtype=float) / M_FUNCTIONS


# ============================================================
# CONFIG
# ============================================================

@dataclass(frozen=True)
class RunConfig:
    n_proteins: int
    mutation_kernel: str
    chi: float

    replicate: int
    ancestry_seed: int
    evolution_seed: int

    ancestor_alpha: float
    network_swaps: int

    effective_pop_size: int
    u_loop: float

    max_fixed: int
    max_generations: float
    record_every: int

    save_history: bool = False


def eta_from_chi(chi: float) -> float:
    """
    chi = eta * (k_ref - 1) * a0 / 2
    """
    if chi < 0:
        raise ValueError("chi must be >= 0.")

    if chi == 0:
        return 0.0

    return (
        2.0 * chi
        / ((K_REF - 1.0) * A0_MEAN)
    )


def redundancy_from_n(n_proteins: int) -> float:
    return (
        n_proteins
        * PROTEIN_FUNCTIONS
        / M_FUNCTIONS
    )


def n_from_redundancy(r0: float) -> int:
    n_float = (
        r0
        * M_FUNCTIONS
        / PROTEIN_FUNCTIONS
    )
    n = int(round(n_float))

    if not np.isclose(
        redundancy_from_n(n),
        r0,
        atol=1e-10,
    ):
        raise ValueError(
            f"Requested redundancy {r0} does not correspond to an integer N "
            f"for M={M_FUNCTIONS}, k_p={PROTEIN_FUNCTIONS}."
        )

    return n


# ============================================================
# MUTATION CATALOGUES
# ============================================================

@lru_cache(maxsize=None)
def mutation_catalogue(kernel: str, k: int):
    if kernel == "loss_biased":
        sign_probs = LOSS_SIGN_PROBS
        neg_probs = NEGATIVE_MAGNITUDE_PROBS
        pos_probs = POSITIVE_MAGNITUDE_PROBS

    elif kernel == "symmetric":
        sign_probs = SYMMETRIC_SIGN_PROBS
        neg_probs = SYMMETRIC_MAGNITUDE_PROBS
        pos_probs = SYMMETRIC_MAGNITUDE_PROBS

    else:
        raise ValueError(
            "mutation_kernel must be 'loss_biased' or 'symmetric'."
        )

    p_map = {
        -1: sign_probs[0],
         0: sign_probs[1],
         1: sign_probs[2],
    }

    all_zero_prob = sign_probs[1] ** k
    normalizer = 1.0 - all_zero_prob

    deltas = []
    probs = []
    n_negative = []
    n_positive = []

    for signs_tuple in product((-1, 0, 1), repeat=k):
        if not any(x != 0 for x in signs_tuple):
            continue

        signs = np.asarray(signs_tuple, dtype=int)

        sign_prob = 1.0
        for x in signs_tuple:
            sign_prob *= p_map[x]
        sign_prob /= normalizer

        has_negative = np.any(signs < 0)
        has_positive = np.any(signs > 0)

        negative_choices = (
            list(zip(MAGNITUDES, neg_probs))
            if has_negative
            else [(0.0, 1.0)]
        )

        positive_choices = (
            list(zip(MAGNITUDES, pos_probs))
            if has_positive
            else [(0.0, 1.0)]
        )

        for m_neg, p_neg in negative_choices:
            for m_pos, p_pos in positive_choices:
                delta = np.zeros(k, dtype=float)

                if has_negative:
                    delta[signs < 0] = -m_neg

                if has_positive:
                    delta[signs > 0] = +m_pos

                deltas.append(delta)
                probs.append(
                    sign_prob
                    * p_neg
                    * p_pos
                )
                n_negative.append(
                    int(np.sum(signs < 0))
                )
                n_positive.append(
                    int(np.sum(signs > 0))
                )

    deltas = np.asarray(deltas, dtype=float)
    probs = np.asarray(probs, dtype=float)

    probs /= probs.sum()

    return {
        "deltas": deltas,
        "probs": probs,
        "n_negative": np.asarray(
            n_negative,
            dtype=int,
        ),
        "n_positive": np.asarray(
            n_positive,
            dtype=int,
        ),
    }


# ============================================================
# STRUCTURAL COST
# ============================================================

def vectorized_loop_structural_cost(
    candidate_matrix: np.ndarray,
    eta: float,
):
    totals = candidate_matrix.sum(axis=1)

    pairwise = 0.5 * (
        totals**2
        - np.sum(
            candidate_matrix**2,
            axis=1,
        )
    )

    return (
        totals
        + eta * pairwise
    )


def loop_structural_cost(
    loop_values: np.ndarray,
    eta: float,
):
    total = float(
        np.sum(loop_values)
    )

    pairwise = 0.5 * (
        total**2
        - float(
            np.sum(loop_values**2)
        )
    )

    return (
        total
        + eta * pairwise
    )


# ============================================================
# RANDOM BALANCED BIPARTITE ANCESTRAL NETWORK
# ============================================================

def balanced_column_degrees(
    n_proteins: int,
):
    edges = (
        n_proteins
        * PROTEIN_FUNCTIONS
    )

    base = edges // M_FUNCTIONS
    remainder = edges % M_FUNCTIONS

    degrees = np.full(
        M_FUNCTIONS,
        base,
        dtype=int,
    )

    if remainder > 0:
        degrees[:remainder] += 1

    return degrees


def make_balanced_simple_network(
    rng: np.random.Generator,
    n_proteins: int,
    network_swaps: int,
):
    """
    Construct a simple bipartite graph with:
      row degree = PROTEIN_FUNCTIONS exactly
      balanced reaction-side degree sequence exactly

    Start from a cyclic simple graph, then randomize by degree-preserving
    double-edge swaps.
    """
    if M_FUNCTIONS < PROTEIN_FUNCTIONS:
        raise ValueError(
            "M_FUNCTIONS must be >= PROTEIN_FUNCTIONS."
        )

    graph = np.zeros(
        (
            n_proteins,
            M_FUNCTIONS,
        ),
        dtype=bool,
    )

    # Cyclic scaffold.
    for i in range(n_proteins):
        funcs = [
            (
                i * PROTEIN_FUNCTIONS
                + q
            )
            % M_FUNCTIONS
            for q
            in range(PROTEIN_FUNCTIONS)
        ]
        graph[i, funcs] = True

    # The cyclic construction produces balanced column degrees for the
    # redundancy values used in the sweep. Verify against target sequence.
    target_cols = np.sort(
        balanced_column_degrees(
            n_proteins
        )
    )

    observed_cols = np.sort(
        graph.sum(axis=0)
    )

    if not np.array_equal(
        observed_cols,
        target_cols,
    ):
        # Fall back to a greedy balanced simple graph if the cyclic scaffold
        # does not exactly match the target degree sequence.
        graph[:] = False

        remaining = (
            balanced_column_degrees(
                n_proteins
            )
            .astype(int)
        )

        for i in range(n_proteins):
            # Random tie-breaking among reactions with most remaining degree.
            jitter = rng.random(
                M_FUNCTIONS
            )
            order = np.lexsort(
                (jitter, -remaining)
            )

            chosen = order[
                :PROTEIN_FUNCTIONS
            ]

            if np.any(
                remaining[chosen]
                <= 0
            ):
                raise RuntimeError(
                    "Could not construct balanced ancestral network."
                )

            graph[i, chosen] = True
            remaining[chosen] -= 1

        if np.any(
            remaining != 0
        ):
            raise RuntimeError(
                "Balanced network construction left unmatched degrees."
            )

    row_before = (
        graph.sum(axis=1)
        .copy()
    )
    col_before = (
        graph.sum(axis=0)
        .copy()
    )

    neighbours = [
        set(
            np.flatnonzero(
                graph[i]
            ).tolist()
        )
        for i
        in range(n_proteins)
    ]

    successful = 0
    attempts = 0
    max_attempts = max(
        1000,
        network_swaps * 50,
    )

    while (
        successful
        < network_swaps
    ):
        attempts += 1

        if attempts > max_attempts:
            raise RuntimeError(
                "Could not complete requested degree-preserving edge swaps."
            )

        p1, p2 = rng.choice(
            n_proteins,
            size=2,
            replace=False,
        )

        f1 = int(
            rng.choice(
                tuple(
                    neighbours[p1]
                )
            )
        )
        f2 = int(
            rng.choice(
                tuple(
                    neighbours[p2]
                )
            )
        )

        if f1 == f2:
            continue

        if (
            f2 in neighbours[p1]
            or f1 in neighbours[p2]
        ):
            continue

        neighbours[p1].remove(f1)
        neighbours[p1].add(f2)

        neighbours[p2].remove(f2)
        neighbours[p2].add(f1)

        successful += 1

    randomized = np.zeros_like(
        graph
    )

    for i, funcs in enumerate(
        neighbours
    ):
        randomized[
            i,
            list(funcs),
        ] = True

    if not np.array_equal(
        randomized.sum(axis=1),
        row_before,
    ):
        raise RuntimeError(
            "Protein degrees changed during randomization."
        )

    if not np.array_equal(
        randomized.sum(axis=0),
        col_before,
    ):
        raise RuntimeError(
            "Reaction degrees changed during randomization."
        )

    return randomized


# ============================================================
# HETEROGENEOUS ANCESTRAL GENERALISTS
# ============================================================

def build_initial_architecture(
    cfg: RunConfig,
    eta: float,
):
    rng = np.random.default_rng(
        cfg.ancestry_seed
    )

    protein_function_mask = (
        make_balanced_simple_network(
            rng,
            cfg.n_proteins,
            cfg.network_swaps,
        )
    )

    L = np.zeros(
        (
            cfg.n_proteins,
            LOOPS_PER_PROTEIN,
            M_FUNCTIONS,
        ),
        dtype=float,
    )

    # Latent loop support T_i,l.
    latent_support = np.zeros_like(
        L,
        dtype=bool,
    )

    A_target = np.zeros(
        (
            cfg.n_proteins,
            M_FUNCTIONS,
        ),
        dtype=float,
    )

    for i in range(
        cfg.n_proteins
    ):
        funcs = rng.permutation(
            np.flatnonzero(
                protein_function_mask[i]
            )
        )

        loops = [
            funcs[0:5],
            funcs[5:9],
            funcs[9:13],
            np.array(
                [
                    funcs[13],
                    funcs[14],
                    funcs[0],
                    funcs[5],
                ],
                dtype=int,
            ),
        ]

        for ell, support in enumerate(
            loops
        ):
            latent_support[
                i,
                ell,
                support,
            ] = True

        # Positive heterogeneous weights.
        weights = rng.gamma(
            shape=cfg.ancestor_alpha,
            scale=1.0,
            size=PROTEIN_FUNCTIONS,
        )

        weights /= weights.sum()

        activities = (
            INITIAL_TOTAL_ACTIVITY
            * weights
        )

        A_target[
            i,
            funcs,
        ] = activities

        # Split protein-level activity equally across every loop whose
        # latent support contains that reaction.
        incidence = (
            latent_support[i]
            .sum(axis=0)
        )

        for ell, support in enumerate(
            loops
        ):
            for j in support:
                L[
                    i,
                    ell,
                    j,
                ] = (
                    A_target[i, j]
                    / incidence[j]
                )

    A = L.sum(axis=1)

    if not np.allclose(
        A,
        A_target,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise RuntimeError(
            "Loop decomposition failed to reproduce ancestral activities."
        )

    if not np.allclose(
        A.sum(axis=1),
        INITIAL_TOTAL_ACTIVITY,
        rtol=1e-10,
        atol=1e-12,
    ):
        raise RuntimeError(
            "Ancestral protein total activity is incorrect."
        )

    if not np.all(
        (
            A > 0
        ).sum(axis=1)
        == PROTEIN_FUNCTIONS
    ):
        raise RuntimeError(
            "Ancestral protein breadth is not k_p."
        )

    for i in range(
        cfg.n_proteins
    ):
        for ell in range(
            LOOPS_PER_PROTEIN
        ):
            support = np.flatnonzero(
                latent_support[
                    i,
                    ell,
                ]
            )

            c = loop_structural_cost(
                L[
                    i,
                    ell,
                    support,
                ],
                eta,
            )

            if (
                c
                > LOOP_STRUCTURAL_BUDGET
                + 1e-12
            ):
                raise RuntimeError(
                    f"Ancestral loop violates structural budget: "
                    f"C={c:.6g}, chi={cfg.chi}, eta={eta:.6g}. "
                    "Increase ancestor_alpha or reduce chi."
                )

    return (
        L,
        A,
        latent_support,
    )


# ============================================================
# FITNESS
# ============================================================

def reaction_performance(
    R: np.ndarray,
    K: np.ndarray,
):
    return (
        R
        / (K + R)
    )


def lambda_vector(
    R: np.ndarray,
    K: np.ndarray,
):
    return (
        W_FUNCTION
        * K
        / (
            np.clip(
                R,
                1e-300,
                None,
            )
            * (K + R)
        )
    )


def exact_s_for_events(
    R: np.ndarray,
    K: np.ndarray,
    support: np.ndarray,
    actual_deltas: np.ndarray,
):
    old_R = R[support]

    old_P = (
        old_R
        / (
            K[support]
            + old_R
        )
    )

    new_R = (
        old_R[None, :]
        + actual_deltas
    )

    valid_R = np.all(
        new_R > 0,
        axis=1,
    )

    s = np.full(
        len(actual_deltas),
        -1.0,
        dtype=float,
    )

    if np.any(
        valid_R
    ):
        new_Rv = new_R[
            valid_R
        ]

        new_P = (
            new_Rv
            / (
                K[support][None, :]
                + new_Rv
            )
        )

        dlogW = np.sum(
            W_FUNCTION[
                support
            ][None, :]
            * (
                np.log(new_P)
                - np.log(
                    old_P[None, :]
                )
            ),
            axis=1,
        )

        s[
            valid_R
        ] = np.expm1(
            dlogW
        )

    return s


def exact_s_single_function(
    R: np.ndarray,
    K: np.ndarray,
    j: int,
    delta: float,
):
    old_R = R[j]
    new_R = (
        old_R
        + delta
    )

    if new_R <= 0:
        return -1.0

    old_P = (
        old_R
        / (
            K[j]
            + old_R
        )
    )

    new_P = (
        new_R
        / (
            K[j]
            + new_R
        )
    )

    return float(
        np.expm1(
            W_FUNCTION[j]
            * (
                np.log(new_P)
                - np.log(old_P)
            )
        )
    )


# ============================================================
# FIXATION PROBABILITY
# ============================================================

def vectorized_fixation_probability(
    s: np.ndarray,
    Ne: int,
):
    s = np.asarray(
        s,
        dtype=float,
    )

    p = np.zeros_like(
        s
    )

    neutral = (
        np.abs(s)
        < 1e-14
    )
    positive = (
        s > 1e-14
    )
    negative = (
        s < -1e-14
    )

    p[
        neutral
    ] = 1.0 / Ne

    if np.any(
        positive
    ):
        sp = s[
            positive
        ]

        numerator = (
            -np.expm1(
                -2.0 * sp
            )
        )

        scaled = (
            2.0
            * Ne
            * sp
        )

        denominator = np.ones_like(
            scaled
        )

        small = (
            scaled
            <= 50.0
        )

        denominator[
            small
        ] = (
            -np.expm1(
                -scaled[small]
            )
        )

        p[
            positive
        ] = (
            numerator
            / denominator
        )

    if np.any(
        negative
    ):
        a = -s[
            negative
        ]

        x_num = (
            2.0 * a
        )
        x_den = (
            2.0
            * Ne
            * a
        )

        log_num = np.empty_like(
            x_num
        )
        log_den = np.empty_like(
            x_den
        )

        num_large = (
            x_num > 50.0
        )
        den_large = (
            x_den > 50.0
        )

        log_num[
            num_large
        ] = (
            x_num[
                num_large
            ]
            + np.log1p(
                -np.exp(
                    -x_num[
                        num_large
                    ]
                )
            )
        )

        log_num[
            ~num_large
        ] = np.log(
            np.expm1(
                x_num[
                    ~num_large
                ]
            )
        )

        log_den[
            den_large
        ] = (
            x_den[
                den_large
            ]
            + np.log1p(
                -np.exp(
                    -x_den[
                        den_large
                    ]
                )
            )
        )

        log_den[
            ~den_large
        ] = np.log(
            np.expm1(
                x_den[
                    ~den_large
                ]
            )
        )

        log_p = (
            log_num
            - log_den
        )

        safe = (
            log_p > -745
        )

        pneg = np.zeros_like(
            log_p
        )

        pneg[
            safe
        ] = np.exp(
            log_p[safe]
        )

        p[
            negative
        ] = pneg

    return np.clip(
        p,
        0.0,
        1.0,
    )


# ============================================================
# DIAGNOSTICS
# ============================================================

def normalized_activity_profiles(
    A: np.ndarray,
):
    return (
        A
        / np.clip(
            A.sum(
                axis=1,
                keepdims=True,
            ),
            1e-300,
            None,
        )
    )


def effective_breadth(
    A: np.ndarray,
):
    p = normalized_activity_profiles(
        A
    )

    return (
        1.0
        / np.clip(
            np.sum(
                p**2,
                axis=1,
            ),
            1e-300,
            None,
        )
    )


def strongest_activity_share(
    A: np.ndarray,
):
    total = np.clip(
        A.sum(axis=1),
        1e-300,
        None,
    )

    return (
        A.max(axis=1)
        / total
    )


def effective_reaction_redundancy(
    A: np.ndarray,
    R: np.ndarray,
):
    shares = (
        A
        / np.clip(
            R[None, :],
            1e-300,
            None,
        )
    )

    return (
        1.0
        / np.clip(
            np.sum(
                shares**2,
                axis=0,
            ),
            1e-300,
            None,
        )
    )


def record_state(
    history: dict,
    time: float,
    fixed: int,
    A: np.ndarray,
    R: np.ndarray,
    K: np.ndarray,
    latent_support: np.ndarray,
):
    history[
        "time"
    ].append(
        float(time)
    )

    history[
        "fixed"
    ].append(
        int(fixed)
    )

    history[
        "mean_breadth"
    ].append(
        float(
            np.mean(
                effective_breadth(
                    A
                )
            )
        )
    )

    history[
        "mean_strongest_share"
    ].append(
        float(
            np.mean(
                strongest_activity_share(
                    A
                )
            )
        )
    )

    history[
        "mean_total_activity"
    ].append(
        float(
            np.mean(
                A.sum(axis=1)
            )
        )
    )

    history[
        "weakest_P"
    ].append(
        float(
            np.min(
                reaction_performance(
                    R,
                    K,
                )
            )
        )
    )

    history[
        "mean_redundancy"
    ].append(
        float(
            np.mean(
                effective_reaction_redundancy(
                    A,
                    R,
                )
            )
        )
    )

    history[
        "mean_nonzero_functions"
    ].append(
        float(
            np.mean(
                np.sum(
                    A > 1e-15,
                    axis=1,
                )
            )
        )
    )

    history[
        "mean_active_functions"
    ].append(
        float(
            np.mean(
                np.sum(
                    A
                    >= ACTIVE_THRESHOLD,
                    axis=1,
                )
            )
        )
    )

    history[
        "mean_latent_loop_support"
    ].append(
        float(
            np.mean(
                latent_support.sum(
                    axis=2
                )
            )
        )
    )


# ============================================================
# ONE CONDITION
# ============================================================

def run_one_condition(
    cfg: RunConfig,
):
    eta = eta_from_chi(
        cfg.chi
    )

    (
        L,
        A,
        latent_support,
    ) = build_initial_architecture(
        cfg,
        eta,
    )

    R = A.sum(
        axis=0
    )

    if np.any(
        R <= 0
    ):
        raise RuntimeError(
            "At least one organismal reaction begins with zero capacity."
        )

    # Condition-specific normalization.
    K = R.copy()

    loop_ids = [
        (i, ell)
        for i
        in range(
            cfg.n_proteins
        )
        for ell
        in range(
            LOOPS_PER_PROTEIN
        )
    ]

    n_loops_total = len(
        loop_ids
    )

    function_to_loops = [
        set()
        for _
        in range(
            M_FUNCTIONS
        )
    ]

    def rebuild_function_to_loops():
        for s in function_to_loops:
            s.clear()

        for lid, (i, ell) in enumerate(
            loop_ids
        ):
            for j in np.flatnonzero(
                latent_support[
                    i,
                    ell,
                ]
            ):
                function_to_loops[
                    j
                ].add(
                    lid
                )

    rebuild_function_to_loops()

    origin_rate_per_loop = (
        cfg.effective_pop_size
        * cfg.u_loop
        / n_loops_total
    )

    ordinary_origin_rate_per_loop = (
        origin_rate_per_loop
        * (
            1.0
            - SUPPORT_EXPANSION_PROB
        )
    )

    expansion_origin_rate_per_loop = (
        origin_rate_per_loop
        * SUPPORT_EXPANSION_PROB
    )

    rng = np.random.default_rng(
        cfg.evolution_seed
    )

    def evaluate_loop_events(
        loop_id: int,
    ):
        i, ell = loop_ids[
            loop_id
        ]

        support = np.flatnonzero(
            latent_support[
                i,
                ell,
            ]
        )

        k = len(
            support
        )

        current_loop = L[
            i,
            ell,
            support,
        ]

        catalogue = mutation_catalogue(
            cfg.mutation_kernel,
            k,
        )

        raw_deltas = catalogue[
            "deltas"
        ]

        ordinary_probs = catalogue[
            "probs"
        ]

        candidate_loops = np.maximum(
            current_loop[
                None,
                :
            ]
            + raw_deltas,
            0.0,
        )

        actual_deltas = (
            candidate_loops
            - current_loop[
                None,
                :
            ]
        )

        state_changing = np.any(
            np.abs(
                actual_deltas
            ) > 1e-15,
            axis=1,
        )

        candidate_A = (
            A[
                i,
                support,
            ][
                None,
                :
            ]
            + actual_deltas
        )

        within_bounds = np.all(
            candidate_A
            <= MAX_ACTIVITY
            + 1e-15,
            axis=1,
        )

        structural_ok = (
            vectorized_loop_structural_cost(
                candidate_loops,
                eta,
            )
            <= LOOP_STRUCTURAL_BUDGET
            + 1e-15
        )

        valid = (
            state_changing
            & within_bounds
            & structural_ok
        )

        s_exact = np.full(
            len(
                raw_deltas
            ),
            -1.0,
            dtype=float,
        )

        pfix = np.zeros(
            len(
                raw_deltas
            ),
            dtype=float,
        )

        if np.any(
            valid
        ):
            s_exact[
                valid
            ] = exact_s_for_events(
                R,
                K,
                support,
                actual_deltas[
                    valid
                ],
            )

            pfix[
                valid
            ] = (
                vectorized_fixation_probability(
                    s_exact[
                        valid
                    ],
                    cfg.effective_pop_size,
                )
            )

        ordinary_hazards = (
            ordinary_origin_rate_per_loop
            * ordinary_probs
            * pfix
        )

        # Rare support expansion.
        expansion_functions = np.array(
            [],
            dtype=int,
        )
        expansion_magnitudes = np.array(
            [],
            dtype=float,
        )
        expansion_s = np.array(
            [],
            dtype=float,
        )
        expansion_hazards = np.array(
            [],
            dtype=float,
        )

        if (
            k
            < MAX_LOOP_SUPPORT
        ):
            candidates = np.flatnonzero(
                ~latent_support[
                    i,
                    ell,
                ]
            )

            if len(
                candidates
            ) > 0:
                functions_out = []
                magnitudes_out = []
                s_out = []
                hazards_out = []

                origin_per_candidate = (
                    expansion_origin_rate_per_loop
                    / len(
                        candidates
                    )
                )

                for j in candidates:
                    for (
                        magnitude,
                        magnitude_prob,
                    ) in zip(
                        SUPPORT_EXPANSION_MAGNITUDES,
                        SUPPORT_EXPANSION_MAGNITUDE_PROBS,
                    ):
                        expanded_values = np.append(
                            current_loop,
                            magnitude,
                        )

                        if (
                            loop_structural_cost(
                                expanded_values,
                                eta,
                            )
                            > LOOP_STRUCTURAL_BUDGET
                            + 1e-15
                        ):
                            continue

                        if (
                            A[
                                i,
                                j,
                            ]
                            + magnitude
                            > MAX_ACTIVITY
                            + 1e-15
                        ):
                            continue

                        s_candidate = (
                            exact_s_single_function(
                                R,
                                K,
                                int(j),
                                float(
                                    magnitude
                                ),
                            )
                        )

                        pfix_candidate = float(
                            vectorized_fixation_probability(
                                np.array(
                                    [
                                        s_candidate
                                    ],
                                    dtype=float,
                                ),
                                cfg.effective_pop_size,
                            )[0]
                        )

                        hazard = (
                            origin_per_candidate
                            * magnitude_prob
                            * pfix_candidate
                        )

                        functions_out.append(
                            int(j)
                        )
                        magnitudes_out.append(
                            float(
                                magnitude
                            )
                        )
                        s_out.append(
                            float(
                                s_candidate
                            )
                        )
                        hazards_out.append(
                            float(
                                hazard
                            )
                        )

                expansion_functions = np.asarray(
                    functions_out,
                    dtype=int,
                )
                expansion_magnitudes = np.asarray(
                    magnitudes_out,
                    dtype=float,
                )
                expansion_s = np.asarray(
                    s_out,
                    dtype=float,
                )
                expansion_hazards = np.asarray(
                    hazards_out,
                    dtype=float,
                )

        return {
            "support": support,
            "candidate_loops": candidate_loops,
            "actual_deltas": actual_deltas,
            "ordinary_hazards": ordinary_hazards,
            "ordinary_s": s_exact,
            "expansion_functions": expansion_functions,
            "expansion_magnitudes": expansion_magnitudes,
            "expansion_s": expansion_s,
            "expansion_hazards": expansion_hazards,
            "total_hazard": float(
                ordinary_hazards.sum()
                + expansion_hazards.sum()
            ),
        }

    history = {
        "time": [],
        "fixed": [],
        "mean_breadth": [],
        "mean_strongest_share": [],
        "mean_total_activity": [],
        "weakest_P": [],
        "mean_redundancy": [],
        "mean_nonzero_functions": [],
        "mean_active_functions": [],
        "mean_latent_loop_support": [],
    }

    record_state(
        history,
        0.0,
        0,
        A,
        R,
        K,
        latent_support,
    )

    initial_breadth = history[
        "mean_breadth"
    ][0]

    initial_strongest_share = history[
        "mean_strongest_share"
    ][0]

    initial_weakest_P = history[
        "weakest_P"
    ][0]

    initial_redundancy = history[
        "mean_redundancy"
    ][0]

    loop_hazards = np.zeros(
        n_loops_total,
        dtype=float,
    )

    for lid in range(
        n_loops_total
    ):
        loop_hazards[
            lid
        ] = evaluate_loop_events(
            lid
        )[
            "total_hazard"
        ]

    time = 0.0
    fixed = 0

    fixed_support_expansions = 0
    fixed_beneficial = 0
    fixed_neutralish = 0
    fixed_deleterious = 0

    stop_reason = (
        "generation horizon"
    )

    while (
        time
        < cfg.max_generations
    ):
        if (
            cfg.max_fixed > 0
            and fixed
            >= cfg.max_fixed
        ):
            stop_reason = (
                "fixed-substitution cap"
            )
            break

        total_rate = float(
            loop_hazards.sum()
        )

        if (
            total_rate <= 0
            or not np.isfinite(
                total_rate
            )
        ):
            stop_reason = (
                "no accessible fixation"
            )
            break

        dt = rng.exponential(
            1.0
            / total_rate
        )

        if (
            time + dt
            >= cfg.max_generations
        ):
            time = float(
                cfg.max_generations
            )
            stop_reason = (
                "generation horizon"
            )
            break

        time += dt

        chosen_loop_id = int(
            rng.choice(
                n_loops_total,
                p=(
                    loop_hazards
                    / total_rate
                ),
            )
        )

        ev = evaluate_loop_events(
            chosen_loop_id
        )

        ordinary_total = float(
            ev[
                "ordinary_hazards"
            ].sum()
        )

        expansion_total = float(
            ev[
                "expansion_hazards"
            ].sum()
        )

        event_total = (
            ordinary_total
            + expansion_total
        )

        if event_total <= 0:
            loop_hazards[
                chosen_loop_id
            ] = 0.0
            continue

        i, ell = loop_ids[
            chosen_loop_id
        ]

        is_expansion = (
            expansion_total > 0
            and rng.random()
            < (
                expansion_total
                / event_total
            )
        )

        affected_functions = set()

        if is_expansion:
            probs = (
                ev[
                    "expansion_hazards"
                ]
                / expansion_total
            )

            eidx = int(
                rng.choice(
                    len(
                        probs
                    ),
                    p=probs,
                )
            )

            j_new = int(
                ev[
                    "expansion_functions"
                ][eidx]
            )

            magnitude = float(
                ev[
                    "expansion_magnitudes"
                ][eidx]
            )

            s_exact = float(
                ev[
                    "expansion_s"
                ][eidx]
            )

            L[
                i,
                ell,
                j_new,
            ] += magnitude

            A[
                i,
                j_new,
            ] += magnitude

            R[
                j_new
            ] += magnitude

            latent_support[
                i,
                ell,
                j_new,
            ] = True

            function_to_loops[
                j_new
            ].add(
                chosen_loop_id
            )

            fixed_support_expansions += 1

            affected_functions.add(
                j_new
            )

        else:
            probs = (
                ev[
                    "ordinary_hazards"
                ]
                / ordinary_total
            )

            oidx = int(
                rng.choice(
                    len(
                        probs
                    ),
                    p=probs,
                )
            )

            support = ev[
                "support"
            ]

            new_loop_values = ev[
                "candidate_loops"
            ][oidx]

            delta = ev[
                "actual_deltas"
            ][oidx]

            s_exact = float(
                ev[
                    "ordinary_s"
                ][oidx]
            )

            # IMPORTANT V2 rule:
            # activity reaching zero does NOT alter latent_support.
            L[
                i,
                ell,
                support,
            ] = new_loop_values

            A[
                i,
                support,
            ] += delta

            R[
                support
            ] += delta

            A[
                i
            ] = np.clip(
                A[
                    i
                ],
                0.0,
                MAX_ACTIVITY,
            )

            L[
                i,
                ell
            ] = np.clip(
                L[
                    i,
                    ell
                ],
                0.0,
                None,
            )

            changed = support[
                np.abs(
                    delta
                ) > 1e-15
            ]

            affected_functions.update(
                changed.tolist()
            )

        fixed += 1

        if s_exact > 1e-12:
            fixed_beneficial += 1
        elif s_exact < -1e-12:
            fixed_deleterious += 1
        else:
            fixed_neutralish += 1

        # Refresh selected-protein and overlapping-support loop hazards.
        affected_loop_ids = {
            chosen_loop_id
        }

        # Any event changes protein i's A profile, which can affect expansion
        # feasibility of every loop on that protein.
        for other_ell in range(
            LOOPS_PER_PROTEIN
        ):
            affected_loop_ids.add(
                i
                * LOOPS_PER_PROTEIN
                + other_ell
            )

        # Any changed R_j alters selection coefficients of all loops whose
        # latent support contains j.
        for j in affected_functions:
            affected_loop_ids.update(
                function_to_loops[
                    j
                ]
            )

        for lid in affected_loop_ids:
            loop_hazards[
                lid
            ] = evaluate_loop_events(
                lid
            )[
                "total_hazard"
            ]

        if (
            fixed
            % cfg.record_every
            == 0
        ):
            record_state(
                history,
                time,
                fixed,
                A,
                R,
                K,
                latent_support,
            )

    if (
        history[
            "fixed"
        ][-1]
        != fixed
    ):
        record_state(
            history,
            time,
            fixed,
            A,
            R,
            K,
            latent_support,
        )

    for key in history:
        history[
            key
        ] = np.asarray(
            history[
                key
            ]
        )

    summary = {
        "n_proteins": cfg.n_proteins,
        "redundancy_r0": redundancy_from_n(
            cfg.n_proteins
        ),
        "mutation_kernel": cfg.mutation_kernel,
        "chi": cfg.chi,
        "eta": eta,
        "replicate": cfg.replicate,
        "ancestry_seed": cfg.ancestry_seed,
        "evolution_seed": cfg.evolution_seed,

        "ancestor_alpha": cfg.ancestor_alpha,

        "fixed_substitutions": fixed,
        "effective_generations": time,
        "stop_reason": stop_reason,

        "initial_mean_breadth": initial_breadth,
        "final_mean_breadth": float(
            history[
                "mean_breadth"
            ][-1]
        ),
        "delta_mean_breadth": float(
            history[
                "mean_breadth"
            ][-1]
            - initial_breadth
        ),

        "initial_strongest_share": initial_strongest_share,
        "final_strongest_share": float(
            history[
                "mean_strongest_share"
            ][-1]
        ),
        "delta_strongest_share": float(
            history[
                "mean_strongest_share"
            ][-1]
            - initial_strongest_share
        ),

        "initial_weakest_P": initial_weakest_P,
        "final_weakest_P": float(
            history[
                "weakest_P"
            ][-1]
        ),
        "delta_weakest_P": float(
            history[
                "weakest_P"
            ][-1]
            - initial_weakest_P
        ),

        "initial_mean_redundancy": initial_redundancy,
        "final_mean_redundancy": float(
            history[
                "mean_redundancy"
            ][-1]
        ),
        "delta_mean_redundancy": float(
            history[
                "mean_redundancy"
            ][-1]
            - initial_redundancy
        ),

        "initial_mean_total_activity": float(
            history[
                "mean_total_activity"
            ][0]
        ),
        "final_mean_total_activity": float(
            history[
                "mean_total_activity"
            ][-1]
        ),

        "final_mean_nonzero_functions": float(
            history[
                "mean_nonzero_functions"
            ][-1]
        ),

        "final_mean_active_functions": float(
            history[
                "mean_active_functions"
            ][-1]
        ),

        "final_mean_latent_loop_support": float(
            history[
                "mean_latent_loop_support"
            ][-1]
        ),

        "fixed_support_expansions": fixed_support_expansions,
        "fixed_beneficial": fixed_beneficial,
        "fixed_neutralish": fixed_neutralish,
        "fixed_deleterious": fixed_deleterious,
    }

    payload = {
        "summary": summary,
    }

    if cfg.save_history:
        payload[
            "history"
        ] = history

        payload[
            "A_final"
        ] = A

    return payload


# ============================================================
# SWEEP
# ============================================================

def build_jobs(
    redundancies,
    kernels,
    chis,
    replicates,
    base_seed,
    ancestor_alpha,
    network_swaps,
    effective_pop_size,
    u_loop,
    max_fixed,
    max_generations,
    record_every,
):
    jobs = []

    condition_index = 0

    for r0 in redundancies:
        n = n_from_redundancy(
            r0
        )

        for replicate in range(
            replicates
        ):
            # Same ancestor across kernel/chi controls within an r0 replicate.
            ancestry_seed = (
                base_seed
                + 1_000_000
                * n
                + replicate
            )

            for kernel in kernels:
                for chi in chis:
                    evolution_seed = (
                        base_seed
                        + 10_000_000
                        * condition_index
                        + 10_000
                        * replicate
                        + n
                    )

                    save_history = (
                        np.isclose(
                            r0,
                            12.5,
                        )
                        and kernel
                        == "loss_biased"
                        and np.isclose(
                            chi,
                            0.325,
                        )
                        and replicate == 0
                    )

                    jobs.append(
                        RunConfig(
                            n_proteins=n,
                            mutation_kernel=kernel,
                            chi=float(
                                chi
                            ),
                            replicate=replicate,
                            ancestry_seed=ancestry_seed,
                            evolution_seed=evolution_seed,
                            ancestor_alpha=ancestor_alpha,
                            network_swaps=network_swaps,
                            effective_pop_size=effective_pop_size,
                            u_loop=u_loop,
                            max_fixed=max_fixed,
                            max_generations=max_generations,
                            record_every=record_every,
                            save_history=save_history,
                        )
                    )

                    condition_index += 1

    return jobs


def summarize_sweep(
    df: pd.DataFrame,
):
    metrics = [
        "delta_mean_breadth",
        "delta_strongest_share",
        "delta_weakest_P",
        "delta_mean_redundancy",
        "final_mean_total_activity",
    ]

    group_cols = [
        "redundancy_r0",
        "mutation_kernel",
        "chi",
    ]

    agg_spec = {}

    for metric in metrics:
        agg_spec[
            f"{metric}_mean"
        ] = (
            metric,
            "mean",
        )
        agg_spec[
            f"{metric}_sd"
        ] = (
            metric,
            "std",
        )

    agg_spec[
        "n_replicates"
    ] = (
        "replicate",
        "count",
    )

    summary = (
        df
        .groupby(
            group_cols,
            as_index=False,
        )
        .agg(
            **agg_spec
        )
    )

    return summary


# ============================================================
# PLOTS
# ============================================================

def plot_response(
    summary: pd.DataFrame,
    kernel: str,
    y_mean: str,
    y_sd: str,
    ylabel: str,
    title: str,
    path: Path,
):
    subset = summary[
        summary[
            "mutation_kernel"
        ]
        == kernel
    ].copy()

    plt.figure(
        figsize=(8, 5)
    )

    for chi in sorted(
        subset[
            "chi"
        ].unique()
    ):
        d = subset[
            np.isclose(
                subset[
                    "chi"
                ],
                chi,
            )
        ].sort_values(
            "redundancy_r0"
        )

        yerr = d[
            y_sd
        ].fillna(
            0.0
        )

        plt.errorbar(
            d[
                "redundancy_r0"
            ],
            d[
                y_mean
            ],
            yerr=yerr,
            marker="o",
            capsize=3,
            label=f"chi={chi:g}",
        )

    plt.axhline(
        0.0,
        linewidth=1,
        linestyle="--",
    )

    plt.xlabel(
        "Initial mean reaction redundancy r0"
    )
    plt.ylabel(
        ylabel
    )
    plt.title(
        title
    )
    plt.legend()
    plt.tight_layout()
    plt.savefig(
        path,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def plot_reference_trajectory(
    history: dict,
    outdir: Path,
):
    plots = [
        (
            "mean_breadth",
            "Mean effective functional breadth",
            "Reference trajectory: effective breadth",
            "reference_breadth.png",
        ),
        (
            "mean_strongest_share",
            "Mean strongest-function activity share",
            "Reference trajectory: catalytic concentration",
            "reference_strongest_share.png",
        ),
        (
            "weakest_P",
            "Weakest reaction performance",
            "Reference trajectory: weakest reaction performance",
            "reference_weakest_P.png",
        ),
        (
            "mean_redundancy",
            "Mean effective reaction redundancy",
            "Reference trajectory: distributed redundancy",
            "reference_redundancy.png",
        ),
        (
            "mean_total_activity",
            "Mean total activity per protein",
            "Reference trajectory: total activity",
            "reference_total_activity.png",
        ),
    ]

    for (
        key,
        ylabel,
        title,
        filename,
    ) in plots:
        plt.figure(
            figsize=(8, 5)
        )

        plt.plot(
            history[
                "fixed"
            ],
            history[
                key
            ],
        )

        plt.xlabel(
            "Fixed substitutions"
        )
        plt.ylabel(
            ylabel
        )
        plt.title(
            title
        )
        plt.tight_layout()
        plt.savefig(
            outdir
            / filename,
            dpi=300,
            bbox_inches="tight",
        )
        plt.close()


def plot_final_matrix(
    A_final: np.ndarray,
    outdir: Path,
):
    plt.figure(
        figsize=(10, 6)
    )

    plt.imshow(
        A_final,
        aspect="auto",
        interpolation="nearest",
    )

    plt.xlabel(
        "Organismal reaction j"
    )
    plt.ylabel(
        "Protein i"
    )
    plt.title(
        "Reference condition final activity matrix"
    )

    plt.colorbar(
        label="activity a_ij"
    )

    plt.tight_layout()
    plt.savefig(
        outdir
        / "reference_final_activity_matrix.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()


def make_all_plots(
    sweep_summary: pd.DataFrame,
    outdir: Path,
):
    metric_specs = [
        (
            "delta_mean_breadth",
            "Change in mean effective breadth",
            "Change in effective breadth",
            "breadth",
        ),
        (
            "delta_strongest_share",
            "Change in strongest-function activity share",
            "Change in catalytic concentration",
            "strongest_share",
        ),
        (
            "delta_weakest_P",
            "Change in weakest reaction performance",
            "Change in weakest-reaction performance",
            "weakest_P",
        ),
        (
            "delta_mean_redundancy",
            "Change in mean effective reaction redundancy",
            "Change in distributed reaction redundancy",
            "redundancy",
        ),
    ]

    for kernel in [
        "loss_biased",
        "symmetric",
    ]:
        for (
            metric,
            ylabel,
            title,
            stem,
        ) in metric_specs:
            plot_response(
                sweep_summary,
                kernel,
                f"{metric}_mean",
                f"{metric}_sd",
                ylabel,
                f"{title}: {kernel.replace('_', ' ')} kernel",
                outdir
                / (
                    f"{stem}_{kernel}.png"
                ),
            )


# ============================================================
# CLI
# ============================================================

def parse_csv_floats(
    value: str,
):
    return [
        float(
            x.strip()
        )
        for x
        in value.split(",")
        if x.strip()
    ]


def parse_csv_strings(
    value: str,
):
    return [
        x.strip()
        for x
        in value.split(",")
        if x.strip()
    ]


def parse_args():
    p = argparse.ArgumentParser(
        description=(
            "Fast parameter sweep for Integrated Reaction Theory V2."
        )
    )

    p.add_argument(
        "--outdir",
        default="irt_v2_sweep_results",
    )

    p.add_argument(
        "--preset",
        choices=[
            "quick",
            "paper",
        ],
        default="quick",
    )

    p.add_argument(
        "--jobs",
        type=int,
        default=max(
            1,
            min(
                4,
                (
                    os.cpu_count()
                    or 2
                )
                - 1,
            ),
        ),
    )

    p.add_argument(
        "--base-seed",
        type=int,
        default=7,
    )

    p.add_argument(
        "--ancestor-alpha",
        type=float,
        default=8.0,
        help=(
            "Gamma/Dirichlet concentration for heterogeneous ancestral "
            "activity weights. Larger values are closer to equal activities."
        ),
    )

    p.add_argument(
        "--network-swaps",
        type=int,
        default=None,
    )

    p.add_argument(
        "--replicates",
        type=int,
        default=None,
    )

    p.add_argument(
        "--redundancies",
        type=str,
        default=None,
        help=(
            "Comma-separated r0 values. Must correspond to integer N."
        ),
    )

    p.add_argument(
        "--chis",
        type=str,
        default=None,
    )

    p.add_argument(
        "--kernels",
        type=str,
        default=None,
        help="Comma-separated: loss_biased,symmetric",
    )

    p.add_argument(
        "--max-fixed",
        type=int,
        default=None,
        help=(
            "0 disables fixed-substitution stopping and uses only "
            "--max-generations."
        ),
    )

    p.add_argument(
        "--max-generations",
        type=float,
        default=MAX_EFFECTIVE_GENERATIONS,
    )

    p.add_argument(
        "--record-every",
        type=int,
        default=10,
    )

    p.add_argument(
        "--Ne",
        type=int,
        default=EFFECTIVE_POP_SIZE,
    )

    p.add_argument(
        "--u-loop",
        type=float,
        default=U_LOOP,
    )

    # Single-condition mode.
    p.add_argument(
        "--single",
        action="store_true",
    )

    p.add_argument(
        "--redundancy",
        type=float,
        default=12.5,
    )

    p.add_argument(
        "--kernel",
        choices=[
            "loss_biased",
            "symmetric",
        ],
        default="loss_biased",
    )

    p.add_argument(
        "--chi",
        type=float,
        default=0.325,
    )

    return p.parse_args()


def resolve_preset(
    args,
):
    if args.preset == "quick":
        defaults = {
            "redundancies": [
                5.0,
                10.0,
                12.5,
            ],
            "chis": [
                0.0,
                0.325,
                1.0,
            ],
            "kernels": [
                "loss_biased",
                "symmetric",
            ],
            "replicates": 1,
            "max_fixed": 40,
            "network_swaps": 5_000,
        }

    else:
        defaults = {
            "redundancies": [
                5.0,
                7.5,
                10.0,
                12.5,
                15.0,
            ],
            "chis": [
                0.0,
                0.1,
                0.325,
                1.0,
            ],
            "kernels": [
                "loss_biased",
                "symmetric",
            ],
            "replicates": 5,
            "max_fixed": 250,
            "network_swaps": 50_000,
        }

    redundancies = (
        parse_csv_floats(
            args.redundancies
        )
        if args.redundancies
        else defaults[
            "redundancies"
        ]
    )

    chis = (
        parse_csv_floats(
            args.chis
        )
        if args.chis
        else defaults[
            "chis"
        ]
    )

    kernels = (
        parse_csv_strings(
            args.kernels
        )
        if args.kernels
        else defaults[
            "kernels"
        ]
    )

    for kernel in kernels:
        if kernel not in {
            "loss_biased",
            "symmetric",
        }:
            raise ValueError(
                f"Unknown kernel: {kernel}"
            )

    replicates = (
        args.replicates
        if args.replicates
        is not None
        else defaults[
            "replicates"
        ]
    )

    max_fixed = (
        args.max_fixed
        if args.max_fixed
        is not None
        else defaults[
            "max_fixed"
        ]
    )

    network_swaps = (
        args.network_swaps
        if args.network_swaps
        is not None
        else defaults[
            "network_swaps"
        ]
    )

    return {
        "redundancies": redundancies,
        "chis": chis,
        "kernels": kernels,
        "replicates": replicates,
        "max_fixed": max_fixed,
        "network_swaps": network_swaps,
    }


# ============================================================
# MAIN
# ============================================================

def main():
    args = parse_args()

    outdir = Path(
        args.outdir
    )

    outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    settings = resolve_preset(
        args
    )

    if args.single:
        n = n_from_redundancy(
            args.redundancy
        )

        cfg = RunConfig(
            n_proteins=n,
            mutation_kernel=args.kernel,
            chi=args.chi,
            replicate=0,
            ancestry_seed=args.base_seed,
            evolution_seed=(
                args.base_seed
                + 100_000
            ),
            ancestor_alpha=args.ancestor_alpha,
            network_swaps=settings[
                "network_swaps"
            ],
            effective_pop_size=args.Ne,
            u_loop=args.u_loop,
            max_fixed=settings[
                "max_fixed"
            ],
            max_generations=args.max_generations,
            record_every=args.record_every,
            save_history=True,
        )

        print(
            "Running single condition:"
        )
        print(
            json.dumps(
                asdict(
                    cfg
                ),
                indent=2,
            )
        )

        result = run_one_condition(
            cfg
        )

        pd.DataFrame(
            [
                result[
                    "summary"
                ]
            ]
        ).to_csv(
            outdir
            / "single_summary.csv",
            index=False,
        )

        history = result[
            "history"
        ]

        pd.DataFrame(
            history
        ).to_csv(
            outdir
            / "single_history.csv",
            index=False,
        )

        plot_reference_trajectory(
            history,
            outdir,
        )

        plot_final_matrix(
            result[
                "A_final"
            ],
            outdir,
        )

        print(
            "\nCompleted."
        )
        print(
            f"Results: {outdir.resolve()}"
        )

        return

    jobs = build_jobs(
        redundancies=settings[
            "redundancies"
        ],
        kernels=settings[
            "kernels"
        ],
        chis=settings[
            "chis"
        ],
        replicates=settings[
            "replicates"
        ],
        base_seed=args.base_seed,
        ancestor_alpha=args.ancestor_alpha,
        network_swaps=settings[
            "network_swaps"
        ],
        effective_pop_size=args.Ne,
        u_loop=args.u_loop,
        max_fixed=settings[
            "max_fixed"
        ],
        max_generations=args.max_generations,
        record_every=args.record_every,
    )

    run_manifest = {
        "preset": args.preset,
        "jobs": args.jobs,
        "redundancies": settings[
            "redundancies"
        ],
        "chis": settings[
            "chis"
        ],
        "kernels": settings[
            "kernels"
        ],
        "replicates": settings[
            "replicates"
        ],
        "max_fixed": settings[
            "max_fixed"
        ],
        "max_generations": args.max_generations,
        "ancestor_alpha": args.ancestor_alpha,
        "network_swaps": settings[
            "network_swaps"
        ],
        "Ne": args.Ne,
        "u_loop": args.u_loop,
        "n_conditions": len(
            jobs
        ),
        "reference_chi_for_eta10": (
            10.0
            * (
                K_REF - 1.0
            )
            * A0_MEAN
            / 2.0
        ),
    }

    with open(
        outdir
        / "run_manifest.json",
        "w",
    ) as f:
        json.dump(
            run_manifest,
            f,
            indent=2,
        )

    print(
        "\nIntegrated Reaction Theory V2 sweep"
    )
    print(
        "==================================="
    )
    print(
        json.dumps(
            run_manifest,
            indent=2,
        )
    )
    print()

    results = []
    reference_payload = None

    if args.jobs <= 1:
        for idx, cfg in enumerate(
            jobs,
            start=1,
        ):
            print(
                f"[{idx}/{len(jobs)}] "
                f"r0={redundancy_from_n(cfg.n_proteins):g}, "
                f"kernel={cfg.mutation_kernel}, "
                f"chi={cfg.chi:g}, "
                f"rep={cfg.replicate}",
                flush=True,
            )

            payload = run_one_condition(
                cfg
            )

            results.append(
                payload[
                    "summary"
                ]
            )

            if cfg.save_history:
                reference_payload = payload

    else:
        with ProcessPoolExecutor(
            max_workers=args.jobs
        ) as executor:
            future_to_cfg = {
                executor.submit(
                    run_one_condition,
                    cfg,
                ): cfg
                for cfg
                in jobs
            }

            completed = 0

            for future in as_completed(
                future_to_cfg
            ):
                cfg = future_to_cfg[
                    future
                ]

                payload = future.result()

                completed += 1

                results.append(
                    payload[
                        "summary"
                    ]
                )

                if cfg.save_history:
                    reference_payload = payload

                print(
                    f"[{completed}/{len(jobs)}] "
                    f"r0={redundancy_from_n(cfg.n_proteins):g}, "
                    f"kernel={cfg.mutation_kernel}, "
                    f"chi={cfg.chi:g}, "
                    f"rep={cfg.replicate} finished",
                    flush=True,
                )

    replicate_df = pd.DataFrame(
        results
    ).sort_values(
        [
            "mutation_kernel",
            "chi",
            "redundancy_r0",
            "replicate",
        ]
    )

    replicate_df.to_csv(
        outdir
        / "replicate_results.csv",
        index=False,
    )

    sweep_summary = summarize_sweep(
        replicate_df
    )

    sweep_summary.to_csv(
        outdir
        / "sweep_summary.csv",
        index=False,
    )

    make_all_plots(
        sweep_summary,
        outdir,
    )

    if reference_payload is not None:
        history = reference_payload[
            "history"
        ]

        pd.DataFrame(
            history
        ).to_csv(
            outdir
            / "reference_history.csv",
            index=False,
        )

        plot_reference_trajectory(
            history,
            outdir,
        )

        plot_final_matrix(
            reference_payload[
                "A_final"
            ],
            outdir,
        )

    print(
        "\nCompleted."
    )
    print(
        f"Results saved to: {outdir.resolve()}"
    )
    print(
        "Primary files:"
    )
    print(
        "  replicate_results.csv"
    )
    print(
        "  sweep_summary.csv"
    )
    print(
        "  run_manifest.json"
    )
    print(
        "  response-surface PNG plots"
    )


if __name__ == "__main__":
    main()
