"""Canonical ancestry calculations; original mathematical operations preserved."""

from __future__ import annotations

import numpy as np

from .parameters import (
    INITIAL_TOTAL_ACTIVITY,
    LOOP_STRUCTURAL_BUDGET,
    LOOPS_PER_PROTEIN,
    M_FUNCTIONS,
    PROTEIN_FUNCTIONS,
    RunConfig,
)
from .structure import loop_structural_cost


def balanced_column_degrees(n_proteins: int):
    edges = n_proteins * PROTEIN_FUNCTIONS
    base = edges // M_FUNCTIONS
    remainder = edges % M_FUNCTIONS
    degrees = np.full(M_FUNCTIONS, base, dtype=int)
    if remainder > 0:
        degrees[:remainder] += 1
    return degrees


def make_balanced_simple_network(rng: np.random.Generator, n_proteins: int, network_swaps: int):
    """
    Construct a simple bipartite graph with:
      row degree = PROTEIN_FUNCTIONS exactly
      balanced reaction-side degree sequence exactly

    Start from a cyclic simple graph, then randomize by degree-preserving
    double-edge swaps.
    """
    if M_FUNCTIONS < PROTEIN_FUNCTIONS:
        raise ValueError("M_FUNCTIONS must be >= PROTEIN_FUNCTIONS.")
    graph = np.zeros((n_proteins, M_FUNCTIONS), dtype=bool)
    for i in range(n_proteins):
        funcs = [(i * PROTEIN_FUNCTIONS + q) % M_FUNCTIONS for q in range(PROTEIN_FUNCTIONS)]
        graph[i, funcs] = True
    target_cols = np.sort(balanced_column_degrees(n_proteins))
    observed_cols = np.sort(graph.sum(axis=0))
    if not np.array_equal(observed_cols, target_cols):
        graph[:] = False
        remaining = balanced_column_degrees(n_proteins).astype(int)
        for i in range(n_proteins):
            jitter = rng.random(M_FUNCTIONS)
            order = np.lexsort((jitter, -remaining))
            chosen = order[:PROTEIN_FUNCTIONS]
            if np.any(remaining[chosen] <= 0):
                raise RuntimeError("Could not construct balanced ancestral network.")
            graph[i, chosen] = True
            remaining[chosen] -= 1
        if np.any(remaining != 0):
            raise RuntimeError("Balanced network construction left unmatched degrees.")
    row_before = graph.sum(axis=1).copy()
    col_before = graph.sum(axis=0).copy()
    neighbours = [set(np.flatnonzero(graph[i]).tolist()) for i in range(n_proteins)]
    successful = 0
    attempts = 0
    max_attempts = max(1000, network_swaps * 50)
    while successful < network_swaps:
        attempts += 1
        if attempts > max_attempts:
            raise RuntimeError("Could not complete requested degree-preserving edge swaps.")
        p1, p2 = rng.choice(n_proteins, size=2, replace=False)
        f1 = int(rng.choice(tuple(neighbours[p1])))
        f2 = int(rng.choice(tuple(neighbours[p2])))
        if f1 == f2:
            continue
        if f2 in neighbours[p1] or f1 in neighbours[p2]:
            continue
        neighbours[p1].remove(f1)
        neighbours[p1].add(f2)
        neighbours[p2].remove(f2)
        neighbours[p2].add(f1)
        successful += 1
    randomized = np.zeros_like(graph)
    for i, funcs in enumerate(neighbours):
        randomized[i, list(funcs)] = True
    if not np.array_equal(randomized.sum(axis=1), row_before):
        raise RuntimeError("Protein degrees changed during randomization.")
    if not np.array_equal(randomized.sum(axis=0), col_before):
        raise RuntimeError("Reaction degrees changed during randomization.")
    return randomized


def build_initial_architecture(cfg: RunConfig, eta: float):
    rng = np.random.default_rng(cfg.ancestry_seed)
    protein_function_mask = make_balanced_simple_network(rng, cfg.n_proteins, cfg.network_swaps)
    L = np.zeros((cfg.n_proteins, LOOPS_PER_PROTEIN, M_FUNCTIONS), dtype=float)
    latent_support = np.zeros_like(L, dtype=bool)
    A_target = np.zeros((cfg.n_proteins, M_FUNCTIONS), dtype=float)
    for i in range(cfg.n_proteins):
        funcs = rng.permutation(np.flatnonzero(protein_function_mask[i]))
        loops = [
            funcs[0:5],
            funcs[5:9],
            funcs[9:13],
            np.array([funcs[13], funcs[14], funcs[0], funcs[5]], dtype=int),
        ]
        for ell, support in enumerate(loops):
            latent_support[i, ell, support] = True
        weights = rng.gamma(shape=cfg.ancestor_alpha, scale=1.0, size=PROTEIN_FUNCTIONS)
        weights /= weights.sum()
        activities = INITIAL_TOTAL_ACTIVITY * weights
        A_target[i, funcs] = activities
        incidence = latent_support[i].sum(axis=0)
        for ell, support in enumerate(loops):
            for j in support:
                L[i, ell, j] = A_target[i, j] / incidence[j]
    A = L.sum(axis=1)
    if not np.allclose(A, A_target, rtol=1e-10, atol=1e-12):
        raise RuntimeError("Loop decomposition failed to reproduce ancestral activities.")
    if not np.allclose(A.sum(axis=1), INITIAL_TOTAL_ACTIVITY, rtol=1e-10, atol=1e-12):
        raise RuntimeError("Ancestral protein total activity is incorrect.")
    if not np.all((A > 0).sum(axis=1) == PROTEIN_FUNCTIONS):
        raise RuntimeError("Ancestral protein breadth is not k_p.")
    for i in range(cfg.n_proteins):
        for ell in range(LOOPS_PER_PROTEIN):
            support = np.flatnonzero(latent_support[i, ell])
            c = loop_structural_cost(L[i, ell, support], eta)
            if c > LOOP_STRUCTURAL_BUDGET + 1e-12:
                raise RuntimeError(
                    f"Ancestral loop violates structural budget: C={c:.6g}, chi={cfg.chi}, eta={eta:.6g}. Increase ancestor_alpha or reduce chi."
                )
    return (L, A, latent_support)
