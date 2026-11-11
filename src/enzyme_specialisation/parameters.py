"""Model parameters and simulation configuration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

M_FUNCTIONS = 120


PROTEIN_FUNCTIONS = 15


LOOPS_PER_PROTEIN = 4


LOOP_SIZES = (5, 4, 4, 4)


MAX_LOOP_SUPPORT = 5


INITIAL_TOTAL_ACTIVITY = 0.3


MAX_ACTIVITY = 1.0


LOOP_STRUCTURAL_BUDGET = 1.0


LOSS_SIGN_PROBS = np.array([0.24, 0.65, 0.11], dtype=float)


MAGNITUDES = np.array([0.002, 0.01, 0.05, 0.2, 0.5], dtype=float)


POSITIVE_MAGNITUDE_PROBS = np.array([0.96, 0.035, 0.004, 0.0008, 0.0002], dtype=float)


NEGATIVE_MAGNITUDE_PROBS = np.array([0.1, 0.2, 0.4, 0.2, 0.1], dtype=float)


SYMMETRIC_SIGN_PROBS = np.array(
    [
        (LOSS_SIGN_PROBS[0] + LOSS_SIGN_PROBS[2]) / 2.0,
        LOSS_SIGN_PROBS[1],
        (LOSS_SIGN_PROBS[0] + LOSS_SIGN_PROBS[2]) / 2.0,
    ],
    dtype=float,
)


SYMMETRIC_MAGNITUDE_PROBS = (POSITIVE_MAGNITUDE_PROBS + NEGATIVE_MAGNITUDE_PROBS) / 2.0


SUPPORT_EXPANSION_PROB = 0.001


SUPPORT_EXPANSION_MAGNITUDES = MAGNITUDES.copy()


SUPPORT_EXPANSION_MAGNITUDE_PROBS = np.array([0.965, 0.03, 0.004, 0.0008, 0.0002], dtype=float)


EFFECTIVE_POP_SIZE = 1000000


U_LOOP = 1e-08


MAX_EFFECTIVE_GENERATIONS = 30000000000.0


ACTIVE_THRESHOLD = 0.005


K_REF = float(np.mean(LOOP_SIZES))


A0_MEAN = INITIAL_TOTAL_ACTIVITY / PROTEIN_FUNCTIONS


W_FUNCTION = np.ones(M_FUNCTIONS, dtype=float) / M_FUNCTIONS


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
    return 2.0 * chi / ((K_REF - 1.0) * A0_MEAN)


def redundancy_from_n(n_proteins: int) -> float:
    return n_proteins * PROTEIN_FUNCTIONS / M_FUNCTIONS


def n_from_redundancy(r0: float) -> int:
    n_float = r0 * M_FUNCTIONS / PROTEIN_FUNCTIONS
    n = int(round(n_float))
    if not np.isclose(redundancy_from_n(n), r0, atol=1e-10):
        raise ValueError(
            f"Requested redundancy {r0} does not correspond to an integer N for M={M_FUNCTIONS}, k_p={PROTEIN_FUNCTIONS}."
        )
    return n
