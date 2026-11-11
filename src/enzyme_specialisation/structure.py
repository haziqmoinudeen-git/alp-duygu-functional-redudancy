"""Loop structural costs and feasibility calculations."""

from __future__ import annotations

import numpy as np


def vectorized_loop_structural_cost(candidate_matrix: np.ndarray, eta: float):
    totals = candidate_matrix.sum(axis=1)
    pairwise = 0.5 * (totals**2 - np.sum(candidate_matrix**2, axis=1))
    return totals + eta * pairwise


def loop_structural_cost(loop_values: np.ndarray, eta: float):
    total = float(np.sum(loop_values))
    pairwise = 0.5 * (total**2 - float(np.sum(loop_values**2)))
    return total + eta * pairwise
