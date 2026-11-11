"""Canonical metrics calculations; original mathematical operations preserved."""

from __future__ import annotations

import numpy as np

from .fitness import reaction_performance
from .parameters import ACTIVE_THRESHOLD


def normalized_activity_profiles(A: np.ndarray):
    return A / np.clip(A.sum(axis=1, keepdims=True), 1e-300, None)


def effective_breadth(A: np.ndarray):
    p = normalized_activity_profiles(A)
    return 1.0 / np.clip(np.sum(p**2, axis=1), 1e-300, None)


def strongest_activity_share(A: np.ndarray):
    total = np.clip(A.sum(axis=1), 1e-300, None)
    return A.max(axis=1) / total


def effective_reaction_redundancy(A: np.ndarray, R: np.ndarray):
    shares = A / np.clip(R[None, :], 1e-300, None)
    return 1.0 / np.clip(np.sum(shares**2, axis=0), 1e-300, None)


def record_state(
    history: dict,
    time: float,
    fixed: int,
    A: np.ndarray,
    R: np.ndarray,
    K: np.ndarray,
    latent_support: np.ndarray,
):
    history["time"].append(float(time))
    history["fixed"].append(int(fixed))
    history["mean_breadth"].append(float(np.mean(effective_breadth(A))))
    history["mean_strongest_share"].append(float(np.mean(strongest_activity_share(A))))
    history["mean_total_activity"].append(float(np.mean(A.sum(axis=1))))
    history["weakest_P"].append(float(np.min(reaction_performance(R, K))))
    history["mean_redundancy"].append(float(np.mean(effective_reaction_redundancy(A, R))))
    history["mean_nonzero_functions"].append(float(np.mean(np.sum(A > 1e-15, axis=1))))
    history["mean_active_functions"].append(float(np.mean(np.sum(A >= ACTIVE_THRESHOLD, axis=1))))
    history["mean_latent_loop_support"].append(float(np.mean(latent_support.sum(axis=2))))
