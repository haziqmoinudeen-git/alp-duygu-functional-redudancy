"""Canonical fitness calculations; original mathematical operations preserved."""

from __future__ import annotations

import numpy as np

from .parameters import W_FUNCTION


def reaction_performance(R: np.ndarray, K: np.ndarray):
    return R / (K + R)


def lambda_vector(R: np.ndarray, K: np.ndarray):
    return W_FUNCTION * K / (np.clip(R, 1e-300, None) * (K + R))


def exact_s_for_events(
    R: np.ndarray, K: np.ndarray, support: np.ndarray, actual_deltas: np.ndarray
):
    old_R = R[support]
    old_P = old_R / (K[support] + old_R)
    new_R = old_R[None, :] + actual_deltas
    valid_R = np.all(new_R > 0, axis=1)
    s = np.full(len(actual_deltas), -1.0, dtype=float)
    if np.any(valid_R):
        new_Rv = new_R[valid_R]
        new_P = new_Rv / (K[support][None, :] + new_Rv)
        dlogW = np.sum(
            W_FUNCTION[support][None, :] * (np.log(new_P) - np.log(old_P[None, :])), axis=1
        )
        s[valid_R] = np.expm1(dlogW)
    return s


def exact_s_single_function(R: np.ndarray, K: np.ndarray, j: int, delta: float):
    old_R = R[j]
    new_R = old_R + delta
    if new_R <= 0:
        return -1.0
    old_P = old_R / (K[j] + old_R)
    new_P = new_R / (K[j] + new_R)
    return float(np.expm1(W_FUNCTION[j] * (np.log(new_P) - np.log(old_P))))
