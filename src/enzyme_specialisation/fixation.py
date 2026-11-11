"""Numerically stable fixation probabilities for origin–fixation dynamics."""

from __future__ import annotations

import numpy as np


def vectorized_fixation_probability(s: np.ndarray, Ne: int):
    s = np.asarray(s, dtype=float)
    p = np.zeros_like(s)
    neutral = np.abs(s) < 1e-14
    positive = s > 1e-14
    negative = s < -1e-14
    p[neutral] = 1.0 / Ne
    if np.any(positive):
        sp = s[positive]
        numerator = -np.expm1(-2.0 * sp)
        scaled = 2.0 * Ne * sp
        denominator = np.ones_like(scaled)
        small = scaled <= 50.0
        denominator[small] = -np.expm1(-scaled[small])
        p[positive] = numerator / denominator
    if np.any(negative):
        a = -s[negative]
        x_num = 2.0 * a
        x_den = 2.0 * Ne * a
        log_num = np.empty_like(x_num)
        log_den = np.empty_like(x_den)
        num_large = x_num > 50.0
        den_large = x_den > 50.0
        log_num[num_large] = x_num[num_large] + np.log1p(-np.exp(-x_num[num_large]))
        log_num[~num_large] = np.log(np.expm1(x_num[~num_large]))
        log_den[den_large] = x_den[den_large] + np.log1p(-np.exp(-x_den[den_large]))
        log_den[~den_large] = np.log(np.expm1(x_den[~den_large]))
        log_p = log_num - log_den
        safe = log_p > -745
        pneg = np.zeros_like(log_p)
        pneg[safe] = np.exp(log_p[safe])
        p[negative] = pneg
    return np.clip(p, 0.0, 1.0)
