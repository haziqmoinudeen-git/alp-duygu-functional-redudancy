"""Finite loop-local mutation catalogues and event probabilities."""

from __future__ import annotations

from functools import lru_cache
from itertools import product

import numpy as np

from .parameters import (
    LOSS_SIGN_PROBS,
    MAGNITUDES,
    NEGATIVE_MAGNITUDE_PROBS,
    POSITIVE_MAGNITUDE_PROBS,
    SYMMETRIC_MAGNITUDE_PROBS,
    SYMMETRIC_SIGN_PROBS,
)


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
        raise ValueError("mutation_kernel must be 'loss_biased' or 'symmetric'.")
    p_map = {-1: sign_probs[0], 0: sign_probs[1], 1: sign_probs[2]}
    all_zero_prob = sign_probs[1] ** k
    normalizer = 1.0 - all_zero_prob
    deltas = []
    probs = []
    n_negative = []
    n_positive = []
    for signs_tuple in product((-1, 0, 1), repeat=k):
        if not any((x != 0 for x in signs_tuple)):
            continue
        signs = np.asarray(signs_tuple, dtype=int)
        sign_prob = 1.0
        for x in signs_tuple:
            sign_prob *= p_map[x]
        sign_prob /= normalizer
        has_negative = np.any(signs < 0)
        has_positive = np.any(signs > 0)
        negative_choices = list(zip(MAGNITUDES, neg_probs)) if has_negative else [(0.0, 1.0)]
        positive_choices = list(zip(MAGNITUDES, pos_probs)) if has_positive else [(0.0, 1.0)]
        for m_neg, p_neg in negative_choices:
            for m_pos, p_pos in positive_choices:
                delta = np.zeros(k, dtype=float)
                if has_negative:
                    delta[signs < 0] = -m_neg
                if has_positive:
                    delta[signs > 0] = +m_pos
                deltas.append(delta)
                probs.append(sign_prob * p_neg * p_pos)
                n_negative.append(int(np.sum(signs < 0)))
                n_positive.append(int(np.sum(signs > 0)))
    deltas = np.asarray(deltas, dtype=float)
    probs = np.asarray(probs, dtype=float)
    probs /= probs.sum()
    return {
        "deltas": deltas,
        "probs": probs,
        "n_negative": np.asarray(n_negative, dtype=int),
        "n_positive": np.asarray(n_positive, dtype=int),
    }
