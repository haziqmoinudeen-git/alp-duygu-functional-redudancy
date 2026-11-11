from dataclasses import replace
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from enzyme_specialisation import model as m
from enzyme_specialisation.experiments import registry


def config():
    return m.RunConfig(
        40, "loss_biased", 0.325, 0, 40000007, 47, 8.0, 100, 1000000, 1e-8, 3, 3e10, 1, True
    )


def test_ancestor():
    c = config()
    L, A, T = m.build_initial_architecture(c, m.eta_from_chi(c.chi))
    np.testing.assert_allclose(A.sum(1), 0.3, atol=1e-15)
    assert np.all((A > 0).sum(1) == 15)
    degrees = (A > 0).sum(0)
    assert degrees.max() - degrees.min() <= 1
    np.testing.assert_array_equal(T.sum(2), np.tile([5, 4, 4, 4], (40, 1)))
    np.testing.assert_allclose(m.reaction_performance(A.sum(0), A.sum(0)), 0.5)
    assert np.all(m.vectorized_loop_structural_cost(L.reshape(-1, 120), m.eta_from_chi(0.325)) <= 1)
    assert np.all((A >= 0) & (A <= 1))


def test_catalogues():
    for k in [4, 5]:
        for kernel in ["loss_biased", "symmetric"]:
            cat = m.mutation_catalogue(kernel, k)
            np.testing.assert_allclose(cat["probs"].sum(), 1, atol=1e-15)
            assert np.all(cat["probs"] >= 0)
            assert not np.any(np.all(cat["deltas"] == 0, axis=1))


def test_fixation():
    s = np.array([-0.01, -1e-8, 0.0, 1e-8, 0.01])
    p = m.vectorized_fixation_probability(s, 1000000)
    assert p[2] == 1e-6
    assert np.all(np.diff(p) >= 0)
    np.testing.assert_allclose(p[3:], -np.expm1(-2 * s[3:]) / -np.expm1(-2e6 * s[3:]), rtol=1e-14)


def test_exact_selection():
    R = np.ones(120)
    K = R.copy()
    d = np.array([[0.01], [-0.02]])
    a = m.exact_s_for_events(R, K, np.array([3]), d)
    b = [m.exact_s_single_function(R, K, 3, x) for x in d[:, 0]]
    np.testing.assert_array_equal(a, b)


def test_metrics():
    A = np.ones((4, 120))
    np.testing.assert_allclose(m.effective_breadth(A), 120)
    np.testing.assert_allclose(m.strongest_activity_share(A), 1 / 120)
    np.testing.assert_allclose(m.effective_reaction_redundancy(A, A.sum(0)), 4)


def test_reproducibility():
    c = config()
    a = m.run_one_condition(c)
    b = m.run_one_condition(c)
    np.testing.assert_array_equal(a["A_final"], b["A_final"])
    for key in a["history"]:
        np.testing.assert_array_equal(a["history"][key], b["history"][key])


def test_seed_registry():
    c = yaml.safe_load(Path("configs/paper.yaml").read_text())
    jobs = registry(c)
    df = pd.read_csv("data/reference/replicate_results.csv")
    assert len(jobs) == 200
    for j in jobs:
        r = df[
            (df.n_proteins == j.n_proteins)
            & (df.mutation_kernel == j.mutation_kernel)
            & (df.chi == j.chi)
            & (df.replicate == j.replicate)
        ].iloc[0]
        assert j.ancestry_seed == r.ancestry_seed and j.evolution_seed == r.evolution_seed
        assert j.max_fixed == max(250, 5 * j.n_proteins)


def test_checkpoint(tmp_path):
    # Simulate interruption after an atomic checkpoint, then continue identically.
    c = replace(config(), max_fixed=27)
    p = tmp_path / "checkpoint"
    import os

    original = os.replace

    def interrupt(a, b):
        original(a, b)
        raise InterruptedError("test interruption")

    os.replace = interrupt
    try:
        try:
            m.run_one_condition(c, p)
        except InterruptedError:
            pass
    finally:
        os.replace = original
    resumed = m.run_one_condition(c, p)
    fresh = m.run_one_condition(c)
    np.testing.assert_array_equal(resumed["A_final"], fresh["A_final"])
    for key in fresh["history"]:
        np.testing.assert_array_equal(resumed["history"][key], fresh["history"][key])


def test_finite_effect_fitness_ratio():
    rng = np.random.default_rng(19)
    R = rng.uniform(0.2, 1.0, 120)
    K = rng.uniform(0.1, 0.5, 120)
    support = np.array([2, 19, 50])
    delta = np.array([[0.01, -0.02, 0.03], [-0.02, 0.01, 0.0]])
    baseline = np.prod(m.reaction_performance(R, K) ** m.W_FUNCTION)
    expected = []
    for d in delta:
        new = R.copy()
        new[support] += d
        expected.append(np.prod(m.reaction_performance(new, K) ** m.W_FUNCTION) / baseline - 1)
    np.testing.assert_allclose(m.exact_s_for_events(R, K, support, delta), expected, atol=3e-15)
    assert m.exact_s_single_function(R, K, 2, -R[2]) == -1


def test_structural_vector_matches_scalar():
    values = np.random.default_rng(11).uniform(0.0, 0.5, (19, 5))
    for chi in [0.0, 0.1, 0.325, 1.0]:
        eta = m.eta_from_chi(chi)
        expected = [m.loop_structural_cost(x, eta) for x in values]
        np.testing.assert_allclose(
            m.vectorized_loop_structural_cost(values, eta), expected, rtol=1e-14
        )
    assert m.eta_from_chi(0.325) == 10.0


def test_magnitude_probability_preservation():
    for kernel, signs, positive, negative in [
        ("loss_biased", m.LOSS_SIGN_PROBS, m.POSITIVE_MAGNITUDE_PROBS, m.NEGATIVE_MAGNITUDE_PROBS),
        (
            "symmetric",
            m.SYMMETRIC_SIGN_PROBS,
            m.SYMMETRIC_MAGNITUDE_PROBS,
            m.SYMMETRIC_MAGNITUDE_PROBS,
        ),
    ]:
        cat = m.mutation_catalogue(kernel, 4)
        for sign, probability, magnitudes in [(-1, signs[0], negative), (1, signs[2], positive)]:
            for magnitude, mp in zip(m.MAGNITUDES, magnitudes):
                mask = np.all(cat["deltas"] == sign * magnitude, axis=1)
                assert mask.sum() == 1
                np.testing.assert_allclose(
                    cat["probs"][mask][0], probability**4 * mp / (1 - signs[1] ** 4), rtol=1e-14
                )


def test_fixation_extremes_and_preserved_boundaries():
    values = m.vectorized_fixation_probability(np.array([-1.0, -1e-14, 0.0, 1e-14, 1.0]), 1000000)
    assert np.isfinite(values).all()
    assert values[0] == 0 and values[2] == 1e-6
    assert values[1] == 0 and values[3] == 0  # Documented inherited exact-boundary masks.
    np.testing.assert_allclose(values[4], -np.expm1(-2.0), rtol=1e-14)
