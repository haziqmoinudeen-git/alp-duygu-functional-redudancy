"""Origin–fixation driver with vectorised expansion and atomic exact checkpoints."""

from __future__ import annotations

import os
from dataclasses import asdict
from pathlib import Path

import numpy as np

from .ancestry import build_initial_architecture
from .catalogue import mutation_catalogue
from .fitness import exact_s_for_events
from .fixation import vectorized_fixation_probability
from .io import scientific_fingerprint
from .metrics import record_state
from .parameters import (
    LOOP_STRUCTURAL_BUDGET,
    LOOPS_PER_PROTEIN,
    M_FUNCTIONS,
    MAX_ACTIVITY,
    MAX_LOOP_SUPPORT,
    SUPPORT_EXPANSION_MAGNITUDE_PROBS,
    SUPPORT_EXPANSION_MAGNITUDES,
    SUPPORT_EXPANSION_PROB,
    W_FUNCTION,
    RunConfig,
    eta_from_chi,
    redundancy_from_n,
)
from .structure import loop_structural_cost, vectorized_loop_structural_cost


def run_one_condition(cfg: RunConfig, checkpoint_path=None):
    eta = eta_from_chi(cfg.chi)
    L, A, latent_support = build_initial_architecture(cfg, eta)
    R = A.sum(axis=0)
    if np.any(R <= 0):
        raise RuntimeError("At least one organismal reaction begins with zero capacity.")
    K = R.copy()
    loop_ids = [(i, ell) for i in range(cfg.n_proteins) for ell in range(LOOPS_PER_PROTEIN)]
    n_loops_total = len(loop_ids)
    function_to_loops = [set() for _ in range(M_FUNCTIONS)]

    def rebuild_function_to_loops():
        for s in function_to_loops:
            s.clear()
        for lid, (i, ell) in enumerate(loop_ids):
            for j in np.flatnonzero(latent_support[i, ell]):
                function_to_loops[j].add(lid)

    rebuild_function_to_loops()
    origin_rate_per_loop = cfg.effective_pop_size * cfg.u_loop / n_loops_total
    ordinary_origin_rate_per_loop = origin_rate_per_loop * (1.0 - SUPPORT_EXPANSION_PROB)
    expansion_origin_rate_per_loop = origin_rate_per_loop * SUPPORT_EXPANSION_PROB
    rng = np.random.default_rng(cfg.evolution_seed)

    def evaluate_loop_events(loop_id: int):
        i, ell = loop_ids[loop_id]
        support = np.flatnonzero(latent_support[i, ell])
        k = len(support)
        current_loop = L[i, ell, support]
        catalogue = mutation_catalogue(cfg.mutation_kernel, k)
        raw_deltas = catalogue["deltas"]
        ordinary_probs = catalogue["probs"]
        candidate_loops = np.maximum(current_loop[None, :] + raw_deltas, 0.0)
        actual_deltas = candidate_loops - current_loop[None, :]
        state_changing = np.any(np.abs(actual_deltas) > 1e-15, axis=1)
        candidate_A = A[i, support][None, :] + actual_deltas
        within_bounds = np.all(candidate_A <= MAX_ACTIVITY + 1e-15, axis=1)
        structural_ok = (
            vectorized_loop_structural_cost(candidate_loops, eta) <= LOOP_STRUCTURAL_BUDGET + 1e-15
        )
        valid = state_changing & within_bounds & structural_ok
        s_exact = np.full(len(raw_deltas), -1.0, dtype=float)
        pfix = np.zeros(len(raw_deltas), dtype=float)
        if np.any(valid):
            s_exact[valid] = exact_s_for_events(R, K, support, actual_deltas[valid])
            pfix[valid] = vectorized_fixation_probability(s_exact[valid], cfg.effective_pop_size)
        ordinary_hazards = ordinary_origin_rate_per_loop * ordinary_probs * pfix
        expansion_functions = np.array([], dtype=int)
        expansion_magnitudes = np.array([], dtype=float)
        expansion_s = np.array([], dtype=float)
        expansion_hazards = np.array([], dtype=float)
        if k < MAX_LOOP_SUPPORT:
            candidates = np.flatnonzero(~latent_support[i, ell])
            if len(candidates) > 0:
                # Same candidate order and arithmetic as the scalar reference.
                magnitudes = SUPPORT_EXPANSION_MAGNITUDES
                structural = np.array(
                    [
                        loop_structural_cost(np.append(current_loop, magnitude), eta)
                        <= LOOP_STRUCTURAL_BUDGET + 1e-15
                        for magnitude in magnitudes
                    ]
                )
                functions = np.repeat(candidates, len(magnitudes))
                values = np.tile(magnitudes, len(candidates))
                probabilities = np.tile(SUPPORT_EXPANSION_MAGNITUDE_PROBS, len(candidates))
                feasible = np.tile(structural, len(candidates)) & (
                    A[i, functions] + values <= MAX_ACTIVITY + 1e-15
                )
                expansion_functions = functions[feasible]
                expansion_magnitudes = values[feasible]
                old_R = R[expansion_functions]
                new_R = old_R + expansion_magnitudes
                old_P = old_R / (K[expansion_functions] + old_R)
                new_P = new_R / (K[expansion_functions] + new_R)
                expansion_s = np.expm1(
                    W_FUNCTION[expansion_functions] * (np.log(new_P) - np.log(old_P))
                )
                expansion_hazards = (
                    (expansion_origin_rate_per_loop / len(candidates))
                    * probabilities[feasible]
                    * vectorized_fixation_probability(expansion_s, cfg.effective_pop_size)
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
            "total_hazard": float(ordinary_hazards.sum() + expansion_hazards.sum()),
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
    record_state(history, 0.0, 0, A, R, K, latent_support)
    counter_keys = (
        "fixed_support_expansions",
        "fixed_beneficial",
        "fixed_neutralish",
        "fixed_deleterious",
    )
    history.update({key: [0] for key in counter_keys})
    initial_breadth = history["mean_breadth"][0]
    initial_strongest_share = history["mean_strongest_share"][0]
    initial_weakest_P = history["weakest_P"][0]
    initial_redundancy = history["mean_redundancy"][0]
    loop_hazards = np.zeros(n_loops_total, dtype=float)
    for lid in range(n_loops_total):
        loop_hazards[lid] = evaluate_loop_events(lid)["total_hazard"]
    time = 0.0
    fixed = 0
    fixed_support_expansions = 0
    fixed_beneficial = 0
    fixed_neutralish = 0
    fixed_deleterious = 0
    stop_reason = "generation horizon"
    # Checkpoints are trusted local files, never load third-party pickle files.
    import pickle

    checkpoint_names = (
        "L",
        "A",
        "R",
        "K",
        "latent_support",
        "loop_hazards",
        "history",
        "time",
        "fixed",
        "fixed_support_expansions",
        "fixed_beneficial",
        "fixed_neutralish",
        "fixed_deleterious",
    )
    if checkpoint_path is not None and Path(checkpoint_path).exists():
        with open(checkpoint_path, "rb") as handle:
            saved = pickle.load(handle)
        if saved.get("scientific_fingerprint") != scientific_fingerprint():
            raise ValueError(
                "Checkpoint scientific source differs; preserve it and restart this trajectory"
            )
        if saved["cfg"] != asdict(cfg):
            raise ValueError("Checkpoint configuration mismatch")
        (
            L,
            A,
            R,
            K,
            latent_support,
            loop_hazards,
            history,
            time,
            fixed,
            fixed_support_expansions,
            fixed_beneficial,
            fixed_neutralish,
            fixed_deleterious,
        ) = [saved[name] for name in checkpoint_names]
        # Earlier checkpoint schemas have scalar counters but no counter time series.
        if not all(key in history for key in counter_keys):
            counter_keys = ()
        rng.bit_generator.state = saved["rng"]
        rebuild_function_to_loops()
    while time < cfg.max_generations:
        if cfg.max_fixed > 0 and fixed >= cfg.max_fixed:
            stop_reason = "fixed-substitution cap"
            break
        total_rate = float(loop_hazards.sum())
        if total_rate <= 0 or not np.isfinite(total_rate):
            stop_reason = "no accessible fixation"
            break
        dt = rng.exponential(1.0 / total_rate)
        if time + dt >= cfg.max_generations:
            time = float(cfg.max_generations)
            stop_reason = "generation horizon"
            break
        time += dt
        chosen_loop_id = int(rng.choice(n_loops_total, p=loop_hazards / total_rate))
        ev = evaluate_loop_events(chosen_loop_id)
        ordinary_total = float(ev["ordinary_hazards"].sum())
        expansion_total = float(ev["expansion_hazards"].sum())
        event_total = ordinary_total + expansion_total
        if event_total <= 0:
            loop_hazards[chosen_loop_id] = 0.0
            continue
        i, ell = loop_ids[chosen_loop_id]
        is_expansion = expansion_total > 0 and rng.random() < expansion_total / event_total
        affected_functions = set()
        if is_expansion:
            probs = ev["expansion_hazards"] / expansion_total
            eidx = int(rng.choice(len(probs), p=probs))
            j_new = int(ev["expansion_functions"][eidx])
            magnitude = float(ev["expansion_magnitudes"][eidx])
            s_exact = float(ev["expansion_s"][eidx])
            L[i, ell, j_new] += magnitude
            A[i, j_new] += magnitude
            R[j_new] += magnitude
            latent_support[i, ell, j_new] = True
            function_to_loops[j_new].add(chosen_loop_id)
            fixed_support_expansions += 1
            affected_functions.add(j_new)
        else:
            probs = ev["ordinary_hazards"] / ordinary_total
            oidx = int(rng.choice(len(probs), p=probs))
            support = ev["support"]
            new_loop_values = ev["candidate_loops"][oidx]
            delta = ev["actual_deltas"][oidx]
            s_exact = float(ev["ordinary_s"][oidx])
            L[i, ell, support] = new_loop_values
            A[i, support] += delta
            R[support] += delta
            A[i] = np.clip(A[i], 0.0, MAX_ACTIVITY)
            L[i, ell] = np.clip(L[i, ell], 0.0, None)
            changed = support[np.abs(delta) > 1e-15]
            affected_functions.update(changed.tolist())
        fixed += 1
        if s_exact > 1e-12:
            fixed_beneficial += 1
        elif s_exact < -1e-12:
            fixed_deleterious += 1
        else:
            fixed_neutralish += 1
        affected_loop_ids = {chosen_loop_id}
        for other_ell in range(LOOPS_PER_PROTEIN):
            affected_loop_ids.add(i * LOOPS_PER_PROTEIN + other_ell)
        for j in affected_functions:
            affected_loop_ids.update(function_to_loops[j])
        for lid in affected_loop_ids:
            loop_hazards[lid] = evaluate_loop_events(lid)["total_hazard"]
        if fixed % cfg.record_every == 0:
            record_state(history, time, fixed, A, R, K, latent_support)
            counts = (
                fixed_support_expansions,
                fixed_beneficial,
                fixed_neutralish,
                fixed_deleterious,
            )
            for key, count in zip(counter_keys, counts):
                history[key].append(count)
        if checkpoint_path is not None and fixed % 25 == 0:
            namespace = locals()
            saved = {name: namespace[name] for name in checkpoint_names}
            saved.update(
                cfg=asdict(cfg),
                rng=rng.bit_generator.state,
                scientific_fingerprint=scientific_fingerprint(),
            )
            temporary = str(checkpoint_path) + ".tmp"
            with open(temporary, "wb") as handle:
                pickle.dump(saved, handle, protocol=5)
            os.replace(temporary, checkpoint_path)
    if history["fixed"][-1] != fixed:
        record_state(history, time, fixed, A, R, K, latent_support)
        counts = (fixed_support_expansions, fixed_beneficial, fixed_neutralish, fixed_deleterious)
        for key, count in zip(counter_keys, counts):
            history[key].append(count)
    for key in history:
        history[key] = np.asarray(history[key])
    summary = {
        "n_proteins": cfg.n_proteins,
        "redundancy_r0": redundancy_from_n(cfg.n_proteins),
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
        "final_mean_breadth": float(history["mean_breadth"][-1]),
        "delta_mean_breadth": float(history["mean_breadth"][-1] - initial_breadth),
        "initial_strongest_share": initial_strongest_share,
        "final_strongest_share": float(history["mean_strongest_share"][-1]),
        "delta_strongest_share": float(
            history["mean_strongest_share"][-1] - initial_strongest_share
        ),
        "initial_weakest_P": initial_weakest_P,
        "final_weakest_P": float(history["weakest_P"][-1]),
        "delta_weakest_P": float(history["weakest_P"][-1] - initial_weakest_P),
        "initial_mean_redundancy": initial_redundancy,
        "final_mean_redundancy": float(history["mean_redundancy"][-1]),
        "delta_mean_redundancy": float(history["mean_redundancy"][-1] - initial_redundancy),
        "initial_mean_total_activity": float(history["mean_total_activity"][0]),
        "final_mean_total_activity": float(history["mean_total_activity"][-1]),
        "final_mean_nonzero_functions": float(history["mean_nonzero_functions"][-1]),
        "final_mean_active_functions": float(history["mean_active_functions"][-1]),
        "final_mean_latent_loop_support": float(history["mean_latent_loop_support"][-1]),
        "fixed_support_expansions": fixed_support_expansions,
        "fixed_beneficial": fixed_beneficial,
        "fixed_neutralish": fixed_neutralish,
        "fixed_deleterious": fixed_deleterious,
    }
    payload = {"summary": summary}
    if cfg.save_history:
        payload["history"] = history
        payload["A_final"] = A
    return payload
