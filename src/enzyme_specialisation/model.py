"""Public API for the enzyme-specialisation model."""

import numpy as np

from .ancestry import (
    balanced_column_degrees,
    build_initial_architecture,
    make_balanced_simple_network,
)
from .catalogue import mutation_catalogue
from .fitness import (
    exact_s_for_events,
    exact_s_single_function,
    lambda_vector,
    reaction_performance,
)
from .fixation import vectorized_fixation_probability
from .metrics import (
    effective_breadth,
    effective_reaction_redundancy,
    normalized_activity_profiles,
    record_state,
    strongest_activity_share,
)
from .parameters import (
    A0_MEAN,
    ACTIVE_THRESHOLD,
    EFFECTIVE_POP_SIZE,
    INITIAL_TOTAL_ACTIVITY,
    K_REF,
    LOOP_SIZES,
    LOOP_STRUCTURAL_BUDGET,
    LOOPS_PER_PROTEIN,
    LOSS_SIGN_PROBS,
    M_FUNCTIONS,
    MAGNITUDES,
    MAX_ACTIVITY,
    MAX_EFFECTIVE_GENERATIONS,
    MAX_LOOP_SUPPORT,
    NEGATIVE_MAGNITUDE_PROBS,
    POSITIVE_MAGNITUDE_PROBS,
    PROTEIN_FUNCTIONS,
    SUPPORT_EXPANSION_MAGNITUDE_PROBS,
    SUPPORT_EXPANSION_MAGNITUDES,
    SUPPORT_EXPANSION_PROB,
    SYMMETRIC_MAGNITUDE_PROBS,
    SYMMETRIC_SIGN_PROBS,
    U_LOOP,
    W_FUNCTION,
    RunConfig,
    eta_from_chi,
    n_from_redundancy,
    redundancy_from_n,
)
from .structure import loop_structural_cost, vectorized_loop_structural_cost


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
        n = n_from_redundancy(r0)
        for replicate in range(replicates):
            ancestry_seed = base_seed + 1000000 * n + replicate
            for kernel in kernels:
                for chi in chis:
                    evolution_seed = base_seed + 10000000 * condition_index + 10000 * replicate + n
                    save_history = (
                        np.isclose(r0, 12.5)
                        and kernel == "loss_biased"
                        and np.isclose(chi, 0.325)
                        and (replicate == 0)
                    )
                    jobs.append(
                        RunConfig(
                            n_proteins=n,
                            mutation_kernel=kernel,
                            chi=float(chi),
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


from .simulation import run_one_condition
