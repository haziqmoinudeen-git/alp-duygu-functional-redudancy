"""Deterministic event-by-event comparison to the preserved source."""

import hashlib
import importlib.util
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

from enzyme_specialisation import model


def load_original():
    spec = importlib.util.spec_from_file_location(
        "reference_model", Path("data/provenance/original_sweep.py")
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_event_catalogue_hazard_equivalence():
    old = load_original()
    cfg = model.RunConfig(
        40, "symmetric", 1.0, 0, 40000007, 47, 8.0, 100, 1000000, 1e-8, 2, 3e10, 1, True
    )
    events = []

    def capture(frame, event, arg):
        if event == "return" and frame.f_code.co_name == "evaluate_loop_events":
            arrays = tuple(
                (key, hashlib.sha256(value.tobytes()).hexdigest())
                for key, value in arg.items()
                if isinstance(value, np.ndarray)
            )
            events.append((arrays, arg["total_hazard"]))

    histories = []
    all_events = []
    for module in [old, model]:
        events = []
        sys.setprofile(capture)
        try:
            result = module.run_one_condition(module.RunConfig(**asdict(cfg)))
        finally:
            sys.setprofile(None)
        all_events.append(events)
        histories.append(result)
    assert len(all_events[0]) > 160
    assert all_events[0] == all_events[1]  # Candidate order, activities, selection and hazards.
    np.testing.assert_array_equal(histories[0]["A_final"], histories[1]["A_final"])


def test_latent_support_survives_zero_activity():
    c = model.mutation_catalogue("loss_biased", 4)
    current = np.zeros(4)
    candidates = np.maximum(current[None, :] + c["deltas"], 0.0)
    # Zero coordinates remain targets for ordinary positive candidates.
    assert np.all(np.any(candidates > 0, axis=0))
    # Truncated negative changes are zero and cannot create negative activity.
    assert np.all(candidates >= 0)
    negative = np.all(c["deltas"] <= 0, axis=1)
    assert np.all(candidates[negative] == 0)
