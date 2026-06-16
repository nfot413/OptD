"""Objective functions for RAMZI simulation optimization."""

from __future__ import annotations

import numpy as np

from .io import load_target_spectrum
from .simulation import simulate_power_db


def mean_squared_error(simulated: np.ndarray, target: np.ndarray) -> float:
    """Return MSE after verifying matching shapes."""
    simulated = np.asarray(simulated, dtype=float)
    target = np.asarray(target, dtype=float)
    if simulated.shape != target.shape:
        raise ValueError(f"Shape mismatch: simulated {simulated.shape}, target {target.shape}.")
    return float(np.mean((simulated - target) ** 2))


class SpectrumObjective:
    """Callable MSE objective for RAMZI coupling optimization."""

    def __init__(self, w: np.ndarray, target_spectrum_db: np.ndarray, **simulation_kwargs):
        self.w = w
        self.target_spectrum_db = np.asarray(target_spectrum_db, dtype=float)
        self.simulation_kwargs = simulation_kwargs

    def __call__(self, params: np.ndarray | list[float]) -> float:
        simulated = simulate_power_db(params, self.w, **self.simulation_kwargs)
        return mean_squared_error(simulated, self.target_spectrum_db)
