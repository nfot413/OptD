"""Optimization helpers for RAMZI simulation."""

from __future__ import annotations

import numpy as np
from scipy.optimize import OptimizeResult, differential_evolution

from .components import validate_power_coupling
from .objective import SpectrumObjective


def optimize_couplings(
    objective: SpectrumObjective,
    *,
    num_params: int,
    coupling_bounds: tuple[float, float],
    max_iterations: int,
    population_size: int,
    random_seed: int | None,
    polish: bool,
    workers: int = 1,
) -> OptimizeResult:
    """Optimize RAMZI ring coupling coefficients with differential evolution."""
    low, high = coupling_bounds
    validate_power_coupling(low, name="coupling lower bound")
    validate_power_coupling(high, name="coupling upper bound")
    if high <= low:
        raise ValueError(f"Invalid coupling bounds: {coupling_bounds}.")
    bounds = [coupling_bounds] * num_params
    return differential_evolution(
        objective,
        bounds=bounds,
        maxiter=max_iterations,
        popsize=population_size,
        seed=random_seed,
        polish=polish,
        workers=workers,
        updating="immediate" if workers == 1 else "deferred",
    )


def optimize_variables(
    objective,
    *,
    bounds: list[tuple[float, float]],
    max_iterations: int,
    population_size: int,
    random_seed: int | None,
    polish: bool,
    workers: int = 1,
) -> OptimizeResult:
    """Optimize an arbitrary vector of bounded variables."""
    for bound in bounds:
        if len(bound) != 2 or bound[1] <= bound[0]:
            raise ValueError(f"Invalid optimization bounds: {bound}.")
    return differential_evolution(
        objective,
        bounds=bounds,
        maxiter=max_iterations,
        popsize=population_size,
        seed=random_seed,
        polish=polish,
        workers=workers,
        updating="immediate" if workers == 1 else "deferred",
    )


def evaluate_initial_guess(
    objective: SpectrumObjective,
    initial_couplings: tuple[float, ...],
) -> tuple[np.ndarray, float]:
    """Evaluate the initial coupling values."""
    params = np.asarray(initial_couplings, dtype=float)
    return params, objective(params)
