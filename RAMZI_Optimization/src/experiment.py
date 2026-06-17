"""High-level workflow for fixed RAMZI square-wave optimization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from .architecture import (
    default_parameters,
    fixed_frequency_grid,
    optimize_fixed_ramzi,
    parameter_vector,
    simulate_fixed_ramzi_db,
    square_target_db,
)
from .config import DEFAULT_CONFIG, ExperimentConfig
from .io import save_spectrum_csv
from .objective import mean_squared_error
from .plotting import save_spectrum_plot


@dataclass(frozen=True)
class ExperimentResult:
    initial_params: np.ndarray
    initial_loss: float
    best_params: np.ndarray
    best_loss: float
    results_dir: Path
    spectrum_csv: Path
    spectrum_plot: Path


def project_root() -> Path:
    """Return the RAMZI_Optimization project root."""
    return Path(__file__).resolve().parents[1]


def resolve_project_path(path: str | Path) -> Path:
    """Resolve a path relative to the RAMZI_Optimization project root."""
    path = Path(path)
    if path.is_absolute():
        return path
    return project_root() / path


def format_values(values) -> list[float]:
    """Format numeric vectors for readable logs."""
    return [round(float(v), 6) for v in values]


def run_experiment(
    config: ExperimentConfig = DEFAULT_CONFIG,
    *,
    log: Callable[[str], None] = print,
) -> ExperimentResult:
    """Run fixed RAMZI optimization against a centered square-wave target."""
    frequency_config = config.frequency
    frequency = fixed_frequency_grid(
        center_hz=frequency_config.center_hz,
        span_hz=2.0 * frequency_config.fsr_span * frequency_config.fsr_hz,
        num_points=frequency_config.num_points,
    )
    target = square_target_db(frequency, passband_width_hz=8e9)

    initial_parameters = default_parameters()
    initial_spectrum = simulate_fixed_ramzi_db(frequency, initial_parameters)
    initial_loss = mean_squared_error(initial_spectrum, target)
    initial_params = parameter_vector(initial_parameters)

    log("固定 RAMZI 方波优化")
    log(f"目标光谱点数: {target.size}")
    log(f"初始参数: {format_values(initial_params)}")
    log(f"初始 MSE: {initial_loss:.6f}")

    optimizer_config = config.optimizer
    result = optimize_fixed_ramzi(
        frequency,
        target,
        initial_parameters=initial_parameters,
        max_iterations=optimizer_config.max_iterations,
        population_size=optimizer_config.population_size,
        random_seed=optimizer_config.random_seed,
        polish=optimizer_config.polish_result,
    )

    final_spectrum = simulate_fixed_ramzi_db(frequency, result.optimized_parameters)
    results_dir = resolve_project_path(config.output.results_dir)
    spectrum_csv = save_spectrum_csv(
        results_dir / "optimized_spectrum.csv",
        frequency_hz=frequency,
        target_db=target,
        simulated_db=final_spectrum,
    )
    spectrum_plot = save_spectrum_plot(
        results_dir / "optimized_spectrum.png",
        frequency_hz=frequency,
        target_db=target,
        simulated_db=final_spectrum,
        ylim=config.output.plot_ylim,
    )

    best_params = parameter_vector(result.optimized_parameters)
    log(f"最优参数: {format_values(best_params)}")
    log(f"最优 MSE: {result.fun:.6f}")
    log(f"结果目录: {results_dir}")

    return ExperimentResult(
        initial_params=initial_params,
        initial_loss=initial_loss,
        best_params=best_params,
        best_loss=float(result.fun),
        results_dir=results_dir,
        spectrum_csv=spectrum_csv,
        spectrum_plot=spectrum_plot,
    )
