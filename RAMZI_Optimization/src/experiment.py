"""High-level experiment workflow for RAMZI optimization."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np

from .config import DEFAULT_CONFIG, ExperimentConfig
from .io import load_target_spectrum, save_spectrum_csv
from .objective import SpectrumObjective
from .optimizer import evaluate_initial_guess, optimize_couplings
from .plotting import save_spectrum_plot
from .simulation import create_frequency_grid, simulate_power_db


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


def ramzi_simulation_kwargs(config: ExperimentConfig) -> dict:
    """Collect RAMZI simulation settings from an experiment config."""
    ramzi = config.ramzi
    return {
        "n_upper": ramzi.n_upper_rings,
        "n_lower": ramzi.n_lower_rings,
        "input_coupling": ramzi.input_coupling,
        "output_coupling": ramzi.output_coupling,
        "round_trip_amplitude": ramzi.round_trip_amplitude,
        "ring_phases": ramzi.ring_phases,
        "delay": ramzi.delay,
        "delay_phase": ramzi.delay_phase,
    }


def format_values(values) -> list[float]:
    """Format numeric vectors for readable logs."""
    return [round(float(v), 6) for v in values]


def run_experiment(
    config: ExperimentConfig = DEFAULT_CONFIG,
    *,
    log: Callable[[str], None] = print,
) -> ExperimentResult:
    """Run one RAMZI optimization experiment and save its outputs."""
    frequency_config = config.frequency
    w, frequency = create_frequency_grid(
        frequency_config.center_hz,
        frequency_config.fsr_hz,
        frequency_config.fsr_span,
        frequency_config.w_min,
        frequency_config.w_max,
        frequency_config.num_points,
    )

    data_config = config.data
    target = load_target_spectrum(
        resolve_project_path(data_config.target_spectrum_file),
        row=data_config.target_row,
        start_column=data_config.target_start_column,
        target_length=w.size,
    )
    objective = SpectrumObjective(w, target, **ramzi_simulation_kwargs(config))

    initial_params, initial_loss = evaluate_initial_guess(
        objective,
        config.ramzi.initial_couplings,
    )
    log("RAMZI 仿真优化")
    log(f"目标光谱点数: {target.size}")
    log(f"初始耦合系数: {format_values(initial_params)}")
    log(f"初始 MSE: {initial_loss:.6f}")

    optimizer_config = config.optimizer
    result = optimize_couplings(
        objective,
        num_params=config.ramzi.n_upper_rings + config.ramzi.n_lower_rings,
        coupling_bounds=config.ramzi.coupling_bounds,
        max_iterations=optimizer_config.max_iterations,
        population_size=optimizer_config.population_size,
        random_seed=optimizer_config.random_seed,
        polish=optimizer_config.polish_result,
    )

    final_spectrum = simulate_power_db(result.x, w, **ramzi_simulation_kwargs(config))
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

    log(f"最优耦合系数: {format_values(result.x)}")
    log(f"最优 MSE: {result.fun:.6f}")
    log(f"结果目录: {results_dir}")

    return ExperimentResult(
        initial_params=initial_params,
        initial_loss=initial_loss,
        best_params=np.asarray(result.x, dtype=float),
        best_loss=float(result.fun),
        results_dir=results_dir,
        spectrum_csv=spectrum_csv,
        spectrum_plot=spectrum_plot,
    )
