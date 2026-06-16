"""RAMZI simulation and optimization configuration."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class DataConfig:
    target_spectrum_file: str = "data/fsr1.xlsx"
    target_row: int = 0
    target_start_column: int = 1


@dataclass(frozen=True)
class RamziConfig:
    n_upper_rings: int = 2
    n_lower_rings: int = 1
    input_coupling: float = 0.5
    output_coupling: float = 0.5
    round_trip_amplitude: float = 0.979888
    ring_phases: tuple[float, ...] = (np.pi, np.pi, np.pi)
    delay: float = 1.0
    delay_phase: float = 0.0
    initial_couplings: tuple[float, ...] = (0.8673, 0.09, 0.393)
    coupling_bounds: tuple[float, float] = (0.0, 1.0)


@dataclass(frozen=True)
class FrequencyGridConfig:
    center_hz: float = 193.1e12
    fsr_hz: float = 100e9
    fsr_span: float = 10
    w_min: float = -20 * np.pi
    w_max: float = 20 * np.pi
    num_points: int = 1401


@dataclass(frozen=True)
class OptimizerConfig:
    max_iterations: int = 200
    population_size: int = 15
    random_seed: int | None = 123
    polish_result: bool = True


@dataclass(frozen=True)
class OutputConfig:
    results_dir: str = "results"
    plot_ylim: tuple[float, float] = (-60, 5)


@dataclass(frozen=True)
class ExperimentConfig:
    data: DataConfig = DataConfig()
    ramzi: RamziConfig = RamziConfig()
    frequency: FrequencyGridConfig = FrequencyGridConfig()
    optimizer: OptimizerConfig = OptimizerConfig()
    output: OutputConfig = OutputConfig()


DEFAULT_CONFIG = ExperimentConfig()

# Backward-compatible constants for small scripts and tests.
TARGET_SPECTRUM_FILE = DEFAULT_CONFIG.data.target_spectrum_file
TARGET_ROW = DEFAULT_CONFIG.data.target_row
TARGET_START_COLUMN = DEFAULT_CONFIG.data.target_start_column

N_UPPER_RINGS = DEFAULT_CONFIG.ramzi.n_upper_rings
N_LOWER_RINGS = DEFAULT_CONFIG.ramzi.n_lower_rings
INPUT_COUPLING = DEFAULT_CONFIG.ramzi.input_coupling
OUTPUT_COUPLING = DEFAULT_CONFIG.ramzi.output_coupling
ROUND_TRIP_AMPLITUDE = DEFAULT_CONFIG.ramzi.round_trip_amplitude
RING_PHASES = DEFAULT_CONFIG.ramzi.ring_phases
DELAY = DEFAULT_CONFIG.ramzi.delay
DELAY_PHASE = DEFAULT_CONFIG.ramzi.delay_phase
INITIAL_COUPLINGS = DEFAULT_CONFIG.ramzi.initial_couplings
COUPLING_BOUNDS = DEFAULT_CONFIG.ramzi.coupling_bounds

F_CENTER = DEFAULT_CONFIG.frequency.center_hz
FSR = DEFAULT_CONFIG.frequency.fsr_hz
FREQUENCY_FSR_SPAN = DEFAULT_CONFIG.frequency.fsr_span
W_MIN = DEFAULT_CONFIG.frequency.w_min
W_MAX = DEFAULT_CONFIG.frequency.w_max
NUM_POINTS = DEFAULT_CONFIG.frequency.num_points

MAX_ITERATIONS = DEFAULT_CONFIG.optimizer.max_iterations
POPULATION_SIZE = DEFAULT_CONFIG.optimizer.population_size
RANDOM_SEED = DEFAULT_CONFIG.optimizer.random_seed
POLISH_RESULT = DEFAULT_CONFIG.optimizer.polish_result

RESULTS_DIR = DEFAULT_CONFIG.output.results_dir
PLOT_YLIM = DEFAULT_CONFIG.output.plot_ylim
