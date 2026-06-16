import numpy as np
import pytest

from src import config
from src.objective import SpectrumObjective, mean_squared_error
from src.simulation import create_frequency_grid, simulate_power_db, split_ring_params


def _simulation_kwargs():
    return {
        "n_upper": config.N_UPPER_RINGS,
        "n_lower": config.N_LOWER_RINGS,
        "input_coupling": config.INPUT_COUPLING,
        "output_coupling": config.OUTPUT_COUPLING,
        "round_trip_amplitude": config.ROUND_TRIP_AMPLITUDE,
        "ring_phases": config.RING_PHASES,
        "delay": config.DELAY,
        "delay_phase": config.DELAY_PHASE,
    }


def test_create_frequency_grid_shapes():
    w, frequency = create_frequency_grid(193.1e12, 100e9, 10, -np.pi, np.pi, 101)
    assert w.shape == (101,)
    assert frequency.shape == (101,)
    assert frequency[0] == pytest.approx(192.1e12)
    assert frequency[-1] == pytest.approx(194.1e12)


def test_split_ring_params_validates_length():
    upper, lower = split_ring_params([0.1, 0.2, 0.3], 2, 1)
    assert upper.tolist() == [0.1, 0.2]
    assert lower.tolist() == [0.3]
    with pytest.raises(ValueError):
        split_ring_params([0.1, 0.2], 2, 1)


def test_simulate_power_db_returns_real_spectrum():
    w = np.linspace(-np.pi, np.pi, 21)
    spectrum = simulate_power_db(config.INITIAL_COUPLINGS, w, **_simulation_kwargs())
    assert spectrum.shape == w.shape
    assert np.isreal(spectrum).all()
    assert np.isfinite(spectrum).all()


def test_spectrum_objective_zero_for_identical_target():
    w = np.linspace(-np.pi, np.pi, 21)
    target = simulate_power_db(config.INITIAL_COUPLINGS, w, **_simulation_kwargs())
    objective = SpectrumObjective(w, target, **_simulation_kwargs())
    assert objective(config.INITIAL_COUPLINGS) == pytest.approx(0.0)


def test_mean_squared_error_rejects_mismatched_shapes():
    with pytest.raises(ValueError):
        mean_squared_error(np.zeros(2), np.zeros(3))
