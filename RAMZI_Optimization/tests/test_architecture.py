import numpy as np
import pytest

from src.architecture import (
    DEFAULT_SPEC,
    default_parameters,
    fixed_frequency_grid,
    parameter_bounds,
    parameter_vector,
    parameters_from_vector,
    ring_fsr_summary,
    simulate_fixed_ramzi,
    simulate_fixed_ramzi_db,
    square_target_db,
)


def test_ring_fsr_comes_from_group_index_and_length():
    rings = ring_fsr_summary()

    assert rings[0]["fsr_hz"] == pytest.approx(3e8 / (350e-6 * 4.3))
    assert rings[1]["fsr_hz"] == pytest.approx(3e8 / (3000e-6 * 4.3))


def test_fixed_ramzi_simulates_matlab_structure_outputs():
    frequency = fixed_frequency_grid(num_points=51)
    transfer = simulate_fixed_ramzi(frequency)
    spectrum = simulate_fixed_ramzi_db(frequency)

    assert set(transfer) == {"bar", "cross", "ring1", "ring2", "ring3", "ring4"}
    assert spectrum.shape == frequency.shape
    assert np.isfinite(spectrum).all()


def test_square_target_is_centered_and_not_flat():
    frequency = fixed_frequency_grid(num_points=101)
    target = square_target_db(frequency, passband_width_hz=8e9)

    assert target.shape == frequency.shape
    assert target.min() < target.max()
    assert target[len(target) // 2] == pytest.approx(0.0)


def test_parameter_vector_round_trip():
    params = default_parameters()
    values = parameter_vector(params)

    assert values.shape == (7,)
    assert parameters_from_vector(values) == params
    assert set(parameter_bounds()) == set(params)
    assert len(DEFAULT_SPEC.rings) == 4
