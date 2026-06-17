import numpy as np
import pytest
from numpy.testing import assert_allclose

from src.components import (
    all_pass_ring_response,
    coupling_power_from_theta,
    delay_line_response,
    directional_coupler_matrix,
    normalize_phase,
)


def test_directional_coupler_matrix_matches_matlab_form():
    matrix = directional_coupler_matrix(0.5)
    expected = np.array(
        [
            [np.sqrt(0.5), -1j * np.sqrt(0.5)],
            [-1j * np.sqrt(0.5), np.sqrt(0.5)],
        ],
        dtype=complex,
    )
    assert_allclose(matrix, expected)


def test_directional_coupler_rejects_non_physical_coupling():
    with pytest.raises(ValueError):
        directional_coupler_matrix(1.2)


def test_all_pass_ring_response_matches_expected_formula():
    w = np.array([0.0, np.pi / 2])
    k = 0.393
    a = 0.979888
    phase = np.pi
    expected = (
        np.sqrt(1 - k) - a * np.exp(-1j * (w + phase))
    ) / (
        1 - np.sqrt(1 - k) * a * np.exp(-1j * (w + phase))
    )
    assert_allclose(all_pass_ring_response(w, k, a, phase), expected)


def test_all_pass_ring_rejects_active_round_trip_amplitude():
    with pytest.raises(ValueError):
        all_pass_ring_response(np.array([0.0]), 0.5, 1.01, 0.0)


def test_delay_line_response_matches_matlab_form():
    w = np.array([0.0, np.pi])
    assert_allclose(delay_line_response(w, 0.9, 1.0, 0.0), 0.9 * np.exp(-1j * w))


def test_phase_normalization_and_theta_conversion():
    assert normalize_phase(3 * np.pi) == pytest.approx(np.pi)
    assert coupling_power_from_theta(np.pi / 2) == pytest.approx(0.5)
