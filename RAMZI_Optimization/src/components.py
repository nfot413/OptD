"""Optical component transfer functions for RAMZI simulation."""

from __future__ import annotations

import numpy as np


def validate_coupling(k: float) -> float:
    """Return k as float after checking that it is a physical coupling value."""
    k = float(k)
    if not 0.0 <= k <= 1.0:
        raise ValueError(f"Coupling coefficient must be in [0, 1], got {k}.")
    return k


def directional_coupler_matrix(k: float) -> np.ndarray:
    """Return the 2x2 directional coupler matrix used in the MATLAB prototype."""
    k = validate_coupling(k)
    through = np.sqrt(1.0 - k)
    cross = -1j * np.sqrt(k)
    return np.array([[through, cross], [cross, through]], dtype=complex)


def all_pass_ring_response(
    w: np.ndarray,
    k: float,
    round_trip_amplitude: float,
    phase: float,
) -> np.ndarray:
    """Return the all-pass MRR response."""
    k = validate_coupling(k)
    a2 = float(round_trip_amplitude) ** 2
    phase_term = np.exp(-1j * (2.0 * w + phase))
    self_coupling = np.sqrt(1.0 - k)
    numerator = self_coupling - a2 * phase_term
    denominator = 1.0 - self_coupling * a2 * phase_term
    return numerator / denominator


def delay_line_response(
    w: np.ndarray,
    amplitude: float,
    delay: float,
    phase: float,
) -> np.ndarray:
    """Return a lower-arm delay-line response."""
    return float(amplitude) * np.exp(-1j * (w * delay + phase))
