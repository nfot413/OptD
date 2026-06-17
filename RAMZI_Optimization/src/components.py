"""Optical component transfer functions for RAMZI simulation."""

from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


def validate_power_coupling(k: float, *, name: str = "coupling") -> float:
    """Return a passive power coupling ratio K in [0, 1]."""
    k = float(k)
    if not 0.0 <= k <= 1.0:
        raise ValueError(f"{name} must be a power coupling ratio in [0, 1], got {k}.")
    return k


def validate_coupling(k: float) -> float:
    """Backward-compatible alias for power coupling validation."""
    return validate_power_coupling(k)


def validate_passive_amplitude(
    amplitude: float,
    *,
    name: str = "amplitude",
    allow_zero: bool = True,
) -> float:
    """Return a passive field-amplitude transmission value."""
    amplitude = float(amplitude)
    lower_ok = amplitude >= 0.0 if allow_zero else amplitude > 0.0
    if not lower_ok or amplitude > 1.0:
        lower = "[0, 1]" if allow_zero else "(0, 1]"
        raise ValueError(f"{name} must be a passive field amplitude in {lower}, got {amplitude}.")
    return amplitude


def normalize_phase(phase: float) -> float:
    """Return a radian phase normalized modulo 2*pi."""
    return float(phase) % TWO_PI


def coupling_power_from_theta(theta: float) -> float:
    """Convert a coupling angle theta in radians to power coupling K."""
    theta = float(theta)
    if not 0.0 <= theta <= np.pi:
        raise ValueError(f"Coupling angle theta must be in [0, pi], got {theta}.")
    return float(np.sin(theta / 2.0) ** 2)


def directional_coupler_matrix(k: float) -> np.ndarray:
    """Return a lossless 2x2 directional-coupler matrix from power coupling K."""
    k = validate_power_coupling(k, name="power coupling K")
    through = np.sqrt(1.0 - k)
    cross = -1j * np.sqrt(k)
    return np.array([[through, cross], [cross, through]], dtype=complex)


def all_pass_ring_response(
    w: np.ndarray,
    k: float,
    round_trip_amplitude: float,
    phase: float,
) -> np.ndarray:
    """Return the all-pass MRR response.

    ``k`` is the power coupling ratio K. ``round_trip_amplitude`` is the
    electric-field amplitude retained after one ring round trip.
    """
    k = validate_power_coupling(k, name="ring power coupling K")
    a = validate_passive_amplitude(
        round_trip_amplitude,
        name="round_trip_amplitude",
        allow_zero=False,
    )
    phase = normalize_phase(phase)
    phase_term = np.exp(-1j * (w + phase))
    self_coupling = np.sqrt(1.0 - k)
    numerator = self_coupling - a * phase_term
    denominator = 1.0 - self_coupling * a * phase_term
    return numerator / denominator


def delay_line_response(
    w: np.ndarray,
    amplitude: float,
    delay: float,
    phase: float,
) -> np.ndarray:
    """Return a lower-arm delay-line response."""
    amplitude = validate_passive_amplitude(amplitude, name="delay amplitude")
    delay = float(delay)
    if delay < 0.0:
        raise ValueError(f"delay must be non-negative, got {delay}.")
    phase = normalize_phase(phase)
    return amplitude * np.exp(-1j * (w * delay + phase))
