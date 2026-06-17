"""RAMZI transfer-function simulation."""

from __future__ import annotations

import numpy as np

from .components import (
    all_pass_ring_response,
    delay_line_response,
    directional_coupler_matrix,
)


def create_frequency_grid(
    f_center: float,
    fsr: float,
    fsr_span: float,
    w_min: float,
    w_max: float,
    num_points: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Create physical-frequency and round-trip phase arrays.

    The normalized angular coordinate is defined as
    ``w = 2*pi*(frequency - f_center)/fsr``. The ``w_min`` and ``w_max``
    arguments are kept for compatibility with older scripts.
    """
    frequency = np.linspace(
        f_center - fsr_span * fsr,
        f_center + fsr_span * fsr,
        num_points,
    )
    w = 2.0 * np.pi * (frequency - f_center) / fsr
    return w, frequency


def split_ring_params(
    params: np.ndarray | list[float],
    n_upper: int,
    n_lower: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Split a flat coupling vector into upper-arm and lower-arm ring values."""
    params = np.asarray(params, dtype=float)
    expected = n_upper + n_lower
    if params.size != expected:
        raise ValueError(f"Expected {expected} ring coupling values, got {params.size}.")
    return params[:n_upper], params[n_upper:]


def ramzi_transfer(
    params: np.ndarray | list[float],
    w: np.ndarray,
    *,
    n_upper: int,
    n_lower: int,
    input_coupling: float,
    output_coupling: float,
    round_trip_amplitude: float,
    ring_phases: tuple[float, ...],
    delay: float,
    delay_phase: float,
) -> dict[str, np.ndarray]:
    """Simulate the four complex RAMZI transfer channels."""
    upper_k, lower_k = split_ring_params(params, n_upper, n_lower)
    expected_phases = n_upper + n_lower
    if len(ring_phases) != expected_phases:
        raise ValueError(f"Expected {expected_phases} ring phases, got {len(ring_phases)}.")

    upper = np.ones_like(w, dtype=complex)
    for k, phase in zip(upper_k, ring_phases[:n_upper]):
        upper *= all_pass_ring_response(w, k, round_trip_amplitude, phase)

    lower = np.ones_like(w, dtype=complex)
    for k, phase in zip(lower_k, ring_phases[n_upper:]):
        lower *= all_pass_ring_response(w, k, round_trip_amplitude, phase)
    lower *= delay_line_response(w, round_trip_amplitude, delay, delay_phase)

    h1 = directional_coupler_matrix(input_coupling)
    h3 = directional_coupler_matrix(output_coupling)
    h2 = np.zeros((w.size, 2, 2), dtype=complex)
    h2[:, 0, 0] = upper
    h2[:, 1, 1] = lower
    h = np.einsum("ab,nbc,cd->nad", h1, h2, h3)

    return {
        "H11": h[:, 0, 0],
        "H12": h[:, 0, 1],
        "H21": h[:, 1, 0],
        "H22": h[:, 1, 1],
    }


def simulate_power_db(
    params: np.ndarray | list[float],
    w: np.ndarray,
    *,
    output_channel: str = "H11",
    floor: float = 1e-12,
    **ramzi_kwargs,
) -> np.ndarray:
    """Return the selected RAMZI channel power response in dB."""
    transfer = ramzi_transfer(params, w, **ramzi_kwargs)
    if output_channel not in transfer:
        raise ValueError(f"Unknown output channel {output_channel!r}.")
    return 20.0 * np.log10(np.maximum(np.abs(transfer[output_channel]), floor))
