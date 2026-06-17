"""Fixed RAMZI structure derived from the MATLAB reference model."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .objective import mean_squared_error

C0 = 3.0e8
TWO_PI = 2.0 * np.pi


@dataclass(frozen=True)
class RingSpec:
    """Physical parameters of one microring."""

    name: str
    length_m: float
    group_index: float = 4.3

    @property
    def fsr_hz(self) -> float:
        return C0 / (self.length_m * self.group_index)


@dataclass(frozen=True)
class FixedRamziSpec:
    """Fixed four-ring RAMZI physical structure."""

    input_coupling: float = 0.5
    output_coupling: float = 0.5
    input_mzi_phase: float = 0.5 * np.pi
    output_mzi_phase: float = 0.5 * np.pi
    loss_db_per_m: float = 15.0
    rings: tuple[RingSpec, ...] = (
        RingSpec("Ring 1", 350e-6),
        RingSpec("Ring 2", 3000e-6),
        RingSpec("Ring 3", 350e-6),
        RingSpec("Ring 4", 3000e-6),
    )


DEFAULT_SPEC = FixedRamziSpec()

DEFAULT_PARAMETERS = {
    "fai1": -0.0468 * np.pi,
    "fai2": -0.6842 * np.pi,
    "fai3": -0.0518 * np.pi,
    "fai4": -0.6198 * np.pi,
    "fait": 0.495 * np.pi,
    "theta13": -0.622 * np.pi,
    "theta24": -0.73 * np.pi,
}

PARAMETER_BOUNDS = {
    "fai1": (-np.pi, np.pi),
    "fai2": (-np.pi, np.pi),
    "fai3": (-np.pi, np.pi),
    "fai4": (-np.pi, np.pi),
    "fait": (-np.pi, np.pi),
    "theta13": (-np.pi, np.pi),
    "theta24": (-np.pi, np.pi),
}


def default_parameters() -> dict[str, float]:
    """Return mutable default RAMZI optimization parameters."""
    return {name: float(value) for name, value in DEFAULT_PARAMETERS.items()}


def parameter_bounds() -> dict[str, tuple[float, float]]:
    """Return optimization bounds for the fixed RAMZI parameters."""
    return {name: (float(low), float(high)) for name, (low, high) in PARAMETER_BOUNDS.items()}


def parameter_vector(parameters: dict[str, float]) -> np.ndarray:
    """Return parameter values in optimizer order."""
    return np.asarray([parameters[name] for name in PARAMETER_BOUNDS], dtype=float)


def parameters_from_vector(values: np.ndarray | list[float]) -> dict[str, float]:
    """Return a parameter dictionary from optimizer-order values."""
    values = np.asarray(values, dtype=float)
    if values.size != len(PARAMETER_BOUNDS):
        raise ValueError(f"Expected {len(PARAMETER_BOUNDS)} parameters, got {values.size}.")
    return {name: float(value) for name, value in zip(PARAMETER_BOUNDS, values)}


def fixed_frequency_grid(
    center_hz: float = 193.5405e12,
    span_hz: float = 35e9,
    num_points: int = 1401,
) -> np.ndarray:
    """Return the default MATLAB-reference frequency grid."""
    if span_hz <= 0.0:
        raise ValueError(f"span_hz must be positive, got {span_hz}.")
    if num_points < 2:
        raise ValueError(f"num_points must be at least 2, got {num_points}.")
    return np.linspace(center_hz - span_hz / 2.0, center_hz + span_hz / 2.0, num_points)


def square_target_db(
    frequency_hz: np.ndarray,
    *,
    passband_width_hz: float = 8e9,
    passband_db: float = 0.0,
    stopband_db: float = -35.0,
) -> np.ndarray:
    """Return a centered square-wave target spectrum in dB."""
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    if passband_width_hz <= 0.0:
        raise ValueError(f"passband_width_hz must be positive, got {passband_width_hz}.")
    center_hz = 0.5 * (frequency_hz[0] + frequency_hz[-1])
    half_width = 0.5 * passband_width_hz
    return np.where(np.abs(frequency_hz - center_hz) <= half_width, passband_db, stopband_db)


def ring_fsr_summary(spec: FixedRamziSpec = DEFAULT_SPEC) -> list[dict[str, float]]:
    """Return ring lengths and FSRs for display."""
    return [
        {
            "name": ring.name,
            "length_m": ring.length_m,
            "group_index": ring.group_index,
            "fsr_hz": ring.fsr_hz,
        }
        for ring in spec.rings
    ]


def _validate_power_coupling(value: float, name: str) -> float:
    value = float(value)
    if not 0.0 <= value <= 1.0:
        raise ValueError(f"{name} must be in [0, 1], got {value}.")
    return value


def _coupler_outputs(k_power: float, upper_in: np.ndarray, lower_in: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    k_power = _validate_power_coupling(k_power, "power coupling")
    through = np.sqrt(1.0 - k_power)
    cross = -1j * np.sqrt(k_power)
    upper_out = through * upper_in + cross * lower_in
    lower_out = cross * upper_in + through * lower_in
    return upper_out, lower_out


def _mzi_3db_outputs(
    upper_in: np.ndarray,
    lower_in: np.ndarray,
    *,
    coupling: float,
    phase: float,
) -> tuple[np.ndarray, np.ndarray]:
    e3, e4 = _coupler_outputs(coupling, upper_in, lower_in)
    e5 = e3
    e6 = np.exp(1j * phase) * e4
    return _coupler_outputs(coupling, e5, e6)


def _ring_loss_amplitude(length_m: float, loss_db_per_m: float) -> float:
    field_loss_db_per_m = loss_db_per_m / 2.0
    alpha = field_loss_db_per_m * np.log(10.0) / 10.0
    return float(np.exp(-alpha * length_m))


def _ring_response(
    frequency_hz: np.ndarray,
    *,
    phase: float,
    theta: float,
    ring: RingSpec,
    loss_db_per_m: float,
) -> np.ndarray:
    dynamic_phase = TWO_PI * frequency_hz * ring.length_m * ring.group_index / C0
    tau = _ring_loss_amplitude(ring.length_m, loss_db_per_m)
    phase_total = dynamic_phase + phase
    numerator = (
        (np.exp(1j * theta) - 1.0) / 2.0
        - tau * np.exp(1j * (phase_total + theta))
    )
    denominator = 1.0 - tau * (1.0 - np.exp(1j * theta)) * np.exp(1j * phase_total) / 2.0
    return numerator / denominator


def simulate_fixed_ramzi(
    frequency_hz: np.ndarray,
    parameters: dict[str, float] | None = None,
    *,
    spec: FixedRamziSpec = DEFAULT_SPEC,
) -> dict[str, np.ndarray]:
    """Simulate the fixed four-ring RAMZI structure."""
    frequency_hz = np.asarray(frequency_hz, dtype=float)
    params = default_parameters()
    if parameters:
        params.update({name: float(value) for name, value in parameters.items()})

    zeros = np.zeros_like(frequency_hz, dtype=complex)
    ones = np.ones_like(frequency_hz, dtype=complex)
    e7, e8 = _mzi_3db_outputs(
        zeros,
        ones,
        coupling=spec.input_coupling,
        phase=spec.input_mzi_phase,
    )

    c1 = _ring_response(
        frequency_hz,
        phase=params["fai1"],
        theta=params["theta13"],
        ring=spec.rings[0],
        loss_db_per_m=spec.loss_db_per_m,
    )
    c2 = _ring_response(
        frequency_hz,
        phase=params["fai2"],
        theta=params["theta24"],
        ring=spec.rings[1],
        loss_db_per_m=spec.loss_db_per_m,
    )
    c3 = _ring_response(
        frequency_hz,
        phase=params["fai3"],
        theta=params["theta13"],
        ring=spec.rings[2],
        loss_db_per_m=spec.loss_db_per_m,
    )
    c4 = _ring_response(
        frequency_hz,
        phase=params["fai4"],
        theta=params["theta24"],
        ring=spec.rings[3],
        loss_db_per_m=spec.loss_db_per_m,
    )

    upper = np.exp(1j * params["fait"]) * c1 * c2 * e7
    lower = np.exp(-1j * params["fait"]) * c3 * c4 * e8

    e11, e12 = _coupler_outputs(spec.output_coupling, upper, lower)
    e15, e16 = _mzi_3db_outputs(
        e11,
        e12,
        coupling=spec.output_coupling,
        phase=spec.output_mzi_phase,
    )
    return {
        "bar": e15,
        "cross": e16,
        "ring1": c1,
        "ring2": c2,
        "ring3": c3,
        "ring4": c4,
    }


def simulate_fixed_ramzi_db(
    frequency_hz: np.ndarray,
    parameters: dict[str, float] | None = None,
    *,
    output_channel: str = "cross",
    floor: float = 1e-12,
    spec: FixedRamziSpec = DEFAULT_SPEC,
) -> np.ndarray:
    """Return selected fixed-RAMZI output power in dB."""
    transfer = simulate_fixed_ramzi(frequency_hz, parameters, spec=spec)
    if output_channel not in transfer:
        raise ValueError(f"Unknown output channel {output_channel!r}.")
    power = np.maximum(np.abs(transfer[output_channel]) ** 2, floor)
    return 10.0 * np.log10(power)


def optimize_fixed_ramzi(
    frequency_hz: np.ndarray,
    target_db: np.ndarray,
    *,
    initial_parameters: dict[str, float] | None = None,
    max_iterations: int = 80,
    population_size: int = 10,
    random_seed: int | None = 123,
    polish: bool = True,
    output_channel: str = "cross",
):
    """Optimize the seven constrained phase parameters of the fixed RAMZI."""
    from .optimizer import optimize_variables

    initial = default_parameters()
    if initial_parameters:
        initial.update({name: float(value) for name, value in initial_parameters.items()})
    bounds = list(parameter_bounds().values())

    def objective(values):
        params = parameters_from_vector(values)
        simulated = simulate_fixed_ramzi_db(frequency_hz, params, output_channel=output_channel)
        return mean_squared_error(simulated, target_db)

    result = optimize_variables(
        objective,
        bounds=bounds,
        max_iterations=max_iterations,
        population_size=population_size,
        random_seed=random_seed,
        polish=polish,
    )
    result.initial_parameters = initial
    result.optimized_parameters = parameters_from_vector(result.x)
    return result


def fixed_ramzi_payload(parameters: dict[str, float] | None = None) -> dict[str, Any]:
    """Return serializable fixed RAMZI defaults."""
    return {
        "parameters": default_parameters() if parameters is None else parameters,
        "bounds": parameter_bounds(),
        "rings": ring_fsr_summary(),
        "output_channel": "cross",
    }
