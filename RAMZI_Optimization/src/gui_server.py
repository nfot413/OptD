"""Local web GUI server for fixed RAMZI simulation and optimization."""

from __future__ import annotations

import argparse
import json
import threading
import webbrowser
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import numpy as np

from .architecture import (
    default_parameters,
    fixed_frequency_grid,
    fixed_ramzi_payload,
    optimize_fixed_ramzi,
    parameter_bounds,
    simulate_fixed_ramzi_db,
    square_target_db,
)
from .config import DEFAULT_CONFIG
from .objective import mean_squared_error

PROJECT_ROOT = Path(__file__).resolve().parents[1]
GUI_ROOT = PROJECT_ROOT / "gui"


def _frequency_grid(payload: dict[str, Any]) -> np.ndarray:
    frequency_config = payload.get("frequency", {})
    defaults = DEFAULT_CONFIG.frequency
    center_hz = float(frequency_config.get("center_hz", defaults.center_hz))
    span_hz = float(
        frequency_config.get(
            "span_hz",
            2.0 * defaults.fsr_span * defaults.fsr_hz,
        )
    )
    num_points = int(frequency_config.get("num_points", defaults.num_points))
    return fixed_frequency_grid(center_hz=center_hz, span_hz=span_hz, num_points=num_points)


def _parameters(payload: dict[str, Any]) -> dict[str, float]:
    parameters = default_parameters()
    parameters.update({name: float(value) for name, value in payload.get("parameters", {}).items()})
    return parameters


def _target(frequency_hz: np.ndarray, payload: dict[str, Any]) -> np.ndarray:
    target_config = payload.get("target", {})
    width_hz = float(target_config.get("passband_width_hz", 8e9))
    stopband_db = float(target_config.get("stopband_db", -35.0))
    passband_db = float(target_config.get("passband_db", 0.0))
    return square_target_db(
        frequency_hz,
        passband_width_hz=width_hz,
        stopband_db=stopband_db,
        passband_db=passband_db,
    )


def _simulate_payload(payload: dict[str, Any]) -> dict[str, Any]:
    frequency = _frequency_grid(payload)
    parameters = _parameters(payload)
    output_channel = payload.get("output_channel", "cross")
    target = _target(frequency, payload)
    simulated = simulate_fixed_ramzi_db(frequency, parameters, output_channel=output_channel)
    return {
        "frequency_hz": frequency.tolist(),
        "wavelength_nm": (3.0e8 / frequency / 1e-9).tolist(),
        "target_db": target.tolist(),
        "simulated_db": simulated.tolist(),
        "mse": mean_squared_error(simulated, target),
        "parameters": parameters,
    }


def _optimize_payload(payload: dict[str, Any]) -> dict[str, Any]:
    frequency = _frequency_grid(payload)
    target = _target(frequency, payload)
    output_channel = payload.get("output_channel", "cross")
    optimizer_payload = payload.get("optimizer", {})
    defaults = DEFAULT_CONFIG.optimizer
    result = optimize_fixed_ramzi(
        frequency,
        target,
        initial_parameters=_parameters(payload),
        max_iterations=int(optimizer_payload.get("max_iterations", 80)),
        population_size=int(optimizer_payload.get("population_size", 10)),
        random_seed=optimizer_payload.get("random_seed", defaults.random_seed),
        polish=bool(optimizer_payload.get("polish_result", True)),
        output_channel=output_channel,
    )
    parameters = result.optimized_parameters
    simulated = simulate_fixed_ramzi_db(frequency, parameters, output_channel=output_channel)
    return {
        "frequency_hz": frequency.tolist(),
        "wavelength_nm": (3.0e8 / frequency / 1e-9).tolist(),
        "target_db": target.tolist(),
        "simulated_db": simulated.tolist(),
        "mse": mean_squared_error(simulated, target),
        "best_loss": float(result.fun),
        "parameters": parameters,
        "variables": list(parameter_bounds()),
    }


class GuiRequestHandler(SimpleHTTPRequestHandler):
    """Serve GUI assets and JSON API endpoints."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, directory=str(GUI_ROOT), **kwargs)

    def do_GET(self):
        path = urlparse(self.path).path
        if path == "/api/default":
            self._send_json(
                {
                    **fixed_ramzi_payload(),
                    "frequency": {
                        "center_hz": DEFAULT_CONFIG.frequency.center_hz,
                        "span_hz": 2.0
                        * DEFAULT_CONFIG.frequency.fsr_span
                        * DEFAULT_CONFIG.frequency.fsr_hz,
                        "num_points": DEFAULT_CONFIG.frequency.num_points,
                    },
                    "target": {
                        "passband_width_hz": 8e9,
                        "passband_db": 0.0,
                        "stopband_db": -35.0,
                    },
                    "optimizer": {
                        "max_iterations": 80,
                        "population_size": 10,
                        "random_seed": DEFAULT_CONFIG.optimizer.random_seed,
                        "polish_result": True,
                    },
                }
            )
            return
        if path == "/":
            self.path = "/index.html"
        super().do_GET()

    def do_POST(self):
        path = urlparse(self.path).path
        try:
            payload = self._read_json()
            if path == "/api/simulate":
                self._send_json(_simulate_payload(payload))
                return
            if path == "/api/optimize":
                self._send_json(_optimize_payload(payload))
                return
            self.send_error(404, "Unknown API endpoint")
        except Exception as exc:
            self._send_json({"error": str(exc)}, status=400)

    def _read_json(self) -> dict[str, Any]:
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length).decode("utf-8")
        return json.loads(raw or "{}")

    def _send_json(self, payload: dict[str, Any], status: int = 200) -> None:
        encoded = json.dumps(payload).encode("utf-8")
        self.send_response(status)
        self.send_header("Content-Type", "application/json; charset=utf-8")
        self.send_header("Content-Length", str(len(encoded)))
        self.end_headers()
        self.wfile.write(encoded)


def run(host: str = "127.0.0.1", port: int = 8765, *, open_browser: bool = True) -> None:
    """Run the local GUI server."""
    server = ThreadingHTTPServer((host, port), GuiRequestHandler)
    url = f"http://{host}:{port}"
    print(f"RAMZI GUI: {url}")
    if open_browser:
        threading.Timer(0.5, lambda: webbrowser.open(url)).start()
    server.serve_forever()


def main() -> None:
    parser = argparse.ArgumentParser(description="Run RAMZI local GUI server.")
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--no-open", action="store_true", help="Start the server without opening a browser.")
    args = parser.parse_args()
    run(args.host, args.port, open_browser=not args.no_open)


if __name__ == "__main__":
    main()
