"""Command-line entry point for RAMZI optimization."""

from __future__ import annotations

from .experiment import run_experiment


def run() -> None:
    """Run the default RAMZI optimization experiment."""
    run_experiment()


if __name__ == "__main__":
    run()
