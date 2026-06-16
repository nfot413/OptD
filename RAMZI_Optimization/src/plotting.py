"""Plotting helpers for RAMZI optimization results."""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

import numpy as np


def configure_matplotlib_cache() -> None:
    """Use writable cache directories before importing matplotlib.pyplot."""
    cache_dir = Path(tempfile.gettempdir()) / "ramzi_optimization_cache"
    cache_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir / "xdg"))


def save_spectrum_plot(
    path: str | Path,
    *,
    frequency_hz: np.ndarray,
    target_db: np.ndarray,
    simulated_db: np.ndarray,
    ylim: tuple[float, float],
) -> Path:
    """Save a target-vs-simulated spectrum plot."""
    configure_matplotlib_cache()
    import matplotlib.pyplot as plt

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(12, 6))
    plt.plot(frequency_hz, target_db, "r--", linewidth=1.5, label="target")
    plt.plot(frequency_hz, simulated_db, "b-", linewidth=1.2, label="optimized")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power (dB)")
    plt.ylim(*ylim)
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()
    return path
