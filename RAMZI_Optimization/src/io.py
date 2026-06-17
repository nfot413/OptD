"""Data loading and result output helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def save_spectrum_csv(
    path: str | Path,
    *,
    frequency_hz: np.ndarray,
    target_db: np.ndarray,
    simulated_db: np.ndarray,
) -> Path:
    """Save target and simulated spectra to CSV."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = np.column_stack([frequency_hz, target_db, simulated_db])
    np.savetxt(
        path,
        data,
        delimiter=",",
        header="frequency_hz,target_db,simulated_db",
        comments="",
    )
    return path
