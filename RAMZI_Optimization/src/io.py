"""Data loading and result output helpers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def load_target_spectrum(
    path: str | Path,
    *,
    row: int = 0,
    start_column: int = 0,
    target_length: int | None = None,
) -> np.ndarray:
    """Load a target spectrum row from Excel and optionally resample it."""
    df = pd.read_excel(path, header=None)
    spectrum = df.iloc[row, start_column:].dropna().to_numpy(dtype=float)
    if target_length is not None and spectrum.size != target_length:
        source_x = np.linspace(0.0, 1.0, spectrum.size)
        target_x = np.linspace(0.0, 1.0, target_length)
        spectrum = np.interp(target_x, source_x, spectrum)
    return spectrum


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
    pd.DataFrame(
        {
            "frequency_hz": frequency_hz,
            "target_db": target_db,
            "simulated_db": simulated_db,
        }
    ).to_csv(path, index=False)
    return path
