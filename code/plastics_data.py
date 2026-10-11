"""Loader for the primary dataset, FTIR-Plastics (Zenodo record version 10736650).

The CSV stores metadata in the first six columns and then interleaved
``Data(x)`` / ``Data(Y)`` pairs, one pair per spectral point. Parsing the
98 MB (C8) and 198 MB (C4) files takes about a minute each, so the parsed
arrays are cached next to the CSV as ``.npz``; delete the cache to force a
re-parse. Nothing here filters, subsamples or transforms the data.
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, "data_raw", "ftir_plastics")
RESULTS = os.path.join(ROOT, "results")

#: Protocol Section 5: three windows fixed a priori from polymer chemistry.
WINDOWS = {"W1": (2800.0, 3000.0), "W2": (1650.0, 1800.0), "W3": (1350.0, 1500.0)}

#: Protocol Section 6: fixed master seed for every random generator.
SEED = 20261010

POLYMERS = ["HDPE", "LDPE", "PET", "PP", "PS", "PVC"]


def _parse(csv_path: str) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    raw = pd.read_csv(csv_path, header=None, skiprows=1, encoding="utf-8-sig")
    meta = raw.iloc[:, :6]
    meta.columns = ["IDE", "Polymer", "Technic", "Sample", "BR", "RST"]
    block = raw.iloc[:, 6:].to_numpy(dtype=float)
    X = block[:, 0::2]
    Y = block[:, 1::2]
    grid = X[0]
    if not np.allclose(X, grid[None, :], atol=1e-6):
        raise ValueError("wavenumber grid is not identical across specimens")
    return grid, Y, meta.reset_index(drop=True)


def load(resolution: str = "c8", use_cache: bool = True):
    """Return ``(grid, Y, meta)`` for ``c8`` (8 cm-1) or ``c4`` (4 cm-1)."""
    csv_path = os.path.join(DATA_DIR, f"FTIR_PLASTIC_{resolution}.csv")
    cache = os.path.join(DATA_DIR, f"_cache_{resolution}.npz")
    if use_cache and os.path.exists(cache):
        z = np.load(cache, allow_pickle=True)
        meta = pd.DataFrame({k: z[f"meta_{k}"] for k in
                             ["IDE", "Polymer", "Technic", "Sample", "BR", "RST"]})
        return z["grid"], z["Y"], meta
    grid, Y, meta = _parse(csv_path)
    np.savez_compressed(cache, grid=grid, Y=Y,
                        **{f"meta_{c}": meta[c].to_numpy().astype(str) for c in meta.columns})
    return grid, Y, meta


def results_path(name: str) -> str:
    os.makedirs(RESULTS, exist_ok=True)
    return os.path.join(RESULTS, name)


def write_csv(df: pd.DataFrame, name: str) -> str:
    """Every result table is written as UTF-8 with BOM, as the deposit requires."""
    path = results_path(name)
    df.to_csv(path, index=False, encoding="utf-8-sig")
    return path
