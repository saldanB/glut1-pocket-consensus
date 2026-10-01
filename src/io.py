"""Paths, parameters and small shared helpers."""
from __future__ import annotations

import os

import numpy as np
import pandas as pd
import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RESULTS = os.path.join(ROOT, "results")


def parameters() -> dict:
    with open(os.path.join(DATA, "parameters.yaml")) as fh:
        return yaml.safe_load(fh)


def results_dir() -> str:
    os.makedirs(RESULTS, exist_ok=True)
    return RESULTS


def pockets(aligned: bool = True) -> pd.DataFrame:
    """The pocket table. ``aligned=True`` requires run/01 to have been executed."""
    name = "pockets_aligned.csv" if aligned else "pockets_master.csv"
    path = os.path.join(RESULTS if aligned else DATA, name)
    if not os.path.exists(path):
        raise SystemExit(f"missing {path}" + ("  — run: python run/01_align.py" if aligned else ""))
    return pd.read_csv(path)


def protein_ca() -> tuple[np.ndarray, np.ndarray]:
    z = np.load(os.path.join(DATA, "protein_ca.npz"))
    return z["ca"], z["frames"]


def published_centres() -> dict[str, np.ndarray]:
    d = pd.read_csv(os.path.join(DATA, "published_consensus_centres.csv"))
    return {ph: g[["x", "y", "z"]].to_numpy() for ph, g in d.groupby("phase")}


def save(frame: pd.DataFrame, name: str) -> str:
    path = os.path.join(results_dir(), name)
    frame.to_csv(path, index=False)
    print(f"  wrote {name}  ({len(frame)} rows)")
    return path


def banner(text: str) -> None:
    print("\n" + text)
    print("-" * len(text))
