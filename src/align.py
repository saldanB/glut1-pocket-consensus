"""Rigid-body superposition of the ensemble onto a reference frame.

The deposited per-frame coordinate files are *not* superposed: they are the
frames as supplied to the detectors. The superposition is applied at the
analysis stage, to the structure together with its pocket pseudo-atoms, before
centroids from different frames are pooled. It is a Kabsch fit on the protein
Calpha atoms with the first frame as reference, equivalent to MDAnalysis
``align.alignto`` with an unweighted CA selection.

Skipping this step changes the result: see run/02, which reports the consensus
with and without it.
"""
from __future__ import annotations

import numpy as np
import pandas as pd


def kabsch(mobile: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotation taking centred `mobile` onto centred `target`."""
    u, _, vt = np.linalg.svd(mobile.T @ target)
    d = np.sign(np.linalg.det(u @ vt))
    return u @ np.diag([1.0, 1.0, d]) @ vt


def superpose_pockets(pockets: pd.DataFrame, ca: np.ndarray, frames: np.ndarray,
                      reference_frame: int = 1) -> pd.DataFrame:
    """Apply, to every pocket centroid, the transform that superposes its frame.

    ``ca`` has shape (n_frames, n_residues, 3) and ``frames`` gives the frame
    number of each row of ``ca``.
    """
    index = {int(f): i for i, f in enumerate(frames)}
    if reference_frame not in index:
        raise ValueError(f"reference frame {reference_frame} not in the CA array")
    ref = ca[index[reference_frame]]
    ref_centre = ref.mean(axis=0)
    ref_centred = ref - ref_centre

    out = pockets.copy()
    xyz = out[["x", "y", "z"]].to_numpy(dtype=float)
    moved = np.empty_like(xyz)
    for frame, rows in out.groupby("frame").indices.items():
        f = int(frame)
        if f not in index:
            moved[rows] = xyz[rows]
            continue
        this = ca[index[f]]
        centre = this.mean(axis=0)
        rot = kabsch(this - centre, ref_centred)
        moved[rows] = (xyz[rows] - centre) @ rot + ref_centre
    out[["x", "y", "z"]] = moved
    return out


def displacement(before: pd.DataFrame, after: pd.DataFrame) -> np.ndarray:
    """Distance each centroid travels under the superposition."""
    key = ["tool", "frame", "pocket_id"]
    a = before.set_index(key)[["x", "y", "z"]]
    b = after.set_index(key)[["x", "y", "z"]]
    common = a.index.intersection(b.index)
    return np.linalg.norm(a.loc[common].to_numpy() - b.loc[common].to_numpy(), axis=1)
