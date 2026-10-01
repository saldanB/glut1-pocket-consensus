"""Volumetric agreement between detectors (Intersection over Union).

Each pocket is rasterised onto a common voxel grid: a voxel is occupied when it
lies within ``radius`` of any pocket-defining element of that tool in that frame
(alpha-sphere centres for geometry-based tools, energy grid points for
energy-based tools, atomic coordinates for the learning-based ones). The
dilation is computed with a Euclidean distance transform, which is equivalent to
dilating with a spherical element but linear in the number of voxels.

IoU is then computed per frame for each pair of tools and averaged over frames.
Note that ``radius`` is the tolerance of the criterion, not a free parameter of
the representation: the mean IoU depends on it monotonically, which is why the
sweep in run/04 reports it explicitly.

Point clouds are not shipped with the repository (31 MB); see data/README.md.
"""
from __future__ import annotations

import itertools
import os

import numpy as np
import pandas as pd
from scipy import ndimage


def load_point_cloud(path: str) -> tuple[np.ndarray, dict[int, list[tuple[int, int]]]]:
    """Load one tool's point cloud archive. Returns (points, frame -> slices)."""
    d = np.load(path, allow_pickle=True)
    points, counts, frames = d["points"], d["counts"], d["frames"]
    offsets = np.concatenate([[0], np.cumsum(counts)])
    by_frame: dict[int, list[tuple[int, int]]] = {}
    for i, f in enumerate(frames):
        by_frame.setdefault(int(f), []).append((offsets[i], offsets[i + 1]))
    return points, by_frame


def frame_points(cloud, frame: int) -> np.ndarray:
    points, by_frame = cloud
    seg = by_frame.get(frame)
    if not seg:
        return np.zeros((0, 3), np.float32)
    return np.concatenate([points[a:b] for a, b in seg])


def occupancy_mask(pts: np.ndarray, origin: np.ndarray, dims: np.ndarray,
                   voxel: float, radius: float) -> np.ndarray:
    grid = np.zeros(dims, bool)
    if len(pts):
        idx = np.clip(((pts - origin) / voxel).astype(int), 0, dims - 1)
        grid[idx[:, 0], idx[:, 1], idx[:, 2]] = True
    if not grid.any():
        return grid
    return ndimage.distance_transform_edt(~grid, sampling=voxel) <= radius


def pairwise_iou(clouds: dict[str, tuple], frames: np.ndarray, tools: list[str],
                 voxel: float, radius: float) -> pd.Series:
    """Frame-averaged IoU for every pair of tools."""
    every = [frame_points(clouds[t], int(f)) for t in tools for f in frames]
    every = np.concatenate([p for p in every if len(p)])
    lo, hi = every.min(0) - 8.0, every.max(0) + 8.0
    dims = np.ceil((hi - lo) / voxel).astype(int) + 1

    acc = {p: [] for p in itertools.combinations(tools, 2)}
    for f in frames:
        masks = {t: occupancy_mask(frame_points(clouds[t], int(f)), lo, dims, voxel, radius)
                 for t in tools}
        for a, b in acc:
            inter = np.count_nonzero(masks[a] & masks[b])
            union = np.count_nonzero(masks[a] | masks[b])
            acc[(a, b)].append(inter / union if union else 0.0)
    return pd.Series({f"{a}|{b}": float(np.mean(v)) for (a, b), v in acc.items()})


def clouds_available(directory: str, tools: list[str], system: str = "GLUT1") -> bool:
    return all(os.path.exists(os.path.join(directory, f"pc_{system}_{t}.npz")) for t in tools)


def load_all(directory: str, tools: list[str], system: str = "GLUT1") -> dict[str, tuple]:
    return {t: load_point_cloud(os.path.join(directory, f"pc_{system}_{t}.npz")) for t in tools}
