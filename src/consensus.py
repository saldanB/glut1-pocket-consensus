"""Spatial consensus across pocket detectors.

Three steps, in this order:

1. per tool, HDBSCAN on the pocket centroids pooled over all frames of a phase,
   with ``min_cluster_size`` taken from ``data/parameters.yaml``;
2. cluster centres belonging to *different* tools and lying within ``cutoff`` of
   one another are merged into a consensus region (connected components of the
   proximity graph);
3. a region is a consensus region when at least ``min_tools`` distinct tools
   contribute to it.

Occupancy is the fraction of frames in which the region is detected. Soft
persistence is the fraction covered by the longest run of frames, tolerating
gaps of at most ``delta`` frames; both are normalised by the number of frames
actually present in the reference set, so that subsampled sets are handled
correctly.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

try:
    import hdbscan as _hdbscan
except ImportError as exc:  # pragma: no cover
    raise SystemExit("hdbscan is required: conda env create -f environment.yml") from exc


def cluster_one_tool(sub: pd.DataFrame, min_cluster_size: int) -> list[dict]:
    """HDBSCAN over the centroids of a single tool. Noise (label -1) is dropped."""
    xyz = sub[["x", "y", "z"]].to_numpy()
    if len(xyz) < min_cluster_size:
        return []
    labels = _hdbscan.HDBSCAN(min_cluster_size=int(min_cluster_size)).fit_predict(xyz)
    frames = sub["frame"].to_numpy()
    out = []
    for lab in sorted(set(labels) - {-1}):
        m = labels == lab
        out.append({"centre": xyz[m].mean(axis=0),
                    "size": int(m.sum()),
                    "frames": np.unique(frames[m])})
    return out


def persistence(frames, all_frames, delta: int = 5) -> tuple[float, float]:
    """Occupancy and soft persistence of one region."""
    all_frames = np.sort(np.asarray(all_frames))
    frames = np.sort(np.unique(np.asarray(frames)))
    if frames.size == 0 or all_frames.size == 0:
        return 0.0, 0.0
    occupancy = frames.size / all_frames.size
    rank = {v: i for i, v in enumerate(all_frames)}
    r = np.array([rank[v] for v in frames if v in rank])
    if r.size == 0:
        return occupancy, 0.0
    best = run = 1
    for a, b in zip(r[:-1], r[1:]):
        run = run + 1 if b - a <= delta + 1 else 1
        best = max(best, run)
    return occupancy, best / all_frames.size


def merge_into_regions(clusters_by_tool: dict[str, list[dict]],
                       cutoff: float = 6.5) -> list[dict]:
    """Connected components of the inter-tool proximity graph."""
    nodes = [{"tool": tool, **c} for tool, cl in clusters_by_tool.items() for c in cl]
    if not nodes:
        return []
    pts = np.array([n["centre"] for n in nodes])
    parent = list(range(len(nodes)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, j in cKDTree(pts).query_pairs(cutoff):
        if nodes[i]["tool"] != nodes[j]["tool"]:            # only across tools
            ri, rj = find(i), find(j)
            if ri != rj:
                parent[ri] = rj

    groups: dict[int, list[int]] = {}
    for i in range(len(nodes)):
        groups.setdefault(find(i), []).append(i)

    regions = []
    for members in groups.values():
        tools = sorted({nodes[i]["tool"] for i in members})
        weights = [nodes[i]["size"] for i in members]
        regions.append({
            "n_tools": len(tools),
            "tools": tools,
            "n_members": len(members),
            "centre": np.average([nodes[i]["centre"] for i in members], axis=0, weights=weights),
            "frames": np.unique(np.concatenate([nodes[i]["frames"] for i in members])),
        })
    regions.sort(key=lambda r: (-r["n_tools"], -r["n_members"]))
    return regions


def run_pipeline(df: pd.DataFrame, phase: str, params: dict, *,
                 min_cluster_size: dict | None = None,
                 tools: list[str] | None = None,
                 frames: np.ndarray | None = None,
                 top_k: int | None = None) -> tuple[dict, list[dict]]:
    """Full pipeline for one phase. Returns (clusters per tool, consensus regions)."""
    tools = tools or params["tools"]
    mcs = min_cluster_size or params["min_cluster_size"][phase]
    sub = df[(df["phase"] == phase) & (df["tool"].isin(tools))]
    if frames is not None:
        sub = sub[sub["frame"].isin(frames)]
    if top_k is not None:
        sub = sub[sub["rank"] <= top_k]
    all_frames = np.sort(sub["frame"].unique())

    clusters = {t: cluster_one_tool(sub[sub["tool"] == t], mcs[t]) for t in tools}
    regions = merge_into_regions(clusters, params["consensus"]["cutoff_A"])
    for r in regions:
        r["occupancy"], r["soft_persistence"] = persistence(
            r["frames"], all_frames, params["consensus"]["delta_frames"])
    return clusters, regions


def select(regions: list[dict], min_tools: int) -> list[dict]:
    return [r for r in regions if r["n_tools"] >= min_tools]


def recovery(reference: np.ndarray, regions: list[dict]) -> np.ndarray:
    """Distance from each reference point to the nearest region centre."""
    if not len(regions):
        return np.full(len(reference), np.inf)
    centres = np.array([r["centre"] for r in regions])
    return cKDTree(centres).query(np.asarray(reference))[0]


def regions_to_frame(regions: list[dict], phase: str, label: str = "") -> pd.DataFrame:
    return pd.DataFrame([{
        "scenario": label, "phase": phase, "region": i,
        "n_tools": r["n_tools"], "tools": "|".join(r["tools"]),
        "x": round(float(r["centre"][0]), 2),
        "y": round(float(r["centre"][1]), 2),
        "z": round(float(r["centre"][2]), 2),
        "occupancy": round(r["occupancy"], 3),
        "soft_persistence": round(r["soft_persistence"], 3),
    } for i, r in enumerate(regions, 1)])
