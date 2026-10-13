"""Spatial consensus across pocket detectors.

Three steps, in this order:

1. per tool, HDBSCAN on the pocket centroids pooled over all frames of a phase,
   with ``min_cluster_size`` taken from ``data/parameters.yaml``;
2. the cluster centres of all tools are grouped by Ward hierarchical clustering
   cut at ``cutoff`` (6.5 A); the centre of a group is the simple mean of its
   member centres;
3. a region is a consensus region when at least ``min_tools`` distinct tools
   contribute to it.

Occupancy is the fraction of frames in which the region is detected by at
least one tool. Persistence follows the Methods: for each contributing tool,
the fraction of the phase covered by runs of detection lasting at least 20 % of
the phase, with gaps of at most ``delta`` frames (soft, delta = 5) or none
(strict); the region value is the mean over tools, with the inter-tool SD.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

try:
    import hdbscan as _hdbscan
except ImportError as exc:  # pragma: no cover
    raise SystemExit("hdbscan is required: conda env create -f environment.yml") from exc


def cluster_labels(xyz: np.ndarray, min_cluster_size: int) -> np.ndarray:
    """HDBSCAN labels for one tool's centroids. Label -1 is noise.

    The single place the clustering is performed, so that anything wanting the
    cluster membership rather than only the cluster centres -- the frame viewer in
    tools/viz -- sees exactly what the pipeline saw.
    """
    return _hdbscan.HDBSCAN(min_cluster_size=int(min_cluster_size)).fit_predict(xyz)


def cluster_one_tool(sub: pd.DataFrame, min_cluster_size: int) -> list[dict]:
    """HDBSCAN over the centroids of a single tool. Noise (label -1) is dropped."""
    xyz = sub[["x", "y", "z"]].to_numpy()
    if len(xyz) < min_cluster_size:
        return []
    labels = cluster_labels(xyz, min_cluster_size)
    frames = sub["frame"].to_numpy()
    out = []
    for lab in sorted(set(labels) - {-1}):
        m = labels == lab
        out.append({"centre": xyz[m].mean(axis=0),
                    "size": int(m.sum()),
                    "frames": np.unique(frames[m])})
    return out


def _runs_fraction(frames, all_frames, delta: int) -> float:
    """Fraction of the phase covered by runs of detection lasting >= 20 % of the phase.

    A run is a sequence of detections in which consecutive frames are separated
    by at most ``delta`` missing frames (delta = 5: soft; delta = 0: strict).
    Frames are indexed on ``all_frames`` so that subsampled sets are handled.
    """
    all_frames = np.sort(np.asarray(all_frames))
    n = all_frames.size
    rank = {v: i for i, v in enumerate(all_frames)}
    r = np.sort(np.unique([rank[v] for v in np.asarray(frames) if v in rank]))
    if r.size == 0 or n == 0:
        return 0.0
    runs, s, p = [], r[0], r[0]
    for x in r[1:]:
        if x - p - 1 <= delta:
            p = x
        else:
            runs.append((s, p)); s = p = x
    runs.append((s, p))
    return sum(e - s + 1 for s, e in runs if e - s + 1 >= 0.2 * n) / n


def persistence(frames, all_frames, delta: int = 5) -> tuple[float, float]:
    """Occupancy and soft persistence of a set of frames (kept for compatibility)."""
    all_frames = np.asarray(all_frames)
    f = np.unique(np.asarray(frames))
    occ = f.size / all_frames.size if all_frames.size else 0.0
    return occ, _runs_fraction(f, all_frames, delta)


def merge_into_regions(clusters_by_tool: dict[str, list[dict]],
                       cutoff: float = 6.5) -> list[dict]:
    """Ward hierarchical clustering of all cluster centres, cut at ``cutoff``.

    Each flat cluster of the dendrogram is a candidate region; its centre is the
    simple (unweighted) mean of the member cluster centres.
    """
    from scipy.cluster.hierarchy import fcluster, linkage
    nodes = [{"tool": tool, **c} for tool, cl in clusters_by_tool.items() for c in cl]
    if not nodes:
        return []
    pts = np.array([n["centre"] for n in nodes])
    if len(nodes) == 1:
        labels = np.array([1])
    else:
        labels = fcluster(linkage(pts, "ward"), cutoff, "distance")
    regions = []
    for lab in np.unique(labels):
        members = np.where(labels == lab)[0]
        tools = sorted({nodes[i]["tool"] for i in members})
        per_tool = {t: np.unique(np.concatenate([nodes[i]["frames"] for i in members
                                                 if nodes[i]["tool"] == t])) for t in tools}
        regions.append({
            "id": int(lab),
            "n_tools": len(tools),
            "tools": tools,
            "n_members": len(members),
            "centre": pts[members].mean(axis=0),
            "frames": np.unique(np.concatenate([nodes[i]["frames"] for i in members])),
            "frames_by_tool": per_tool,
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
    delta = params["consensus"]["delta_frames"]
    for r in regions:
        r["occupancy"] = np.unique(r["frames"]).size / all_frames.size
        soft = [_runs_fraction(f, all_frames, delta) for f in r["frames_by_tool"].values()]
        strict = [_runs_fraction(f, all_frames, 0) for f in r["frames_by_tool"].values()]
        r["soft_persistence"], r["soft_sd"] = float(np.mean(soft)), float(np.std(soft))
        r["strict_persistence"], r["strict_sd"] = float(np.mean(strict)), float(np.std(strict))
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
        "soft_sd": round(r.get("soft_sd", float("nan")), 3),
        "strict_persistence": round(r.get("strict_persistence", float("nan")), 3),
        "strict_sd": round(r.get("strict_sd", float("nan")), 3),
    } for i, r in enumerate(regions, 1)])
