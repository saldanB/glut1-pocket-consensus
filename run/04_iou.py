#!/usr/bin/env python3
"""Volumetric agreement and its parameter sensitivity (Editor E1, R3 #8, R4 #1).

Needs the per-tool point clouds in data/pointclouds/, one archive per tool; see
data/README.md. The script exits with a clear message if they are absent.

Writes results/04_iou_sensitivity.csv, results/04_iou_pair_ranking.csv
"""
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.iou import clouds_available, load_all, pairwise_iou

p = io.parameters()
TOOLS = p["tools"]; C = p["iou"]
clouds_dir = os.path.join(io.DATA, "pointclouds")
if not clouds_available(clouds_dir, TOOLS):
    sys.exit(f"point clouds not found in {clouds_dir}\n"
             "all six pc_GLUT1_<tool>.npz must be there; see data/README.md\n"
             "every other step of the pipeline runs without them.")

io.banner("04  IoU sensitivity over voxel size, expansion radius and threshold")
clouds = load_all(clouds_dir, TOOLS)
rng = np.random.default_rng(C["seed"])
lo1, hi1 = p["ensemble"]["phase1"]; lo2, hi2 = p["ensemble"]["phase2"]
sel = {"Phase1": rng.choice(np.arange(lo1, hi1 + 1), C["n_frames_sampled"], replace=False),
       "Phase2": rng.choice(np.arange(lo2, hi2 + 1), C["n_frames_sampled"], replace=False)}

rows, ranks = [], {}
for voxel, radius in itertools.product(C["sweep_voxel_A"], C["sweep_radius_A"]):
    for phase in ["Phase1", "Phase2"]:
        pair = pairwise_iou(clouds, sel[phase], TOOLS, voxel, radius)
        v = pair.to_numpy()
        for thr in C["sweep_threshold"]:
            rows.append(dict(voxel_A=voxel, radius_A=radius, phase=phase, iou_threshold=thr,
                             pairs_above=int((v >= thr).sum()), pairs_total=len(v),
                             iou_mean=round(float(v.mean()), 4),
                             iou_median=round(float(np.median(v)), 4),
                             iou_max=round(float(v.max()), 4), top_pair=pair.idxmax()))
        ranks[f"v{voxel}_R{radius}_{phase}"] = pair.rank(ascending=False)
        print(f"  voxel {voxel}  R {radius}  {phase}:  mean {v.mean():.3f}  max {v.max():.3f} ({pair.idxmax()})")

io.save(pd.DataFrame(rows), "04_iou_sensitivity.csv")
R = pd.DataFrame(ranks)
summary = pd.DataFrame({"median_rank": R.median(axis=1), "min_rank": R.min(axis=1),
                        "max_rank": R.max(axis=1), "excursion": R.max(axis=1) - R.min(axis=1)})
summary.index.name = "tool_pair"
io.save(summary.sort_values("median_rank").reset_index(), "04_iou_pair_ranking.csv")
print(summary.sort_values("median_rank").to_string())
