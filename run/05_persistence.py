#!/usr/bin/env python3
"""Sensitivity of the persistence estimates (R4 #1).

  a  minimum run length (10, 20, 30% of the phase) x gap tolerance delta (0, 2, 5, 10)
  b  minimum soft/strict persistence for each min_cluster_size scale of run/03a
  c  phase 2 subsampled to 100 frames (phase 1 min_cluster_size), 40 replicates

Writes results/05_persistence_sensitivity.csv, 05b_persistence_by_mcs.csv,
       05c_phase2_100frames.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
import src.consensus as C

p = io.parameters(); S = p["sensitivity"]
df = io.pockets(aligned=True)
MT = p["consensus"]["min_tools"]


def runs_fraction(min_run):
    """Same definition as src.consensus._runs_fraction, with a variable minimum run length."""
    def f(frames, all_frames, delta):
        all_frames = np.sort(np.asarray(all_frames)); n = all_frames.size
        rank = {v: i for i, v in enumerate(all_frames)}
        r = np.sort(np.unique([rank[v] for v in np.asarray(frames) if v in rank]))
        if r.size == 0:
            return 0.0
        runs, s, q = [], r[0], r[0]
        for x in r[1:]:
            if x - q - 1 <= delta: q = x
            else: runs.append((s, q)); s = q = x
        runs.append((s, q))
        return sum(e - s + 1 for s, e in runs if e - s + 1 >= min_run * n) / n
    return f


io.banner("05a  persistence: minimum run length x delta")
orig, rows = C._runs_fraction, []
for ph in ["Phase1", "Phase2"]:
    for frac in [0.1, 0.2, 0.3]:
        C._runs_fraction = runs_fraction(frac)
        for delta in [0, 2, 5, 10]:
            p2 = dict(p); p2["consensus"] = dict(p["consensus"], delta_frames=delta)
            _, reg = C.run_pipeline(df, ph, p2); s = C.select(reg, MT)
            sp = [r["soft_persistence"] for r in s]
            rows.append(dict(phase=ph, run_min=frac, delta=delta, regions=len(s),
                             pers_min=round(100 * min(sp), 1), pers_median=round(100 * float(np.median(sp)), 1),
                             pers_max=round(100 * max(sp), 1), regions_ge50=sum(x >= 0.5 for x in sp)))
C._runs_fraction = orig
t = pd.DataFrame(rows); io.save(t, "05_persistence_sensitivity.csv"); print(t.to_string(index=False))

io.banner("05b  minimum persistence for each min_cluster_size scale")
rows = []
for ph in ["Phase1", "Phase2"]:
    base = p["min_cluster_size"][ph]
    for s in S["min_cluster_size_scales"]:
        mcs = {t: max(5, int(round(base[t] * s))) for t in p["tools"]}
        _, reg = C.run_pipeline(df, ph, p, min_cluster_size=mcs); sel = C.select(reg, MT)
        rows.append(dict(phase=ph, scale=s, regions=len(sel),
                         soft_min=round(min(r["soft_persistence"] for r in sel), 3),
                         strict_min=round(min(r["strict_persistence"] for r in sel), 3)))
t = pd.DataFrame(rows); io.save(t, "05b_persistence_by_mcs.csv"); print(t.to_string(index=False))

io.banner("05c  phase 2 subsampled to 100 frames")
rng = np.random.default_rng(S["seed"]); fr = np.sort(df[df.phase == "Phase2"].frame.unique()); v = []
for k in range(S["subsample_replicates"]):
    pick = rng.choice(fr, 100, replace=False)
    _, reg = C.run_pipeline(df, "Phase2", p, min_cluster_size=p["min_cluster_size"]["Phase1"], frames=pick)
    v.append(len(C.select(reg, MT)))
v = np.array(v)
io.save(pd.DataFrame({"replicate": range(1, len(v) + 1), "regions_ge4": v}), "05c_phase2_100frames.csv")
print(f"regions >= 4 tools: mean {v.mean():.2f}, sd {v.std(ddof=1):.2f}, range {v.min()}-{v.max()}")
