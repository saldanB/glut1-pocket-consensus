#!/usr/bin/env python3
"""Phase 1 versus phase 2 on frame-level quantities (R3 #2, R4 #2).

  a  pockets per frame and radius of gyration of the pocket centroids, per tool,
     compared with a block permutation test (block = autocorrelation time)
  b  mean pairwise Calpha RMSD after optimal superposition, at matched sample
     size (80 frames, same replicates as run/03e)

Writes results/10_per_frame_test.csv, results/10_pairwise_rmsd.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io

p = io.parameters(); S = p["sensitivity"]
m = io.pockets(aligned=True)
P1 = p["ensemble"]["phase1"]; N = p["ensemble"]["n_frames"]


def tau(x):
    x = x - x.mean(); v = x.var()
    if v == 0: return 1
    s = 1
    for k in range(1, len(x) // 4):
        r = np.dot(x[:-k], x[k:]) / ((len(x) - k) * v)
        if r <= 0: break
        s += 2 * r
    return max(1, int(np.ceil(s)))


io.banner("10a  per-frame block permutation test")
rng = np.random.default_rng(S["seed"]); rows = []; n1 = P1[1]
for t, g in m.groupby("tool"):
    cnt = g.groupby("frame").size().reindex(range(1, N + 1), fill_value=0).values.astype(float)
    rg = g.groupby("frame").apply(lambda d: np.sqrt(((d[["x", "y", "z"]] - d[["x", "y", "z"]].mean()) ** 2).sum(1).mean())).reindex(range(1, N + 1)).values
    for name, x in [("pockets_per_frame", cnt), ("rg_centroids_A", rg)]:
        a, b = x[:n1], x[n1:]; obs = np.nanmean(b) - np.nanmean(a)
        L = tau(np.nan_to_num(x, nan=np.nanmean(x)))
        blocks = [x[i:i + L] for i in range(0, N, L)]; nb1 = int(round(n1 / L)); ext = 0; R = 10000
        for _ in range(R):
            q = rng.permutation(len(blocks))
            A = np.concatenate([blocks[i] for i in q[:nb1]]); B = np.concatenate([blocks[i] for i in q[nb1:]])
            ext += abs(np.nanmean(B) - np.nanmean(A)) >= abs(obs)
        rows.append(dict(tool=t, metric=name, phase1=round(np.nanmean(a), 2), phase2=round(np.nanmean(b), 2),
                         diff=round(obs, 2), block_len=L, p=round((ext + 1) / (R + 1), 4)))
t = pd.DataFrame(rows); io.save(t, "10_per_frame_test.csv"); print(t.to_string(index=False))

io.banner("10b  pairwise Calpha RMSD at matched sample size")
ca, frames = io.protein_ca()


def rmsd(a, b):
    a = a - a.mean(0); b = b - b.mean(0); U, s, Vt = np.linalg.svd(a.T @ b)
    s[-1] *= np.sign(np.linalg.det(U @ Vt))
    return np.sqrt(max(0.0, ((a ** 2).sum() + (b ** 2).sum() - 2 * s.sum()) / len(a)))


rng = np.random.default_rng(S["seed"]); rows = []
for ph, idx in [("Phase1", np.arange(0, n1)), ("Phase2", np.arange(n1, N))]:
    v = []
    for r in range(S["subsample_replicates"]):
        s = rng.choice(idx, S["subsample_n_frames"], replace=False)
        v.append(np.mean([rmsd(ca[i], ca[j]) for k, i in enumerate(s) for j in s[k + 1:]]))
    rows.append(dict(phase=ph, n_frames=S["subsample_n_frames"], replicates=len(v),
                     mean_pairwise_rmsd_A=round(float(np.mean(v)), 2), sd=round(float(np.std(v, ddof=1)), 3)))
t = pd.DataFrame(rows); io.save(t, "10_pairwise_rmsd.csv"); print(t.to_string(index=False))
