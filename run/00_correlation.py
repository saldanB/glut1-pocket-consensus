#!/usr/bin/env python3
"""Inter-tool Pearson correlation (Table 2).

For every frame and every tool, the distance from the geometric centre of each
protein residue (data/residue_centres.npz, frames as supplied to the detectors)
to the nearest pocket centroid of that tool in that frame is computed. The
distances of all residues in all 500 frames form one vector per tool, and the
Pearson coefficient is computed for every pair of tools.

Writes results/00_correlation.csv
"""
import os, sys, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from scipy.spatial import cKDTree
from src import io
io.banner('00  inter-tool correlation (Table 2)')
p = io.parameters(); m = io.pockets(aligned=False)
rc = np.load(os.path.join(io.DATA, 'residue_centres.npz'))['prot']          # frames x residues x 3
T = sorted(p['tools']); D = {t: [] for t in T}
for k in range(rc.shape[0]):
    f = k + 1
    for t in T:
        P = m[(m.tool == t) & (m.frame == f)][['x', 'y', 'z']].values
        D[t].append(cKDTree(P).query(rc[k])[0] if len(P) else np.full(rc.shape[1], np.nan))
D = {t: np.concatenate(v) for t, v in D.items()}
rows = []
for a, b in itertools.combinations(T, 2):
    ok = ~np.isnan(D[a]) & ~np.isnan(D[b])
    rows.append(dict(tool_a=a, tool_b=b, pearson_r=round(float(np.corrcoef(D[a][ok], D[b][ok])[0, 1]), 3)))
t = pd.DataFrame(rows); io.save(t, '00_correlation.csv')
print(t.pivot(index='tool_a', columns='tool_b', values='pearson_r').to_string())
