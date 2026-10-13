#!/usr/bin/env python3
"""DBSCAN versus HDBSCAN cluster counts, per tool and phase (Results, Table S16).

DBSCAN is run with min_samples equal to the HDBSCAN min_cluster_size of
data/parameters.yaml and eps from 2 to 6 A, on the superposed centroids.
D<eps> = number of clusters, Dn<eps> = % of points labelled as noise.

Writes results/08_dbscan_vs_hdbscan.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import pandas as pd, hdbscan
from sklearn.cluster import DBSCAN
from src import io
io.banner('08  DBSCAN versus HDBSCAN')
p = io.parameters(); m = io.pockets(aligned=True); rows = []
for ph in ('Phase1', 'Phase2'):
    for t in p['tools']:
        X = m[(m.tool == t) & (m.phase == ph)][['x', 'y', 'z']].values; mcs = p['min_cluster_size'][ph][t]
        l = hdbscan.HDBSCAN(min_cluster_size=mcs).fit_predict(X)
        r = dict(phase=ph, tool=t, n=len(X), mcs=mcs, H=len(set(l)) - (-1 in l), H_noise=round((l == -1).mean() * 100, 1))
        for e in (2, 3, 4, 5, 6):
            d = DBSCAN(eps=e, min_samples=mcs).fit_predict(X)
            r[f'D{e}'] = len(set(d)) - (-1 in d); r[f'Dn{e}'] = round((d == -1).mean() * 100, 1)
        rows.append(r)
t = pd.DataFrame(rows); io.save(t, '08_dbscan_vs_hdbscan.csv'); print(t.to_string(index=False))
