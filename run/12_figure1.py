#!/usr/bin/env python3
"""Figure 1: HDBSCAN versus DBSCAN on the AutoSite centroids (superposed), both phases.

HDBSCAN: min_cluster_size of data/parameters.yaml. DBSCAN: eps and min_samples of
the `figure1` block of data/parameters.yaml, identical for the two phases. The
figure is an illustrative example of the dependence of DBSCAN on eps; the full
sweep is produced by run/08.

Writes results/12_figure1.png, results/12_figure1.svg, results/12_figure1_counts.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd, hdbscan
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
from src import io
io.banner('12  Figure 1 (HDBSCAN versus DBSCAN, AutoSite)')
p = io.parameters(); f1 = p['figure1']; tool = f1['tool']
m = io.pockets(aligned=True)
plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 11})
fig, axes = plt.subplots(2, 2, figsize=(13.5, 14.5))
rows = []
for i, ph in enumerate(('Phase1', 'Phase2')):
    X = m[(m.tool == tool) & (m.phase == ph)][['x', 'y', 'z']].values
    mcs = p['min_cluster_size'][ph][tool]
    labs = {'HDBSCAN': hdbscan.HDBSCAN(min_cluster_size=mcs).fit_predict(X),
            'DBSCAN': DBSCAN(eps=f1['dbscan_eps_A'], min_samples=f1['dbscan_min_samples']).fit_predict(X)}
    for j, (meth, l) in enumerate(labs.items()):
        ax = axes[i, j]; ids = sorted(set(l) - {-1}); pal = [c for k, c in enumerate(plt.get_cmap('tab10').colors) if k != 7]  # tab10 without its grey, as in the original
        noise = l == -1
        ax.scatter(X[noise, 0], X[noise, 1], s=3, c='#c7c7c7', label='-1', linewidths=0)
        for k, c in enumerate(ids):
            mm = l == c
            ax.scatter(X[mm, 0], X[mm, 1], s=3, color=pal[k % len(pal)], label=str(k), linewidths=0)
        ax.set_xlim(-35, 35); ax.set_ylim(-35, 35); ax.set_xticks(range(-30, 31, 10)); ax.set_yticks(range(-30, 31, 10))
        ax.set_aspect('equal'); ax.grid(True, ls='--', lw=0.6, color='#b0b0b0')
        ax.set_title(f'{meth} clustering - Autosite (phase {i + 1})', fontweight='bold')
        lg = ax.legend(title='Cluster', loc='upper right', markerscale=2.5, fontsize=9, title_fontsize=10,
                       handletextpad=0.2, borderpad=0.4, labelspacing=0.15, framealpha=0.9)
        rows.append(dict(phase=ph, method=meth, n_points=len(X), clusters=len(ids),
                         noise_pct=round(100 * noise.mean(), 1),
                         parameters=(f'min_cluster_size={mcs}' if meth == 'HDBSCAN' else
                                     f"eps={f1['dbscan_eps_A']} A, min_samples={f1['dbscan_min_samples']}")))
    axes[i, 0].text(-0.08, 1.05, 'AB'[i], transform=axes[i, 0].transAxes, fontsize=15, fontweight='bold')
fig.tight_layout(h_pad=3)
for ext in ('png', 'svg'):
    fig.savefig(os.path.join(io.RESULTS, f'12_figure1.{ext}'), dpi=300 if ext == 'png' else None)
t = pd.DataFrame(rows); io.save(t, '12_figure1_counts.csv'); print(t.to_string(index=False))
