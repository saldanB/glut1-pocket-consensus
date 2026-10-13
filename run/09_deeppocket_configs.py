#!/usr/bin/env python3
"""Consensus with alternative DeepPocket configurations.

For each table in data/deeppocket_configs/ (frame, pocket, x, y, z of the
segmented pockets, centroid = mean of the pocket atoms) the DeepPocket pockets
are replaced and the consensus is recomputed, everything else unchanged.
seg0_<classifier>_t<threshold>_m<mask_dist>: phase 1 only, segmentation model
seg0_best_test_IOU_91; supplementary_cmdline_*: the command line of the original
Supplementary, all 500 frames.

Writes results/09_deeppocket_configs.csv
"""
import os, sys, glob
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.align import superpose_pockets
from src.consensus import run_pipeline, select, recovery
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(io.DATA, 'deeppocket_configs')
p = io.parameters(); raw = io.pockets(aligned=False); ca, frames = io.protein_ca()
ref = pd.read_csv(os.path.join(io.RESULTS, '02_consensus_regions.csv')); ref = ref[ref.scenario == 'with superposition']
rows = []
for f in sorted(glob.glob(f'{D}/*.csv')):
    name = os.path.basename(f)[:-4]
    dp = pd.read_csv(f)
    dp = pd.DataFrame(dict(system='GLUT1', tool='DeepPocket', frame=dp.frame, pocket_id=dp.pocket, rank=dp.pocket,
                           x=dp.x, y=dp.y, z=dp.z, n_points=np.nan, score=np.nan,
                           phase=np.where(dp.frame <= 100, 'Phase1', 'Phase2')))
    data = pd.concat([raw[raw.tool != 'DeepPocket'], dp], ignore_index=True)
    al = superpose_pockets(data, ca, frames, p['alignment']['reference_frame'])
    for ph in ('Phase1', 'Phase2'):
        if (dp.phase == ph).sum() == 0: continue
        cl, reg = run_pipeline(al, ph, p); sel = select(reg, 4)
        r = ref[ref.phase == ph][['x', 'y', 'z']].values
        d = recovery(r, sel)
        rows.append(dict(config=name, phase=ph, dp_pockets=int((dp.phase == ph).sum()),
                         dp_per_frame=round((dp.phase == ph).sum() / dp[dp.phase == ph].frame.nunique(), 2),
                         dp_clusters=len(cl['DeepPocket']), regions=len(sel),
                         regions_with_dp=sum('DeepPocket' in str(x.get('tools')) for x in sel),
                         final_regions_recovered=f'{int((d <= 6.5).sum())}/{len(r)}', median_dist=round(float(np.median(d)), 2)))
        print(rows[-1], flush=True)
io.save(pd.DataFrame(rows), '09_deeppocket_configs.csv')
