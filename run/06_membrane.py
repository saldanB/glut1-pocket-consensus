#!/usr/bin/env python3
"""Membrane-band retention per tool, and consensus without membrane pockets.

A pocket is membrane-embedded when every one of its pocket-defining points lies
within +/- 15.6 A of the protein's Calpha centre along z (the membrane normal),
in the coordinates of its own frame.

Writes results/06_membrane_fraction.csv, results/06_membrane_consensus.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.consensus import run_pipeline, select, recovery

HALF = 15.6
p = io.parameters()
ca, frames = io.protein_ca()
zc = {int(f): float(ca[i][:, 2].mean()) for i, f in enumerate(frames)}
FAMILY = dict(AutoSite="energy", SiteHound="energy", Fpocket="geometry", SiteFerret="geometry",
              P2Rank="machine learning", DeepPocket="deep learning")

io.banner("06  membrane-band retention")
rows, members = [], set()
for t in p["tools"]:
    z = np.load(os.path.join(io.DATA, "pointclouds", f"pc_GLUT1_{t}.npz"), allow_pickle=True)
    cnt, fr, pid = z["counts"], z["frames"], z["pids"]
    dz = z["points"][:, 2] - np.repeat([zc[int(f)] for f in fr], cnt)
    seg = np.repeat(np.arange(len(cnt)), cnt)
    outside = np.zeros(len(cnt), bool); np.logical_or.at(outside, seg, np.abs(dz) > HALF)
    members |= {(t, int(f), int(k)) for f, k, o in zip(fr, pid, outside) if not o}
    rows.append(dict(tool=t, family=FAMILY[t], pockets=len(cnt), membrane=int((~outside).sum()),
                     membrane_pct=round(100 * (~outside).mean(), 1)))
t = pd.DataFrame(rows); io.save(t, "06_membrane_fraction.csv"); print(t.to_string(index=False))

al = io.pockets(aligned=True)
keep = al[[k not in members for k in zip(al.tool, al.frame, al.pocket_id)]]
ref = pd.read_csv(os.path.join(io.RESULTS, "02_consensus_regions.csv"))
ref = ref[ref.scenario == "with superposition"]
rows = []
for label, data in [("all pockets", al), ("membrane pockets removed", keep)]:
    for ph in ["Phase1", "Phase2"]:
        _, reg = run_pipeline(data, ph, p); s = select(reg, p["consensus"]["min_tools"])
        d = recovery(ref[ref.phase == ph][["x", "y", "z"]].values, s)
        rows.append(dict(scenario=label, phase=ph, regions_ge4=len(s),
                         regions_ge5=sum(r["n_tools"] >= 5 for r in s),
                         reference_regions_recovered=f"{int((d <= p['consensus']['cutoff_A']).sum())}/{len(d)}"))
t = pd.DataFrame(rows); io.save(t, "06_membrane_consensus.csv"); print(t.to_string(index=False))
