#!/usr/bin/env python3
"""Superpose the ensemble and quantify what the superposition does (Editor E3, R3 #4).

Reads  data/pockets_master.csv, data/protein_ca.npz
Writes results/pockets_aligned.csv, results/01_alignment_effect.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.align import superpose_pockets, displacement

p = io.parameters()
io.banner("01  rigid-body superposition on protein CA, reference = first frame")
raw = io.pockets(aligned=False)
ca, frames = io.protein_ca()
aligned = superpose_pockets(raw, ca, frames, p["alignment"]["reference_frame"])
io.save(aligned, "pockets_aligned.csv")

rows = []
d = displacement(raw, aligned)
phase = raw.set_index(["tool", "frame", "pocket_id"]).loc[
    aligned.set_index(["tool", "frame", "pocket_id"]).index, "phase"].to_numpy()
cut = p["consensus"]["cutoff_A"]
for name, mask in [("all", np.ones(len(d), bool)), ("Phase1", phase == "Phase1"),
                   ("Phase2", phase == "Phase2")]:
    v = d[mask]
    rows.append(dict(subset=name, n=len(v),
                     median_A=round(float(np.median(v)), 2), mean_A=round(float(v.mean()), 2),
                     p90_A=round(float(np.percentile(v, 90)), 2),
                     p99_A=round(float(np.percentile(v, 99)), 2),
                     max_A=round(float(v.max()), 2),
                     beyond_cutoff=int((v > cut).sum())))
t = pd.DataFrame(rows)
io.save(t, "01_alignment_effect.csv")
print(t.to_string(index=False))
