#!/usr/bin/env python3
"""Consensus regions, with and without superposition.

Writes results/02_consensus_summary.csv, results/02_consensus_regions.csv
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.consensus import run_pipeline, select, recovery, regions_to_frame

p = io.parameters()
io.banner("02  consensus regions")
aligned = io.pockets(aligned=True)
raw = io.pockets(aligned=False)
pub = io.published_centres()
cut, mt = p["consensus"]["cutoff_A"], p["consensus"]["min_tools"]

summary, regions_out = [], []
for label, data in [("without superposition", raw), ("with superposition", aligned)]:
    for phase in ["Phase1", "Phase2"]:
        clusters, regions = run_pipeline(data, phase, p)
        sel = select(regions, mt)
        d = recovery(pub[phase], sel)
        summary.append(dict(
            scenario=label, phase=phase,
            **{f"clusters_{t}": len(clusters[t]) for t in p["tools"]},
            regions_ge4=len(sel),
            regions_ge5=sum(r["n_tools"] >= 5 for r in sel),
            regions_ge6=sum(r["n_tools"] == 6 for r in sel),
            published_centres=len(d),
            recovered_within_cutoff=int((d <= cut).sum()),
            median_distance_A=round(float(np.median(d)), 2)))
        regions_out.append(regions_to_frame(sel, phase, label))

t = pd.DataFrame(summary)
io.save(t, "02_consensus_summary.csv")
io.save(pd.concat(regions_out, ignore_index=True), "02_consensus_regions.csv")
print(t[["scenario", "phase", "regions_ge4", "regions_ge5", "regions_ge6",
         "recovered_within_cutoff", "published_centres", "median_distance_A"]].to_string(index=False))
