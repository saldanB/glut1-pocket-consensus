#!/usr/bin/env python3
"""Robustness of the consensus regions.

Five analyses, each writing one table:

  a  min_cluster_size rescaled 0.5x - 2x, with recovery of the published centres
  b  inter-tool agreement threshold from 3/6 to 6/6
  c  leave-one-tool-out, at the strict and the proportionally equivalent threshold
  d  harmonised output retention, top K = 5 / 10 / 15
  e  matched-sample-size subsampling of the two phases

Usage:  python run/03_sensitivity.py [a b c d e]   (default: all)
        --quick   fewer subsampling replicates, for a smoke test
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io
from src.consensus import run_pipeline, select, recovery

p = io.parameters()
TOOLS, PHASES = p["tools"], ["Phase1", "Phase2"]
CUT, MT = p["consensus"]["cutoff_A"], p["consensus"]["min_tools"]
S = p["sensitivity"]
args = [a for a in sys.argv[1:] if not a.startswith("-")] or list("abcde")
quick = "--quick" in sys.argv
rng = np.random.default_rng(S["seed"])

df = io.pockets(aligned=True)
pub = io.published_centres()
ref = {}
for ph in PHASES:
    _, reg = run_pipeline(df, ph, p)
    ref[ph] = select(reg, MT)


def frac(d, t):
    return int((d <= t).sum())


# ---------------------------------------------------------------- a --------
if "a" in args:
    io.banner("03a  min_cluster_size, 0.5x to 2x")
    rows = []
    for ph in PHASES:
        base = p["min_cluster_size"][ph]
        for s in S["min_cluster_size_scales"]:
            mcs = {t: max(5, int(round(base[t] * s))) for t in TOOLS}
            _, reg = run_pipeline(df, ph, p, min_cluster_size=mcs)
            sel = select(reg, MT)
            d = recovery(pub[ph], sel)
            rows.append(dict(phase=ph, scale=s,
                             min_cluster_size="|".join(f"{t}:{mcs[t]}" for t in TOOLS),
                             regions_ge4=len(sel),
                             regions_ge5=sum(r["n_tools"] >= 5 for r in sel),
                             regions_ge6=sum(r["n_tools"] == 6 for r in sel),
                             published_centres=len(d),
                             recovered_3A=frac(d, 3.0), recovered_cutoff=frac(d, CUT)))
    t = pd.DataFrame(rows); io.save(t, "03a_min_cluster_size.csv")
    print(t[["phase", "scale", "regions_ge4", "regions_ge5", "regions_ge6",
             "recovered_cutoff", "published_centres"]].to_string(index=False))

# ---------------------------------------------------------------- b --------
if "b" in args:
    io.banner("03b  inter-tool agreement threshold, 3/6 to 6/6")
    rows = []
    for ph in PHASES:
        _, reg = run_pipeline(df, ph, p)
        for k in S["consensus_thresholds"]:
            sel = select(reg, k)
            occ = [r["occupancy"] for r in sel] or [0]
            sp = [r["soft_persistence"] for r in sel] or [0]
            rows.append(dict(phase=ph, threshold=f"{k}/6", n_regions=len(sel),
                             occupancy_median=round(float(np.median(occ)), 3),
                             occupancy_min=round(float(np.min(occ)), 3),
                             soft_persistence_min=round(float(np.min(sp)), 3),
                             regions_below_0_9=int(sum(o < 0.9 for o in occ))))
    t = pd.DataFrame(rows); io.save(t, "03b_agreement_threshold.csv")
    print(t.to_string(index=False))

# ---------------------------------------------------------------- c --------
if "c" in args:
    io.banner("03c  leave-one-tool-out")
    rows = []
    for ph in PHASES:
        centres = np.array([r["centre"] for r in ref[ph]])
        for drop in [None] + TOOLS:
            keep = TOOLS if drop is None else [t for t in TOOLS if t != drop]
            _, reg = run_pipeline(df, ph, p, tools=keep)
            r = {}
            # with five tools, 4/5 is stricter than 4/6 and 3/5 is the equivalent
            # fraction; both are reported, the second is the fair comparison
            for k in ([MT] if drop is None else [MT, 3]):
                sel = select(reg, k)
                d = recovery(centres, sel)
                r[k] = (len(sel), frac(d, 3.0), frac(d, CUT))
            strict, fair = r[MT], r.get(3, r[MT])
            n = len(centres)
            rows.append(dict(phase=ph, excluded=drop or "(none)", n_tools=len(keep),
                             strict=f"{MT}/6" if drop is None else f"{MT}/5",
                             strict_regions=strict[0],
                             strict_recovered_3A=f"{strict[1]}/{n}",
                             strict_recovered_cutoff=f"{strict[2]}/{n}",
                             fair=f"{MT}/6" if drop is None else "3/5",
                             fair_regions=fair[0],
                             fair_recovered_3A=f"{fair[1]}/{n}",
                             fair_recovered_cutoff=f"{fair[2]}/{n}"))
    t = pd.DataFrame(rows); io.save(t, "03c_leave_one_tool_out.csv")
    print(t[["phase", "excluded", "strict_regions", "strict_recovered_3A",
             "fair_regions", "fair_recovered_3A"]].to_string(index=False))

# ---------------------------------------------------------------- d --------
if "d" in args:
    io.banner("03d  harmonised output retention, top K")
    rows = []
    for ph in PHASES:
        centres = np.array([r["centre"] for r in ref[ph]])
        total = {t: len(df[(df.phase == ph) & (df.tool == t)]) for t in TOOLS}
        for K in S["top_k"] + [None]:
            sub = df if K is None else df[df["rank"] <= K]
            kept = {t: len(sub[(sub.phase == ph) & (sub.tool == t)]) for t in TOOLS}
            # retaining fewer predictions also changes the effective granularity
            # of the clustering, so the parameter is reported both ways
            variants = {
                "rescaled": {t: max(5, int(round(p["min_cluster_size"][ph][t]
                                                 * kept[t] / max(1, total[t])))) for t in TOOLS},
                "unscaled": dict(p["min_cluster_size"][ph]),
            }
            for name, mcs in variants.items():
                _, reg = run_pipeline(df, ph, p, min_cluster_size=mcs, top_k=K)
                sel = select(reg, MT)
                d = recovery(centres, sel)
                rows.append(dict(phase=ph, K=K or "all", min_cluster_size=name,
                                 pockets=len(sub[sub.phase == ph]),
                                 retained="|".join(f"{t}:{kept[t]/max(1,total[t]):.2f}" for t in TOOLS),
                                 regions_ge4=len(sel),
                                 recovered_3A=f"{frac(d, 3.0)}/{len(centres)}",
                                 recovered_cutoff=f"{frac(d, CUT)}/{len(centres)}"))
    t = pd.DataFrame(rows); io.save(t, "03d_top_k_retention.csv")
    print(t[["phase", "K", "min_cluster_size", "regions_ge4", "recovered_3A"]].to_string(index=False))

# ---------------------------------------------------------------- e --------
if "e" in args:
    io.banner("03e  matched-sample-size subsampling of the two phases")
    n, reps = S["subsample_n_frames"], (5 if quick else S["subsample_replicates"])
    rows = []
    for rep in range(reps):
        for ph in PHASES:
            frames = np.sort(df[df.phase == ph].frame.unique())
            pick = rng.choice(frames, size=min(n, len(frames)), replace=False)
            _, reg = run_pipeline(df, ph, p, min_cluster_size=p["min_cluster_size"]["Phase1"],
                                  frames=pick)
            rows.append(dict(replicate=rep + 1, phase=ph,
                             regions_ge4=len(select(reg, 4)),
                             regions_ge5=len(select(reg, 5))))
    t = pd.DataFrame(rows); io.save(t, "03e_matched_subsampling.csv")
    summary = []
    for metric in ["regions_ge4", "regions_ge5"]:
        a = t[t.phase == "Phase1"][metric].to_numpy()
        b = t[t.phase == "Phase2"][metric].to_numpy()
        obs = b.mean() - a.mean()
        pool = np.concatenate([a, b])
        cnt = sum(abs(np.mean(q[len(a):]) - np.mean(q[:len(a)])) >= abs(obs)
                  for q in (rng.permutation(pool) for _ in range(20000)))
        summary.append(dict(metric=metric, n_replicates=reps, n_frames=n,
                            mean_Phase1=round(float(a.mean()), 2), sd_Phase1=round(float(a.std(ddof=1)), 2),
                            mean_Phase2=round(float(b.mean()), 2), sd_Phase2=round(float(b.std(ddof=1)), 2),
                            difference=round(float(obs), 2),
                            p_permutation=round((cnt + 1) / 20001, 4)))
    s = pd.DataFrame(summary); io.save(s, "03e_matched_subsampling_test.csv")
    print(s.to_string(index=False))
