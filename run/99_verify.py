#!/usr/bin/env python3
"""Check the regenerated tables against the values reported in the paper.

Exit status 0 when every check passes. Checks marked TOLERANT allow a small
deviation: HDBSCAN cluster counts can differ by one between library builds at
borderline densities, which propagates to a few derived counts. Everything that
is build-independent is checked exactly.
"""
import os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import pandas as pd
from src import io

R = io.RESULTS
ok = warn = bad = 0


def check(name, got, want, tol=0):
    global ok, warn, bad
    if got == want:
        print(f"  PASS  {name}: {got}"); ok += 1
    elif tol and abs(got - want) <= tol:
        print(f"  TOL   {name}: {got} (paper: {want}, within +/-{tol})"); warn += 1
    else:
        print(f"  FAIL  {name}: {got}, paper says {want}"); bad += 1


def need(name):
    path = os.path.join(R, name)
    if not os.path.exists(path):
        print(f"  SKIP  {name} not produced yet"); return None
    return pd.read_csv(path)

io.banner("99  verification against the published values")

d = need("01_alignment_effect.csv")
if d is not None:
    r = d.set_index("subset")
    check("centroids superposed", int(r.loc["all", "n"]), 43405)
    check("median displacement (A)", float(r.loc["all", "median_A"]), 1.11)
    check("maximum displacement (A)", float(r.loc["all", "max_A"]), 5.20)
    check("centroids beyond the 6.5 A cutoff", int(r.loc["all", "beyond_cutoff"]), 0)
    check("median displacement, phase 1", float(r.loc["Phase1", "median_A"]), 0.58)
    check("median displacement, phase 2", float(r.loc["Phase2", "median_A"]), 1.32)

d = need("02_consensus_summary.csv")
if d is not None:
    a = d[d.scenario == "with superposition"].set_index("phase")
    b = d[d.scenario == "without superposition"].set_index("phase")
    check("consensus regions, phase 1", int(a.loc["Phase1", "regions_ge4"]), 10)
    check("consensus regions, phase 2", int(a.loc["Phase2", "regions_ge4"]), 11)
    check("regions with >=5 tools, phase 1", int(a.loc["Phase1", "regions_ge5"]), 6)
    check("regions with >=5 tools, phase 2", int(a.loc["Phase2", "regions_ge5"]), 9)
    check("regions with 6 tools, phase 1", int(a.loc["Phase1", "regions_ge6"]), 2)
    check("regions with 6 tools, phase 2", int(a.loc["Phase2", "regions_ge6"]), 4)
    check("original centres recovered, phase 1", int(a.loc["Phase1", "recovered_within_cutoff"]), 8)
    check("original centres recovered, phase 2", int(a.loc["Phase2", "recovered_within_cutoff"]), 9)
    check("regions without superposition, phase 1", int(b.loc["Phase1", "regions_ge4"]), 9)
    check("regions without superposition, phase 2", int(b.loc["Phase2", "regions_ge4"]), 8)
    check("phase 2 original centres recovered without superposition",
          int(b.loc["Phase2", "recovered_within_cutoff"]), 7)

d = need("03a_min_cluster_size.csv")
if d is not None:
    r = d.set_index(["phase", "scale"])
    for ph, want in [("Phase1", [11, 9, 10, 10, 6, 2]), ("Phase2", [17, 13, 11, 11, 8, 6])]:
        for s, w in zip([0.5, 0.75, 1.0, 1.25, 1.5, 2.0], want):
            check(f"regions at {s}x, {ph}", int(r.loc[(ph, s), "regions_ge4"]), w, tol=1)
    for ph, want in [("Phase1", [8, 7, 8, 8, 6, 2]), ("Phase2", [11, 10, 9, 9, 8, 5])]:
        for s, w in zip([0.5, 0.75, 1.0, 1.25, 1.5, 2.0], want):
            check(f"original centres recovered at {s}x, {ph}", int(r.loc[(ph, s), "recovered_cutoff"]), w, tol=1)

d = need("03b_agreement_threshold.csv")
if d is not None:
    r = d.set_index(["phase", "threshold"])
    for ph, want in [("Phase1", [13, 10, 6, 2]), ("Phase2", [18, 11, 9, 4])]:
        for k, w in zip(["3/6", "4/6", "5/6", "6/6"], want):
            check(f"regions at {k}, {ph}", int(r.loc[(ph, k), "n_regions"]), w, tol=1)
    check("minimum occupancy at 3/6, phase 2",
          round(float(r.loc[("Phase2", "3/6"), "occupancy_min"]), 2), 0.72, tol=0.01)
    check("minimum occupancy at 4/6, phase 2",
          round(float(r.loc[("Phase2", "4/6"), "occupancy_min"]), 2), 0.89, tol=0.01)

d = need("03c_leave_one_tool_out.csv")
if d is not None:
    r = d.set_index(["phase", "excluded"])
    for tool in ["AutoSite", "Fpocket", "P2Rank", "SiteHound", "SiteFerret", "DeepPocket"]:
        check(f"phase 1, excluding {tool}: recovered at the equivalent threshold",
              r.loc[("Phase1", tool), "fair_recovered_3A"], "10/10")
    for tool in ["AutoSite", "Fpocket", "P2Rank", "SiteFerret", "DeepPocket"]:
        check(f"phase 2, excluding {tool}: recovered at the equivalent threshold",
              r.loc[("Phase2", tool), "fair_recovered_3A"], "11/11")
    check("phase 2, excluding SiteHound: recovered within 3 A",
          r.loc[("Phase2", "SiteHound"), "fair_recovered_3A"], "10/11")
    check("phase 2, excluding SiteHound: recovered within the cutoff",
          r.loc[("Phase2", "SiteHound"), "fair_recovered_cutoff"], "11/11")

d = need("03d_top_k_retention.csv")
if d is not None:
    r = d[d.min_cluster_size == "rescaled"].set_index(["phase", "K"])
    for (ph, k), w in {("Phase1", "15"): "9/10", ("Phase2", "15"): "10/11",
                       ("Phase1", "10"): "8/10", ("Phase2", "10"): "10/11",
                       ("Phase1", "5"): "5/10", ("Phase2", "5"): "7/11"}.items():
        check(f"top-{k} retention, {ph}", r.loc[(ph, k), "recovered_cutoff"], w)

d = need("03e_matched_subsampling_test.csv")
if d is not None:
    r = d.set_index("metric")
    if int(r.loc["regions_ge4", "n_replicates"]) >= 40:
        check("matched subsampling >=4 tools, phase 1 mean", float(r.loc["regions_ge4", "mean_Phase1"]), 9.50, tol=0.3)
        check("matched subsampling >=4 tools, phase 2 mean", float(r.loc["regions_ge4", "mean_Phase2"]), 5.45, tol=0.3)
        check("matched subsampling >=5 tools, phase 1 mean", float(r.loc["regions_ge5", "mean_Phase1"]), 4.53, tol=0.3)
        check("matched subsampling >=5 tools, phase 2 mean", float(r.loc["regions_ge5", "mean_Phase2"]), 2.50, tol=0.3)
    else:
        print("  SKIP  matched subsampling run with --quick, means not comparable")

d = need("04_iou_sensitivity.csv")
if d is not None:
    m = d[(d.voxel_A == 1.0) & (d.radius_A == 5.0)]
    check("maximum pairwise IoU at the manuscript parameters", round(float(m.iou_max.max()), 3), 0.458)
    check("maximum pairwise IoU over the sweep", round(float(d.iou_max.max()), 3), 0.519)
    for ph, w in [("Phase1", [15, 15, 11, 0]), ("Phase2", [15, 15, 12, 0])]:
        q = m[m.phase == ph].set_index("iou_threshold")
        for thr, ww in zip([0.1, 0.2, 0.3, 0.5], w):
            check(f"pairs above IoU {thr}, {ph}", int(q.loc[thr, "pairs_above"]), ww)

d = need("05_persistence_sensitivity.csv")
if d is not None:
    r = d[d.delta == 5].set_index(["phase", "run_min"])
    check("minimum soft persistence, phase 1 (run 20%)", float(r.loc[("Phase1", 0.2), "pers_min"]), 85.5)
    for frac, w in [(0.1, 50.8), (0.2, 43.9), (0.3, 33.8)]:
        check(f"minimum soft persistence, phase 2 (run {int(frac*100)}%)", float(r.loc[("Phase2", frac), "pers_min"]), w)

d = need("05c_phase2_100frames.csv")
if d is not None:
    check("phase 2 at 100 frames, mean regions", round(float(d.regions_ge4.mean()), 2), 8.15, tol=0.3)

d = need("06_membrane_fraction.csv")
if d is not None:
    r = d.set_index("tool")
    for t, w in dict(SiteHound=33.2, P2Rank=42.4, Fpocket=37.9, AutoSite=33.5, SiteFerret=9.4, DeepPocket=19.8).items():
        check(f"membrane-band fraction, {t} (%)", float(r.loc[t, "membrane_pct"]), w)

d = need("06_membrane_consensus.csv")
if d is not None:
    r = d.set_index(["scenario", "phase"])
    check("regions without membrane pockets, phase 1", int(r.loc[("membrane pockets removed", "Phase1"), "regions_ge4"]), 5, tol=1)
    check("regions without membrane pockets, phase 2", int(r.loc[("membrane pockets removed", "Phase2"), "regions_ge4"]), 5, tol=1)

d = need("07_iou_vs_hdbscan.csv")
if d is not None:
    for ph, w in [("Phase1", 10), ("Phase2", 10)]:
        check(f"HDBSCAN regions with an IoU trajectory, {ph}", int(d[(d.phase == ph)].iou_occ.notna().sum()), w)

d = need("09_deeppocket_configs.csv")
if d is not None:
    r = d[d.config.str.startswith("seg0_")]
    check("DeepPocket alternative configurations recovering 10/10 phase 1 regions",
          int((r.final_regions_recovered == "10/10").sum()), 8)

d = need("10_pairwise_rmsd.csv")
if d is not None:
    r = d.set_index("phase")
    check("pairwise Calpha RMSD, phase 1 (A)", float(r.loc["Phase1", "mean_pairwise_rmsd_A"]), 0.79)
    check("pairwise Calpha RMSD, phase 2 (A)", float(r.loc["Phase2", "mean_pairwise_rmsd_A"]), 1.42)

d = need("00_correlation.csv")
if d is not None:
    r = d.set_index(["tool_a", "tool_b"]).pearson_r
    for k, w in {("DeepPocket", "P2Rank"): 0.544, ("AutoSite", "DeepPocket"): 0.357, ("AutoSite", "Fpocket"): 0.356,
                 ("AutoSite", "P2Rank"): 0.354, ("AutoSite", "SiteHound"): 0.334, ("Fpocket", "SiteHound"): 0.096,
                 ("DeepPocket", "Fpocket"): 0.198}.items():
        check(f"Pearson r, {k[0]}-{k[1]}", round(float(r.loc[k]), 3), w)

d = need("12_figure1_counts.csv")
if d is not None:
    r = d.set_index(["phase", "method"]).clusters
    for k, w in {("Phase1", "HDBSCAN"): 14, ("Phase1", "DBSCAN"): 18, ("Phase2", "HDBSCAN"): 18, ("Phase2", "DBSCAN"): 23}.items():
        check(f"Figure 1, {k[1]} clusters, {k[0]}", int(r.loc[k]), w, tol=1)

print(f"\n  {ok} passed, {warn} within tolerance, {bad} failed")
sys.exit(1 if bad else 0)
