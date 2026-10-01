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
    check("centroids superposed", int(r.loc["all", "n"]), 43826)
    check("median displacement (A)", float(r.loc["all", "median_A"]), 1.10)
    check("maximum displacement (A)", float(r.loc["all", "max_A"]), 5.20)
    check("centroids beyond the 6.5 A cutoff", int(r.loc["all", "beyond_cutoff"]), 0)
    check("median displacement, phase 1", float(r.loc["Phase1", "median_A"]), 0.57)
    check("median displacement, phase 2", float(r.loc["Phase2", "median_A"]), 1.31)

d = need("02_consensus_summary.csv")
if d is not None:
    a = d[d.scenario == "with superposition"].set_index("phase")
    b = d[d.scenario == "without superposition"].set_index("phase")
    check("consensus regions, phase 1", int(a.loc["Phase1", "regions_ge4"]), 11)
    check("consensus regions, phase 2", int(a.loc["Phase2", "regions_ge4"]), 13)
    check("regions with >=5 tools, phase 1", int(a.loc["Phase1", "regions_ge5"]), 7)
    check("regions with >=5 tools, phase 2", int(a.loc["Phase2", "regions_ge5"]), 7)
    check("published centres recovered, phase 1", int(a.loc["Phase1", "recovered_within_cutoff"]), 9)
    check("published centres recovered, phase 2", int(a.loc["Phase2", "recovered_within_cutoff"]), 10)
    check("phase 2 regions without superposition", int(b.loc["Phase2", "regions_ge4"]), 10)
    check("phase 2 centres recovered without superposition",
          int(b.loc["Phase2", "recovered_within_cutoff"]), 8)

d = need("03a_min_cluster_size.csv")
if d is not None:
    r = d.set_index(["phase", "scale"])
    for ph, want in [("Phase1", [14, 12, 11, 9, 7, 5]), ("Phase2", [16, 16, 13, 10, 8, 4])]:
        for s, w in zip([0.5, 0.75, 1.0, 1.25, 1.5, 2.0], want):
            check(f"regions at {s}x, {ph}", int(r.loc[(ph, s), "regions_ge4"]), w, tol=1)
    for ph, want in [("Phase1", [9, 9, 9, 7, 5, 4]), ("Phase2", [11, 11, 10, 8, 7, 3])]:
        for s, w in zip([0.5, 0.75, 1.0, 1.25, 1.5, 2.0], want):
            check(f"centres recovered at {s}x, {ph}", int(r.loc[(ph, s), "recovered_cutoff"]), w, tol=1)

d = need("03b_agreement_threshold.csv")
if d is not None:
    r = d.set_index(["phase", "threshold"])
    for ph, want in [("Phase1", [13, 11, 7, 3]), ("Phase2", [17, 13, 7, 5])]:
        for k, w in zip(["3/6", "4/6", "5/6", "6/6"], want):
            check(f"regions at {k}, {ph}", int(r.loc[(ph, k), "n_regions"]), w, tol=1)
    check("minimum occupancy at 3/6, phase 2",
          round(float(r.loc[("Phase2", "3/6"), "occupancy_min"]), 2), 0.71, tol=0.01)

d = need("03c_leave_one_tool_out.csv")
if d is not None:
    r = d.set_index(["phase", "excluded"])
    for tool in ["AutoSite", "Fpocket", "P2Rank", "SiteHound", "SiteFerret", "DeepPocket"]:
        check(f"phase 1, excluding {tool}: recovered at the equivalent threshold",
              r.loc[("Phase1", tool), "fair_recovered_3A"], "11/11")
    check("phase 2, excluding DeepPocket: recovered within 3 A",
          r.loc[("Phase2", "DeepPocket"), "fair_recovered_3A"], "12/13")
    check("phase 2, excluding DeepPocket: recovered within the cutoff",
          r.loc[("Phase2", "DeepPocket"), "fair_recovered_cutoff"], "13/13")

d = need("03e_matched_subsampling_test.csv")
if d is not None:
    r = d.set_index("metric")
    if int(r.loc["regions_ge4", "n_replicates"]) >= 40:
        check("matched subsampling, phase 1 mean", float(r.loc["regions_ge4", "mean_Phase1"]), 9.20, tol=0.3)
        check("matched subsampling, phase 2 mean", float(r.loc["regions_ge4", "mean_Phase2"]), 5.60, tol=0.5)
    else:
        print("  SKIP  matched subsampling run with --quick, means not comparable")

print(f"\n  {ok} passed, {warn} within tolerance, {bad} failed")
sys.exit(1 if bad else 0)
