"""Compare the MDpocket frequency map with the six-tool consensus (R3 #9).

The comparison is made on lining RESIDUES rather than on distances between
points, for two reasons. First, a residue lines the wall of a cavity, not its
interior: on GLUT1 the lining residue closest to the centre of its own pocket
lies at a median of 4.45 A and never below 2.5 A, so a 3 A residue-to-cavity
criterion is geometrically unsatisfiable and measures nothing. Second, residues
make the two systems commensurable and use the same 5 A convention as Table 4.

    python3 compare_mdpocket.py glut1 <mdpocket_out_dir> <topology.pdb> <dest> --table4 <Table4.csv>
    python3 compare_mdpocket.py ar    <mdpocket_out_dir> <topology.pdb> <dest>

Reported, at several frequency thresholds:
  - recovery of the known site (AR: the four residues of the specificity pocket);
  - overlap between our lining residues and MDpocket's (GLUT1, Jaccard index per
    consensus region);
  - the reverse direction: MDpocket regions with no counterpart in the consensus.
"""
import os, sys, argparse

import numpy as np
import pandas as pd
from scipy import ndimage
from scipy.spatial import cKDTree

THRESHOLDS = [0.2, 0.3, 0.5, 0.7, 0.9]
RADIUS = 5.0            # same lining-residue convention as Table 4
AR_SITE = [298, 299, 300, 301]

ap = argparse.ArgumentParser()
ap.add_argument("system", choices=["glut1", "ar"])
ap.add_argument("mdpocket_out")
ap.add_argument("topology")
ap.add_argument("dest")
ap.add_argument("--table4", default=None)
A = ap.parse_args()
os.makedirs(A.dest, exist_ok=True)


def read_dx(path):
    """Minimal OpenDX reader for mdpout_freq.dx (orthogonal grid)."""
    counts = origin = None
    delta, values, in_data = [], [], False
    for l in open(path, errors="replace"):
        if not in_data:
            if l.startswith("object 1"):
                counts = tuple(int(x) for x in l.split()[-3:])
            elif l.startswith("origin"):
                origin = np.array([float(x) for x in l.split()[1:4]])
            elif l.startswith("delta"):
                delta.append([float(x) for x in l.split()[1:4]])
            elif "data follows" in l:
                in_data = True
            continue
        for t in l.split():
            try:
                values.append(float(t))
            except ValueError:
                pass
    grid = np.array(values[:counts[0] * counts[1] * counts[2]]).reshape(counts)
    return grid, origin, np.array([delta[0][0], delta[1][1], delta[2][2]])


# ------------------------------------------------------------------- protein
xyz, labels = [], []
for l in open(A.topology, errors="replace"):
    if l.startswith("ATOM"):
        xyz.append([float(l[30:38]), float(l[38:46]), float(l[46:54])])
        labels.append("%s%s" % (l[17:20].strip(), l[22:26].strip()))
xyz = np.array(xyz)
labels = np.array(labels)
tree = cKDTree(xyz)
print("protein: %d atoms, %d residues" % (len(xyz), len(set(labels))))

grid, origin, step = read_dx(os.path.join(A.mdpocket_out, "mdpout_freq.dx"))
voxel = float(step.prod())
print("grid %s, spacing %s A, maximum frequency %.2f"
      % (grid.shape, np.round(step, 2), grid.max()))


def regions(threshold):
    """Connected components of the frequency map above `threshold`, each with
    its lining residues."""
    mask = grid >= threshold
    if not mask.any():
        return []
    lab, n = ndimage.label(mask, structure=np.ones((3, 3, 3)))
    out = []
    for i in range(1, n + 1):
        sel = lab == i
        idx = np.argwhere(sel)
        P = idx * step + origin
        w = grid[sel]
        lining = set()
        for j in tree.query_ball_point(P, RADIUS):
            lining.update(labels[j])
        out.append(dict(label=i, n_voxels=len(idx), volume_A3=round(len(idx) * voxel, 1),
                        freq_max=round(float(w.max()), 3),
                        freq_mean=round(float(w.mean()), 3),
                        x=round(float((P * w[:, None]).sum(0)[0] / w.sum()), 2),
                        y=round(float((P * w[:, None]).sum(0)[1] / w.sum()), 2),
                        z=round(float((P * w[:, None]).sum(0)[2] / w.sum()), 2),
                        residues=lining, points=P))
    out.sort(key=lambda r: -r["n_voxels"])
    for k, r in enumerate(out, 1):
        r["region"] = k
    return out


def jaccard(a, b):
    return len(a & b) / len(a | b) if (a or b) else 0.0


T4 = pd.read_csv(A.table4) if (A.system == "glut1" and A.table4) else None
rows, per_region, reverse = [], [], []

for s in THRESHOLDS:
    R = regions(s)
    if not R:
        rows.append(dict(threshold=s, n_regions=0))
        continue
    row = dict(threshold=s, n_regions=len(R),
               total_volume_A3=round(sum(r["volume_A3"] for r in R), 1),
               largest_volume_A3=max(r["volume_A3"] for r in R),
               residues_covered=len(set().union(*[r["residues"] for r in R])),
               protein_residues=len(set(labels)))
    row["fraction_residues_covered"] = round(row["residues_covered"] / row["protein_residues"], 3)

    if A.system == "ar":
        target = {e for e in set(labels) if e[3:].isdigit() and int(e[3:]) in AR_SITE}
        found = [e for e in target if any(e in r["residues"] for r in R)]
        row["known_site_residues"] = len(target)
        row["known_site_recovered"] = len(found)
        if s in (0.2, 0.5):
            print("\n--- threshold %.1f: the four residues of the specificity pocket" % s)
            for e in sorted(target, key=lambda x: int(x[3:])):
                q = [r for r in R if e in r["residues"]]
                if q:
                    q = max(q, key=lambda r: r["volume_A3"])
                    print("   %-8s  lines region #%d (%.0f A3, max frequency %.2f)"
                          % (e, q["region"], q["volume_A3"], q["freq_max"]))
                else:
                    print("   %-8s  lines no region" % e)
    else:
        scores = []
        union = set().union(*[r["residues"] for r in R])
        for _, r4 in T4.iterrows():
            ours = {x for x in str(r4.lining_residues).split(";") if x}
            # Two different questions, two different denominators:
            #  - coverage:   do our residues fall inside the MDpocket map at all,
            #                in any of its regions?
            #  - resolution: is there a SINGLE region of its map corresponding to
            #                one of ours?
            best = max(R, key=lambda r: jaccard(ours, r["residues"]))
            j = jaccard(ours, best["residues"])
            scores.append(j)
            per_region.append(dict(threshold=s, id=r4["id"], n_tools=r4.n_tools,
                                   jaccard=round(j, 3),
                                   our_residues=len(ours),
                                   shared_with_best_region=len(ours & best["residues"]),
                                   inside_map_anywhere=len(ours & union),
                                   mdpocket_region=best["region"],
                                   volume_A3=best["volume_A3"]))
        scores = np.array(scores)
        row.update(reference_regions=len(scores),
                   median_jaccard=round(float(np.median(scores)), 3),
                   with_jaccard_ge_0_3=int((scores >= 0.3).sum()),
                   with_jaccard_ge_0_5=int((scores >= 0.5).sum()),
                   with_no_overlap=int((scores < 0.05).sum()))

        # reverse direction: MDpocket regions with no counterpart of ours
        all_ours = [{x for x in str(r4.lining_residues).split(";") if x}
                    for _, r4 in T4.iterrows()]
        orphan = 0
        for r in R:
            jb = max(jaccard(n, r["residues"]) for n in all_ours)
            reverse.append(dict(threshold=s, region=r["region"], volume_A3=r["volume_A3"],
                                freq_max=r["freq_max"], best_jaccard=round(jb, 3)))
            if jb < 0.05:
                orphan += 1
        row["mdpocket_regions_without_counterpart"] = orphan
    rows.append(row)

T = pd.DataFrame(rows)
T.to_csv(os.path.join(A.dest, "mdpocket_%s_comparison.csv" % A.system), index=False)
if per_region:
    pd.DataFrame(per_region).to_csv(
        os.path.join(A.dest, "mdpocket_glut1_per_region.csv"), index=False)
if reverse:
    pd.DataFrame(reverse).to_csv(
        os.path.join(A.dest, "mdpocket_glut1_reverse.csv"), index=False)

pd.set_option("display.width", 250)
print("\n=== %s ===" % A.system.upper())
print(T.to_string(index=False))
print("\nwritten to", A.dest)
