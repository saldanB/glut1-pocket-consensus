"""SiteHound / EasyMIFs -> output.pdb.N_pockets.pdb + centroid table

    python3 convert_sitehound.py <native_dir> <prepared_dir> <dest_dir> <n_frames>
                                 [--probe CMET] [--centroids centroids_SiteHound.csv]

Output format: the protein plus one HETATM with residue name POC per retained
cluster (its centre), in rank order, preceded by REMARK records carrying the
interaction energy.

Cluster retention is adaptive rather than a fixed energy cutoff or a fixed
top-N, because the absolute energies returned by EasyMIFs are not comparable
between frames. For each frame the cluster energies (column 2 of
<name>_<probe>_summary.dat) are taken in absolute value and sorted in increasing
order; the knee of that curve is located with

    kneed.KneeLocator(curve="convex", direction="increasing")

and the clusters with energy below -|knee| are retained. On the GLUT1 ensemble
this yields between 2 and 24 clusters per frame (mean 10.61 in phase 1, 9.83 in
phase 2).

Only the CMET probe map is used; the OP map is computed by the pipeline but does
not contribute to the reported pockets.
"""
import os, sys, glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import protein_lines, write_pdb, summarise, write_centroids

try:
    from kneed import KneeLocator
except ImportError:
    sys.exit("this converter requires the kneed package:  pip install kneed")

argv = [a for a in sys.argv[1:] if not a.startswith("--")]
NAT, PREP, DEST, N = argv[0], argv[1], argv[2], int(argv[3])
PROBE = (sys.argv[sys.argv.index("--probe") + 1] if "--probe" in sys.argv
         else os.environ.get("SITEHOUND_PROBE", "CMET"))
OUT = (sys.argv[sys.argv.index("--centroids") + 1]
       if "--centroids" in sys.argv else os.path.join(DEST, "centroids_SiteHound.csv"))
os.makedirs(DEST, exist_ok=True)


def knee_cutoff(energies):
    v = np.sort(np.abs(np.asarray(energies, dtype=float)))
    if len(v) < 3:
        return None
    k = KneeLocator(np.arange(len(v)), v, curve="convex", direction="increasing")
    return None if k.knee is None else float(v[int(k.knee)])


done, rows, per_frame = 0, [], []
for i in range(1, N + 1):
    base = os.path.join(NAT, "f%d" % i)
    if not os.path.isdir(base):
        continue
    summary = (glob.glob(os.path.join(base, "**", "f%d_%s_summary.dat" % (i, PROBE)), recursive=True)
               or glob.glob(os.path.join(base, "**", "*_%s_summary.dat" % PROBE), recursive=True))
    if not summary:
        continue

    data = []
    with open(summary[0], errors="replace") as fh:
        for l in fh:
            c = l.split()
            if len(c) >= 6:
                try:
                    data.append((int(float(c[0])), float(c[1]), float(c[2]),
                                 float(c[3]), float(c[4]), float(c[5])))
                except ValueError:
                    pass
    if not data:
        continue

    cut = knee_cutoff([d[1] for d in data])
    kept = [d for d in data if cut is None or d[1] < -abs(cut)]
    if not kept:
        continue
    kept.sort(key=lambda d: d[1])          # most negative first

    remarks = ["REMARK   --- SITEHOUND POCKET CENTERS ---"]
    groups = []
    for rank, (k, energy, npts, x, y, z) in enumerate(kept, start=1):
        remarks.append("REMARK   POCKET RANK pocket_%d ENERGY: %.3f" % (rank, energy))
        groups.append((rank, [(x, y, z)]))
        rows.append((i, rank, k, x, y, z, energy, npts))

    write_pdb(os.path.join(DEST, "output.pdb.%d_pockets.pdb" % i),
              protein_lines(PREP, i), groups, resn="POC", remarks=remarks)
    per_frame.append(len(kept))
    done += 1

print("[sitehound] %d frames converted (probe %s)" % (done, PROBE))
if per_frame:
    print("   pockets per frame: mean %.2f, min %d, max %d   (deposited: 10.61 / 2 / 24)"
          % (float(np.mean(per_frame)), min(per_frame), max(per_frame)))
summarise(DEST, "POC", "(deposited files: 2-24 per frame)")

write_centroids(
    OUT, rows, "SiteHound",
    "system,tool,frame,rank,cluster_id,x,y,z,energy,n_points,phase",
    lambda r, ph: ("GLUT1,SiteHound,%d,%d,%d,%.3f,%.3f,%.3f,%.3f,%d,%s"
                   % (r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], ph)))
