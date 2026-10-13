"""SiteFerret -> output.pdb.N.pdb + centroid table

SiteFerret, run through scripts/loop.py, writes per frame an
output/<name>_pockets directory (older versions: *_Pfiles) containing for each pocket

    p<n>_atm.pqr            the surface atoms lining the pocket
    clusterPocket<n>.pqr    the probe spheres filling it

and the corresponding p<n>_sub*/clusterPocket<n>_sub* for the subpockets, which
are excluded here: the analysis uses main pockets only.

The combined PDB is built from the p<n>_atm.pqr atoms, which is what defines the
pocket envelope. The pocket centre written to the centroid table is the centre
of geometry of those same atoms; the centres of the probe spheres are written to
a second table for comparison, since the two differ by a few angstrom and the
choice affects nothing downstream but should be explicit.

    python3 convert_siteferret.py <native_dir> <prepared_dir> <dest_dir> <n_frames>
                                  [--centroids centroids_SiteFerret.csv]
"""
import os, sys, glob, re

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import protein_lines, write_pdb, points_from_pqr, summarise, write_centroids

argv = [a for a in sys.argv[1:] if not a.startswith("--")]
NAT, PREP, DEST, N = argv[0], argv[1], argv[2], int(argv[3])
OUT = (sys.argv[sys.argv.index("--centroids") + 1]
       if "--centroids" in sys.argv else os.path.join(DEST, "centroids_SiteFerret.csv"))
os.makedirs(DEST, exist_ok=True)


def index(path, pattern):
    m = re.search(pattern, os.path.basename(path))
    return int(m.group(1)) if m else 0


lining, probes, done = [], [], 0
for i in range(1, N + 1):
    base = os.path.join(NAT, "f%d" % i)
    if not os.path.isdir(base):
        continue
    pfiles = glob.glob(os.path.join(base, "*_Pfiles")) or glob.glob(os.path.join(base, "*_pockets"))
    if not pfiles:
        continue
    pfiles = pfiles[0]

    atm = [p for p in glob.glob(os.path.join(pfiles, "p*_atm.pqr"))
           if "sub" not in os.path.basename(p)]
    clu = [p for p in glob.glob(os.path.join(pfiles, "clusterPocket*.pqr"))
           if "_sub" not in os.path.basename(p)]
    if not atm and not clu:
        continue

    groups = []
    for p in sorted(atm, key=lambda s: index(s, r"p(\d+)_atm")):
        pts = points_from_pqr(p)
        if pts:
            groups.append((index(p, r"p(\d+)_atm"), pts))
    if groups:
        write_pdb(os.path.join(DEST, "output.pdb.%d.pdb" % i), protein_lines(PREP, i), groups)
        done += 1

    for k, pts in groups:
        c = np.mean(pts, axis=0)
        lining.append((i, k, c[0], c[1], c[2], len(pts)))

    for p in sorted(clu, key=lambda s: index(s, r"clusterPocket(\d+)")):
        pts = points_from_pqr(p)
        if pts:
            c = np.mean(pts, axis=0)
            probes.append((i, index(p, r"clusterPocket(\d+)"), c[0], c[1], c[2], len(pts)))

print("[siteferret] %d frames converted" % done)
summarise(DEST, "STP")

HEADER = "system,tool,frame,pocket_id,rank,x,y,z,n_points,phase"
row = lambda r, ph: ("GLUT1,SiteFerret,%d,%d,%d,%.3f,%.3f,%.3f,%d,%s"
                     % (r[0], r[1], r[1], r[2], r[3], r[4], r[5], ph))
write_centroids(OUT, lining, "SiteFerret", HEADER, row)
write_centroids(OUT.replace(".csv", "_probe_spheres.csv"), probes, "SiteFerret", HEADER, row)
