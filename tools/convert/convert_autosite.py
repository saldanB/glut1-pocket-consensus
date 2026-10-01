"""AutoSite -> output.pdb.N_autosite_merged.pdb

AutoSite (ADFRsuite) writes one file per cluster. They are merged into a single
PDB per frame in which the pocket points are ATOM records with residue name STP
and a residue number equal to the cluster index.

    python3 convert_autosite.py <native_dir> <prepared_dir> <dest_dir> <n_frames>

<native_dir> holds one subdirectory f1 ... fN per frame, as produced by AS.py.
"""
import os, sys, glob, re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import protein_lines, write_pdb, points_from_pdb, summarise

NAT, PREP, DEST, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
os.makedirs(DEST, exist_ok=True)

# Depending on the ADFRsuite version the cluster files are named
# <name>_cl_001.pdb or <name>_CL_1.pdb. The *_fp_* files are the feature points
# and must not be collected.
RX = re.compile(r"_cl_0*(\d+)\.pdb$", re.IGNORECASE)

done = 0
for i in range(1, N + 1):
    base = os.path.join(NAT, "f%d" % i)
    if not os.path.isdir(base):
        continue
    clusters = [f for f in glob.glob(os.path.join(base, "**", "*.pdb"), recursive=True)
                if RX.search(os.path.basename(f))]
    if not clusters:
        continue

    groups = []
    for p in sorted(clusters, key=lambda s: int(RX.search(os.path.basename(s)).group(1))):
        pts = points_from_pdb(p)
        if pts:
            groups.append((int(RX.search(os.path.basename(p)).group(1)), pts))
    if not groups:
        continue

    write_pdb(os.path.join(DEST, "output.pdb.%d_autosite_merged.pdb" % i),
              protein_lines(PREP, i), groups, rec="ATOM")
    done += 1

print("[autosite] %d frames converted" % done)
summarise(DEST, "STP", "(deposited files: about 2300 per frame)")
