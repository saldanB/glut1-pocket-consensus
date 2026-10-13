"""DeepPocket -> output.pdb.N.pdb

    python3 convert_deeppocket.py <native_dir> <prepared_dir> <dest_dir> <n_frames>

Which files are the segmented pockets
-------------------------------------
At the end of a run, predict.py calls segment_pockets.test(), which writes, for
each ranked pocket,

    <protein>_nowat_pocket<rank>.pdb

containing the residues of the predicted site. When the segmentation network
predicts no residue for a candidate, that file is not created at all, so the
indices are not contiguous: with -r 15 a frame typically yields 7 to 11 files
numbered non-consecutively. With the configuration used in the paper
(coach420 models, -r 15 -t 0.5 --mask_dist 3.5) a frame yields about nine.
The pocket centroid in the table is the mean of the atoms of the segmented pocket.

The frame's working directory also contains <protein>_nowat_out/pockets/ with
pocket0_atm.pdb, pocket1_atm.pdb and so on. Those are the *Fpocket* candidates
that DeepPocket ranks, not the segmented output; collecting every file matching
pocket*.pdb would therefore return Fpocket's pocket set (about 31 per frame)
under DeepPocket's name. Only files matching <something>_pocket<N>.pdb at the
frame level are read here, and anything below a *_out/ directory is ignored.
"""
import os, sys, glob, re

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import protein_lines, write_pdb, points_from_pdb, summarise

NAT, PREP, DEST, N = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
os.makedirs(DEST, exist_ok=True)

RX = re.compile(r"_pocket(\d+)\.pdb$", re.IGNORECASE)

done, empty, per_frame = 0, 0, []
for i in range(1, N + 1):
    base = os.path.join(NAT, "f%d" % i)
    if not os.path.isdir(base):
        continue

    segmented = [p for p in glob.glob(os.path.join(base, "*_pocket*.pdb")) if RX.search(p)]
    groups = []
    for p in sorted(segmented, key=lambda s: int(RX.search(s).group(1))):
        pts = points_from_pdb(p)
        if pts:
            groups.append((int(RX.search(p).group(1)), pts))

    if not groups:
        empty += 1
        continue

    write_pdb(os.path.join(DEST, "output.pdb.%d.pdb" % i), protein_lines(PREP, i), groups)
    per_frame.append(len(groups))
    done += 1

print("[deeppocket] %d frames converted" % done)
if empty:
    print("   %d frames with no segmented pocket" % empty)
if per_frame:
    print("   pockets per frame: mean %.2f, min %d, max %d   (deposited table: 8.81 / 8.97 per phase)"
          % (sum(per_frame) / float(len(per_frame)), min(per_frame), max(per_frame)))
summarise(DEST, "STP")
