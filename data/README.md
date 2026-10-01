# Data

## In this repository

| file | what it is |
|---|---|
| `pockets_master.csv` | every pocket predicted by the six detectors on the 500 GLUT1 frames: tool, frame, pocket id, rank, centroid, number of pocket-defining points, phase. 43,826 rows. This is the input to the whole analysis. |
| `protein_ca.npz` | protein Calpha coordinates per frame (500 x 447 x 3) plus the frame numbers, used for the rigid-body superposition. |
| `published_consensus_centres.csv` | the consensus centres reported in the paper (9 for phase 1, 11 for phase 2), used as the reference in the recovery analyses. |
| `parameters.yaml` | every parameter of the pipeline. |

The coordinates in `pockets_master.csv` are those of the frames as supplied to
the detectors, i.e. **not** superposed. `run/01_align.py` produces the superposed
table; the difference is quantified in the paper.

## Not in this repository

**Point clouds** (`pointclouds/pc_GLUT1_<tool>.npz`, 31 MB) hold the
pocket-defining elements of each pocket — alpha-sphere centres, energy grid
points or atomic coordinates depending on the tool — and are needed only by
`run/04_iou.py`. Download them from the deposit and put them in
`data/pointclouds/`.

**The 500 frames** themselves (PDB) are in the deposit. They are needed only to
re-run the detectors, which is documented in `tools/`.

Both are at the DOI given in the Data Availability statement of the paper.
