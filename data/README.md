# Data

## In this repository

| file | what it is |
|---|---|
| `pockets_master.csv` | every pocket predicted by the six detectors on the 500 GLUT1 frames: tool, frame, pocket id, rank, centroid, number of pocket-defining points, phase. 43,826 rows. This is the input to the whole analysis. |
| `protein_ca.npz` | protein Calpha coordinates per frame (500 x 447 x 3) plus the frame numbers, used for the rigid-body superposition. |
| `published_consensus_centres.csv` | the consensus centres reported in the paper (9 for phase 1, 11 for phase 2), used as the reference in the recovery analyses. |
| `pointclouds/pc_GLUT1_<tool>.npz` | the pocket-defining elements of every pocket — alpha-sphere or probe-sphere centres, energy grid points or atomic coordinates, depending on the tool. One archive per detector, 31 MB in total. Read only by `run/04_iou.py`. |
| `parameters.yaml` | every parameter of the pipeline. |

The coordinates in `pockets_master.csv` are those of the frames as supplied to
the detectors, i.e. **not** superposed. `run/01_align.py` produces the superposed
table; the difference is quantified in the paper.

## The point clouds

`run/04_iou.py` reads `pointclouds/` directly: nothing is downloaded and no
network access is involved. All six archives must be present, named exactly

```
pointclouds/pc_GLUT1_AutoSite.npz     pointclouds/pc_GLUT1_SiteHound.npz
pointclouds/pc_GLUT1_Fpocket.npz      pointclouds/pc_GLUT1_SiteFerret.npz
pointclouds/pc_GLUT1_P2Rank.npz       pointclouds/pc_GLUT1_DeepPocket.npz
```

The tool name is capitalised as in `parameters.yaml` and the system prefix is
`GLUT1`; a differently cased name is not found rather than reported, and
`run/04_iou.py` stops if any of the six is missing.

Each archive holds three arrays — `points` (N x 3), `counts` (points per pocket)
and `frames` (frame number per pocket) — plus `pids` and `ranks`, which identify
the pocket each run of points belongs to. The archives and `pockets_master.csv`
describe the same objects: the pocket counts agree per tool and in total, and the
mean of each pocket's points reproduces its centroid in the table.

One asymmetry matters when reading the IoU results. SiteHound contributes exactly
one point per pocket, as its `n_points` column in `pockets_master.csv` also
records, so its occupancy is a single sphere per pocket rather than a cloud. Its
pairs therefore sit at the bottom of the pair ranking partly by construction, not
only because it disagrees with the other detectors.

## Not in this repository

**The 500 frames** themselves (PDB) are not included: they are needed only to
re-run the detectors or the MDpocket cross-check, both documented in `tools/`,
and every analysis in `run/` works from the files above. They are deposited at
doi:[10.5281/zenodo.21238931](https://doi.org/10.5281/zenodo.21238931) as
`frames-GLUT1_MD.zip`, which unpacks to `frames/frame.<N>.pdb` for N = 1 to 500.
