# Data

## In this repository

| file | what it is |
|---|---|
| `pockets_master.csv` | every pocket predicted by the six detectors on the 500 GLUT1 frames: tool, frame, pocket id, rank, centroid, number of pocket-defining points, phase. 43,405 rows (8,732 in phase 1, 34,673 in phase 2). This is the input to the whole analysis. |
| `protein_ca.npz` | protein Calpha coordinates per frame (500 x 447 x 3) plus the frame numbers, used for the rigid-body superposition. |
| `reference_frame1.pdb` | frame 1 (protein and D-glucose; the substrate is ignored by `run/11_table4.py`, which lists the protein residues lining each region). |
| `residue_centres.npz` | geometric centre of every protein residue in every frame (`prot`, 500 x 447 x 3), used by `run/00_correlation.py` for Table 2. |
| `published_consensus_centres.csv` | the consensus centres of the original submission (9 for phase 1, 11 for phase 2), used as a reference in the recovery analyses. |
| `pointclouds/pc_GLUT1_<tool>.npz` | the pocket-defining elements of every pocket, depending on the tool (see below). One archive per detector. Read by `run/04`, `06` and `07`. |
| `deeppocket_configs/*.csv` | DeepPocket pockets (frame, pocket, centroid) under eight alternative configurations on phase 1 and under the command line of the original Supplementary on all 500 frames. Read only by `run/09`. |
| `parameters.yaml` | every parameter of the pipeline. |

The D-glucose substrate was excluded from the input of every detector. The
coordinates in `pockets_master.csv` are those of the frames as supplied to
the detectors, i.e. **not** superposed. `run/01_align.py` produces the superposed
table; the difference is quantified in the paper.

Pockets per tool:

| tool | phase 1 | phase 2 | centroid |
|---|---|---|---|
| Fpocket | 2,989 | 12,356 | mean of the alpha-sphere centres |
| AutoSite | 1,556 | 5,640 | mean of the cluster points (`*_cl_*.pdb`) |
| SiteHound | 1,228 | 4,702 | cluster centre from `*_CMET_summary.dat` |
| P2Rank | 1,078 | 4,386 | mean of the SAS points (pocket 0 excluded) |
| DeepPocket | 881 | 3,589 | mean of the atoms of the segmented pocket |
| SiteFerret | 1,000 | 4,000 | mean of the lining atoms (`p<n>_atm.pqr`) |

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
describe the same objects: the pocket counts agree per tool and in total, and
`(frame, pid)` identifies the same pocket in both.

What the points are:

| tool | points |
|---|---|
| Fpocket | alpha-sphere centres |
| SiteFerret | atoms lining the pocket |
| AutoSite | cluster grid points |
| SiteHound | grid points of the cluster (`*_CMET_clusters.dat`) |
| P2Rank | solvent-accessible surface points |
| DeepPocket | atoms of the segmented pocket |

For every tool except SiteHound the mean of a pocket's points is its centroid in
the table. For SiteHound the centroid is the cluster centre reported in the
summary file, which lies 0.15 A (median) from the mean of the grid points.

## Not in this repository

**The 500 frames** themselves (PDB) are not included: they are needed only to
re-run the detectors or the MDpocket cross-check, both documented in `tools/`,
and every analysis in `run/` works from the files above. They are deposited at
doi:[10.5281/zenodo.21238931](https://doi.org/10.5281/zenodo.21238931) as
`frames-GLUT1_MD.zip`, which unpacks to `frames/frame.<N>.pdb` for N = 1 to 500.
