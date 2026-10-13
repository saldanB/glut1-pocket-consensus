# Converters

Each detector reports pockets in its own format. These scripts turn one tool's
native per-frame output into the common representation used by everything
downstream: a PDB per frame containing the protein plus the pocket-defining
points, and a centroid table with the columns of `data/pockets_master.csv`.

They are provided so that the step between the detectors and the deposited
pocket table is inspectable. Reproducing the analyses in the paper does not
require them — `data/pockets_master.csv` is deposited.

| script | reads | writes |
|---|---|---|
| `convert_autosite.py` | `f<N>/*_cl_*.pdb` | `output.pdb.N_autosite_merged.pdb` |
| `convert_siteferret.py` | `f<N>/*_Pfiles/p*_atm.pqr` | `output.pdb.N.pdb` + centroids |
| `convert_sitehound.py` | `f<N>/*_CMET_summary.dat`, `*_CMET_clusters.dat` | `output.pdb.N_pockets.pdb` + centroids |
| `convert_deeppocket.py` | `f<N>/*_pocket<N>.pdb` | `output.pdb.N.pdb` |

Fpocket and P2Rank are not listed because their native output is already read
directly (`*_out/pockets/pocket*_vert.pqr` and `*_predictions.csv`).

The four converters take the same arguments:

```bash
python3 convert_<tool>.py <native_dir> <prepared_dir> <dest_dir> <n_frames>
```

where `<native_dir>` contains one directory `f1 ... fN` per frame, and
`<prepared_dir>` contains the input structures `output.pdb.<N>.pdb` that were
given to the detector.

Three choices in these scripts are decisions rather than mechanics, and are the
ones worth checking:

* **SiteFerret** — subpockets are excluded, and the pocket centre is the centre
  of geometry of the lining atoms (`p<n>_atm.pqr`), not of the probe spheres.
  The probe-sphere centres are written to a second table for comparison.
* **SiteHound** — the number of clusters kept per frame is set by the knee of
  the sorted energy distribution (clusters at or below -|knee|), not by a fixed
  cutoff. The pocket centre is the cluster centre of the summary file; the grid
  points of each cluster are written as its point cloud. See the script header.
* **DeepPocket** — only the segmentation output is read, and the pocket centre
  is the mean of the atoms of the segmented pocket. The Fpocket candidates
  that DeepPocket ranks sit in the same working directory and must not be
  collected; see the script header.
