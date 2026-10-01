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

Both of the following are in the deposit named by the `deposit` block of
`parameters.yaml`.

**Set `deposit.doi` yourself before fetching them.** It ships as `null`, because
the DOI is minted only when the deposit is published, so the download cannot work
until you put it there:

```yaml
deposit:
  doi: 10.5281/zenodo.1234567    # the DOI in the Data Availability statement
```

Use the **version** DOI, not the concept DOI. The concept DOI always resolves to
the newest version of the deposit, so a reproduction run against it would silently
fetch whatever has been uploaded since; the version DOI pins the files the paper
was written from. `deposit.api` is Zenodo's real endpoint and needs no change.

Then

```bash
python run/00_fetch.py            # the point clouds
python run/00_fetch.py --frames   # and the frames as well, see deposit.frames_archive
```

which downloads into the right place, checks each file against the checksum the
deposit publishes for it, and skips what is already present. Downloading by hand
works too; the filenames below are what the code looks for.

**Point clouds** (`pointclouds/pc_GLUT1_<tool>.npz`, 31 MB) hold the
pocket-defining elements of each pocket — alpha-sphere centres, energy grid
points or atomic coordinates depending on the tool — and are needed only by
`run/04_iou.py`. All six must be present, with exactly these names: the tool
name is capitalised as in `parameters.yaml` and the system prefix is `GLUT1`.

**The 500 frames** themselves (PDB) are needed only to re-run the detectors,
which is documented in `tools/`.
