# Pocket detection across MD ensembles — analysis code

Code and parameters for *Dynamic consensus pocket detection across molecular
dynamics ensembles reveals persistent and transient druggable sites*
([preprint](https://doi.org/10.64898/2026.06.27.734992)).

Everything reported in the paper about clustering, consensus, persistence,
volumetric agreement and the robustness analyses can be regenerated from the
deposited pocket table **without installing any of the six detectors**. That is
deliberate: two of them are no longer straightforward to obtain (EasyMIFs is not
currently distributed, and the DeepPocket checkpoints sit behind a link that has
been dead at least once), so a reproduction path that depends on all six would
not survive. Re-running the detectors themselves is documented separately in
[`tools/`](tools/README.md), with the exact command lines.

## Quick start

```bash
conda env create -f environment.yml
conda activate glut1-pockets

make all          # every step below, then the verification
```

or step by step:

```bash
python run/00_correlation.py
python run/01_align.py
python run/02_consensus.py
python run/03_sensitivity.py
python run/04_iou.py
python run/05_persistence.py
python run/06_membrane.py
python run/07_iou_trajectories.py          # 0.5 / 5 A; the paper also uses 0.3 and 0.7
python run/08_dbscan_vs_hdbscan.py
python run/09_deeppocket_configs.py
python run/10_phase_comparison.py
python run/11_table4.py
python run/12_figure1.py
python run/99_verify.py
```

Every step reads only what is in `data/`: no network access is involved at any
point, and the only thing you have to install is the environment above.

`run/99_verify.py` compares the regenerated tables with the values printed in the
paper and exits non-zero if any of them disagrees. It is the fastest way to tell
whether your environment reproduces our results.

## What each step does

| step | question it answers | paper | output |
|---|---|---|---|
| `00_correlation` | Pearson correlation between tools on residue-to-pocket distances | Table 2 | `00_correlation.csv` |
| `01_align` | how much does the rigid-body superposition move the centroids? | Fig. S7, Table S5 | `01_alignment_effect.csv` |
| `02_consensus` | the consensus regions, with and without superposition | Results, Table 4 | `02_consensus_*.csv` |
| `03_sensitivity a` | does the answer depend on `min_cluster_size`? | Fig. S9, Table S7 | `03a_*.csv` |
| `03_sensitivity b` | on the 4/6 agreement rule? | Fig. S8, Table S8 | `03b_*.csv` |
| `03_sensitivity c` | on any single tool? | Fig. S10, Table S9 | `03c_*.csv` |
| `03_sensitivity d` | on a top-K retention rule (indicative only)? | Fig. S11, Table S11 | `03d_*.csv` |
| `03_sensitivity e` | on the number of frames per phase? | Fig. S6 | `03e_*.csv` |
| `04_iou` | on the voxel size and expansion radius of the IoU criterion? | Fig. 3, Table S10 | `04_iou_*.csv` |
| `05_persistence` | on the minimum run length and gap tolerance of persistence? | Fig. S1, Table S14 | `05*_*.csv` |
| `06_membrane` | how much does the membrane-band filter remove, per tool? | Fig. S12, Table S10 | `06_membrane_*.csv` |
| `07_iou_trajectories` | volumetric trajectories and their match with the HDBSCAN regions | Table 4, Table S15 | `07_*.csv` |
| `08_dbscan_vs_hdbscan` | DBSCAN cluster counts as a function of eps | Table S16 | `08_*.csv` |
| `09_deeppocket_configs` | does the DeepPocket configuration change the consensus? | Table S12 | `09_*.csv` |
| `10_phase_comparison` | per-frame differences between the phases; pairwise RMSD | Fig. S5, Table S6 | `10_*.csv` |
| `11_table4` | Table 4 of the paper | Table 4 | `11_table4*.csv` |
| `12_figure1` | HDBSCAN versus DBSCAN on AutoSite, with declared parameters | Figure 1 | `12_figure1.png`, `12_figure1_counts.csv` |

Individual analyses can be run on their own, for example
`python run/03_sensitivity.py c e`. `--quick` reduces the subsampling
replicates of `03` for a smoke test.

## Looking at the frames

`tools/viz/` draws one frame of the ensemble with the predicted pockets projected
onto it, in the browser via NGL. It is an inspection tool, not part of the
reproduction path: nothing in `run/` depends on it.

```bash
python tools/viz/view_frame.py 42                     # one frame, all six detectors
python tools/viz/view_frame.py 42 --tools Fpocket P2Rank --top-k 5
python tools/viz/view_frame.py 1-24 --out ensemble.html --open
python tools/viz/view_frame.py 1-500:25 --out sweep.html    # every 25th frame
python tools/viz/view_frame.py 42 --clouds                  # pocket-defining points
python tools/viz/view_frame.py 42 --clouds --colour-by cluster   # by HDBSCAN cluster (default)
python tools/viz/view_frame.py 42 --clouds --colour-by tool      # by detector
python tools/viz/view_frame.py 42 --hulls                   # one surface per cluster
python tools/viz/view_frame.py 42 --consensus --published    # regions and paper centres
```

A frame argument is a single frame, an interval (`1-24`), a list (`1,5,9`) or an
interval with a step (`1-500:25`); more than one frame gives a slider. The output
is a standalone HTML page under `results/`.

The frames themselves are not in this repository: the viewer fetches the ones it
needs from the Zenodo deposit into `data/frames/`, which is not tracked. Only the
requested members are read out of the archive, so looking at one frame does not
download the whole ensemble. This is the one part of the project that uses the
network, and nothing in `run/` ever calls it.

`--clouds` and `--hulls` read `data/pointclouds/`. `--colour-by cluster`, the
default, gives every HDBSCAN cluster its own hue, so a cluster keeps its colour as
you step through the frames and the points that belong to no cluster are grey;
cluster numbers restart per detector, so cluster 7 of AutoSite and cluster 7 of
Fpocket are different clusters with different colours. `--colour-by tool` gives one
colour per detector instead, matching the pocket markers. `--color-by` is accepted
as well.

Colouring by cluster is the only option that needs the superposed centroids, so it
implies `--align` and therefore `results/pockets_aligned.csv` from
`run/01_align.py`; `--colour-by tool` does not. `--consensus` needs
`results/02_consensus_regions.csv` from `run/02_consensus.py`. The tool says so if
a table is missing.

See [`tools/viz/README.md`](tools/viz/README.md) for the layers, the colouring and
the superposition.

## Layout

```
data/                    pocket table, Calpha coordinates, residue centres, reference frame,
                         published centres, parameters
data/pointclouds/        per-tool pocket-defining points, read by run/04, 06, 07
data/deeppocket_configs/ DeepPocket pockets under alternative configurations, read by run/09
src/                     the pipeline: consensus.py, align.py, iou.py
run/                     one script per analysis, plus 99_verify.py
tools/                   how the six detectors were run, with the exact command lines
tools/convert/           native tool output -> the common pocket table
tools/mdpocket/          the trajectory-based cross-check
tools/viz/               frame viewer: pockets projected on the structure (NGL)
results/                 produced by the scripts; not tracked
```

**All parameters live in `data/parameters.yaml`.** Nothing is hard-coded in the
code, so an user can change a cutoff and re-run without reading the source.

## A note on reproducibility

HDBSCAN cluster counts can differ by one between library builds where a cluster
sits at the density threshold. `run/99_verify.py` marks the checks that can be
affected `TOL` and treats everything else as exact.

The pocket table itself is the output of the six detectors and is deposited as
data: the numbers in the paper follow from it deterministically, apart from the
caveat above. It was regenerated during the revision by re-running the six
detectors on the deposited frames with the command lines given in `tools/`, so
that the table, the point clouds and the paper describe the same run. The
D-glucose substrate present in the frames was excluded from the input of every
detector (removed from the PDB files for AutoSite, P2Rank and EasyMIFs/SiteHound;
Fpocket, DeepPocket and SiteFerret exclude it through their own input handling).

## Citation

Marigliani G, Petrizzelli F, Mangoni M, Bianco SD, Orzella I, Guzzi PH, Caputo V,
Biagini T, Mazza T. *Dynamic consensus pocket detection across molecular dynamics
ensembles reveals persistent and transient druggable sites*. bioRxiv, 2 July 2026.
doi:[10.64898/2026.06.27.734992](https://doi.org/10.64898/2026.06.27.734992)

`CITATION.cff` carries the same reference in machine-readable form; GitHub renders
it as BibTeX and APA under *Cite this repository*.

The ensemble frames are deposited separately at
doi:[10.5281/zenodo.21238931](https://doi.org/10.5281/zenodo.21238931).
