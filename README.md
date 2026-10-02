# Pocket detection across MD ensembles — analysis code

Code and parameters for *From Static Prediction to Dynamic Consensus: Binding
Pocket Detection Across Protein Molecular Dynamics Ensembles*.

Everything reported in the paper about clustering, consensus, persistence,
volumetric agreement and the robustness analyses can be regenerated from the
deposited pocket table **without installing any of the six detectors**. That is
deliberate: two of them are no longer straightforward to obtain (EasyMIFs is not
currently distributed, and the DeepPocket checkpoints sit behind a link that has
been dead at least once), so a reproduction path that depends on all six would
not survive. Re-running the detectors themselves is documented separately in
[`tools/`](tools/README.md).

## Quick start

```bash
conda env create -f environment.yml
conda activate glut1-pockets

python run/01_align.py         # superposition, ~1 min
python run/02_consensus.py     # consensus regions, ~3 min
python run/03_sensitivity.py   # robustness analyses, ~25 min
python run/04_iou.py           # volumetric agreement, ~1 min
python run/99_verify.py        # check the output against the paper
```

or simply `make all`. Every step reads only what is in `data/`: no network access
is involved at any point, and the only thing you have to install is the
environment above.

`run/99_verify.py` compares the regenerated tables with the values printed in the
paper and exits non-zero if any of them disagrees. It is the fastest way to tell
whether your environment reproduces our results.

## What each step does

| step | question it answers | reviewer item | output |
|---|---|---|---|
| `01_align` | how much does the rigid-body superposition move the centroids? | Editor E3, R3 #4 | `01_alignment_effect.csv` |
| `02_consensus` | the consensus regions, with and without superposition | Editor E3, R1 #2 | `02_consensus_*.csv` |
| `03_sensitivity a` | does the answer depend on `min_cluster_size`? | R3 #3, R4 #1 | `03a_*.csv` |
| `03_sensitivity b` | on the 4/6 agreement rule? | R2 #5, R3 #6 | `03b_*.csv` |
| `03_sensitivity c` | on any single tool? | R3 #6, R4 #3 | `03c_*.csv` |
| `03_sensitivity d` | on how much output each tool is allowed to contribute? | Editor E2, R3 #7 | `03d_*.csv` |
| `03_sensitivity e` | on the number of frames per phase? | R3 #2, R4 #2 | `03e_*.csv` |
| `04_iou` | on the voxel size and expansion radius of the IoU criterion? | R3 #8, R4 #1 | `04_iou_*.csv` |

Individual analyses can be run on their own, for example
`python run/03_sensitivity.py c e`. `--quick` reduces the subsampling
replicates for a smoke test.

## Layout

```
data/            pocket table, Calpha coordinates, published centres, parameters
data/pointclouds/
                 per-tool pocket-defining points, read by run/04
src/             the pipeline: consensus.py, align.py, iou.py
run/             one script per analysis, plus 99_verify.py
tools/           how the six detectors were run, with the exact command lines
tools/convert/   native tool output -> the common pocket table
tools/mdpocket/  the trajectory-based cross-check
results/         produced by the scripts; not tracked
```

**All parameters live in `data/parameters.yaml`.** Nothing is hard-coded in the
code, so a reviewer can change a cutoff and re-run without reading the source.

## A note on reproducibility

HDBSCAN cluster counts can differ by one between library builds where a cluster
sits at the density threshold. In our testing this affects one cell of the
`min_cluster_size` sweep and shifts a subsampling mean by less than 0.2 regions;
it does not change any conclusion. `run/99_verify.py` marks those checks
`TOL` and treats everything else as exact.

The pocket table itself is the output of the six detectors and is deposited as
data: the numbers in the paper follow from it deterministically, apart from the
caveat above.
