# Running the six detectors

This directory documents how `data/pockets_master.csv` was produced. Re-running
it is **not** required to reproduce any analysis in the paper: the table is
deposited, and everything downstream of it runs from `run/`.

We document it rather than ship a turnkey installer because two of the six are
no longer reliably obtainable — EasyMIFs is not currently distributed, and the
DeepPocket checkpoints are behind a link that has been dead at least once. A
script that pretended otherwise would break.

## Input preparation

Each frame is stripped of the POPC bilayer, waters and neutralising ions. One
D-glucose molecule (residue `GLC 448`) is present in every frame; it is excluded
from the analysis, so that every detector sees the protein alone. For the three
detectors that would otherwise read it as part of the receptor it is removed from
the PDB first (`grep -v " GLC " frame.pdb > frame_noGLC.pdb`); the other three
exclude it through their own input handling:

| tool | input | D-glucose |
|---|---|---|
| Fpocket | PDB (glucose written as `HETATM`) | ignored by Fpocket |
| P2Rank | PDB without glucose | removed (predictions are identical with and without it) |
| AutoSite | PDBQT without glucose, via `prepare_receptor4.py -A none` | removed |
| EasyMIFs/SiteHound | PDB without glucose, read by Open Babel inside `Easymifs.py` | removed |
| DeepPocket | PDB | removed with the other non-standard residues before candidate generation |
| SiteFerret | PQR, via PDB2PQR with the PARSE force field | dropped by PDB2PQR |

## Command lines, one frame at a time

```bash
# Fpocket 4
fpocket -f frame.pdb

# P2Rank 2.5.1   (pocket 0 is discarded downstream)
prank predict -f frame.pdb -o out/

# AutoSite  (ADFRsuite)
pythonsh prepare_receptor4.py -r frame.pdb -A none -o frame.pdbqt
autosite -r frame.pdbqt -o out/

# SiteFerret   (top 10 pockets per structure, through the package's loop script)
pdb2pqr30 --ff=PARSE --keep-chain frame.pdb structures/frame.pqr
python scripts/loop.py structures

# EasyMIFs + SiteHound   (Easymifs.py: Open Babel 2.4.1, CMET probe,
#                         grid spacing 0.5 A, box = protein + 5 A)
python Easymifs.py frame.pdb CMET            # writes input_CMET.dx
sitehound -f=input_CMET.dx -t=easymifs -e=-8.9 -l=average -s=7.8

# DeepPocket   (needs fpocket on PATH; it generates the candidates)
python predict.py -c first_model_coach420_best_test_auc_91001.pth.tar \
                  -s coach420_best_test_IOU_44.pth.tar \
                  -p frame.pdb -r 15 -t 0.5 --mask_dist 3.5
```

DeepPocket is run with its default segmentation threshold (0.5) and
residue-to-mask distance (3.5 A). The only non-default setting is the number of
ranked candidates (15 instead of 1), needed to obtain more than one pocket per
structure. The same command line reproduces the deposited Aldose Reductase
predictions. Eight alternative configurations (two classification models,
threshold 0.01 or 0.5, mask distance 3.5 or 5.0 A, segmentation model
`seg0_best_test_IOU_91`) and the command line of the original Supplementary are
in `data/deeppocket_configs/`; `run/09_deeppocket_configs.py` shows that none of
them changes the phase 1 consensus.

## Post-processing

Each tool's native output is converted to the common centroid table by the
converters in [`convert/`](convert/README.md), which also document the three
non-obvious choices in that step. For SiteHound the per-frame energy cutoff is
chosen at the knee of the sorted energy distribution
(`kneed.KneeLocator(curve="convex", direction="increasing")`) rather than fixed,
because the absolute energies are not comparable between frames.

The consensus analysis uses all predictions. Predictions falling entirely within
the membrane band are removed only for visualisation; the filter is strongly
tool-dependent and is quantified by `run/06_membrane.py`, because it turns out to
be a first-order determinant of apparent inter-tool disagreement.

## MDpocket

The trajectory-based cross-check reported in the revised manuscript is in
[`mdpocket/`](mdpocket/README.md): input preparation, the de novo pass, and the
residue-level comparison against the consensus regions. It runs on the
trajectory rather than on the deposited pocket table.

## Versions

Recorded in `environment.yml` for the analysis environment. The detectors
themselves were run at the versions given in Table 1 of the paper.
