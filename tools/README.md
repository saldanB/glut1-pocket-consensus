# Running the six detectors

This directory documents how `data/pockets_master.csv` was produced. Re-running
it is **not** required to reproduce any analysis in the paper: the table is
deposited, and everything downstream of it runs from `run/`.

We document it rather than ship a turnkey installer because two of the six are
no longer reliably obtainable — EasyMIFs is not currently distributed, and the
DeepPocket checkpoints are behind a link that has been dead at least once. A
script that pretended otherwise would break.

## Input preparation

Each frame is stripped of the POPC bilayer, waters and neutralising ions, so
that the structure given to the detectors is the protein plus the substrate
where applicable. D-glucose is retained for four of the six tools; the two
exceptions are consequences of the tools themselves, not choices:

* the EasyMIFs/SiteHound pipeline needs a GROMOS-typed structure and cannot
  parametrise the ligand (`prepare_for_easymifs.py` removes `GLC` and renames
  `HID` to `HIS`, both required by `pdb2gmx`);
* DeepPocket's `predict.py` strips heteroatoms as its first step.

Each tool then receives the format it requires:

| tool | input |
|---|---|
| P2Rank, DeepPocket | PDB |
| AutoSite | PDBQT, via `prepare_receptor4.py -A none` |
| Fpocket, SiteFerret, EasyMIFs/SiteHound | PQR, via PDB2PQR with the AMBER force field |

Protonation therefore differs between the PDB- and PQR-derived inputs, which is
stated in the Methods.

## Command lines, one frame at a time

```bash
# Fpocket
fpocket -f frame.pqr

# P2Rank
prank predict -f frame.pdb -o out/

# AutoSite  (ADFRsuite)
pythonsh prepare_receptor4.py -r frame.pdb -A none -o frame.pdbqt
autosite -r frame.pdbqt -o out/

# SiteFerret   (the input filename must contain no dots)
python -m siteFerret frame.pqr

# EasyMIFs + SiteHound
python prepare_pdb.py frame.pdb      # needs atom_types.txt, ffG43b1nb.params, pdb2gmx/
./easymifs ...                       # then SiteHound on the resulting maps

# DeepPocket   (needs fpocket on PATH; it generates the candidates)
python predict.py -c first_model_fold6_best_test_auc_196001.pth.tar \
                  -s seg0_best_test_IOU_91.pth.tar \
                  -p frame.pdb -r 15 -t 0.01 --mask_dist 5.0
```

The three DeepPocket settings that differ from the defaults — 15 ranked pockets
instead of 1, segmentation threshold 0.01 instead of 0.5, residue-to-mask
distance 5.0 Å instead of 3.5 Å — were adopted jointly to operate the tool in a
high-recall regime. Their effect is measured in the paper: the segmentation
threshold was swept over seven values across the whole ensemble and changes the
pocket set only marginally.

## Post-processing

Each tool's native output is converted to the common centroid table by the
converters in [`convert/`](convert/README.md), which also document the three
non-obvious choices in that step. For SiteHound the per-frame energy cutoff is
chosen at the knee of the sorted energy distribution
(`kneed.KneeLocator(curve="convex", direction="increasing")`) rather than fixed,
because the absolute energies are not comparable between frames.

Predictions falling entirely within the membrane band are removed before the
downstream analysis. That filter is strongly tool-dependent and is quantified in
the paper, because it turns out to be a first-order determinant of apparent
inter-tool disagreement.

## MDpocket

The trajectory-based cross-check reported in the revised manuscript is in
[`mdpocket/`](mdpocket/README.md): input preparation, the de novo pass, and the
residue-level comparison against the consensus regions. It runs on the
trajectory rather than on the deposited pocket table.

## Versions

Recorded in `environment.yml` for the analysis environment. The detectors
themselves were run at the versions given in Table 1 of the paper.
