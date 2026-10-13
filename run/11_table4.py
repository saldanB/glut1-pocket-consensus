#!/usr/bin/env python3
"""Table 4 of the paper: consensus regions with identifiers, persistence, IoU
occupancy, lining residues and cross-phase correspondence.

Lining residues: residues with at least one atom within 6 A of the region centre
in the reference frame (data/reference_frame1.pdb), in human GLUT1 numbering
(model numbering + 8). Needs run/02 and run/07.

Writes results/11_table4.csv, and results/11_table4_mdpocket.csv (model numbering,
       input of tools/mdpocket/compare_mdpocket.py --table4)
"""
import os, sys, re
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import numpy as np, pandas as pd
from src import io

OFFSET, CUT = 8, 6.0
R = pd.read_csv(os.path.join(io.RESULTS, "02_consensus_regions.csv")); R = R[R.scenario == "with superposition"]
iou_path = os.path.join(io.RESULTS, "07_iou_vs_hdbscan.csv")
IO = pd.read_csv(iou_path).set_index(["phase", "region"]) if os.path.exists(iou_path) else None
atoms = [(l[17:20].strip(), int(l[22:26]), np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])]))
         for l in open(os.path.join(io.DATA, "reference_frame1.pdb")) if l.startswith(("ATOM", "HETATM"))]
X = np.array([a[2] for a in atoms])
NAME = {"HID": "His", "HIE": "His", "HIP": "His", "CYX": "Cys"}


def lining_model(c):
    """Protein residues in model numbering, as read by tools/mdpocket/compare_mdpocket.py."""
    d = np.linalg.norm(X - c, axis=1)
    return ";".join(sorted({f"{a[0]}{a[1]}" for a, dd in zip(atoms, d) if dd <= CUT and a[0] != "GLC"}))


def lining(c):
    d = np.linalg.norm(X - c, axis=1)
    res = sorted({(a[1], a[0]) for a, dd in zip(atoms, d) if dd <= CUT and a[0] != "GLC"})  # the substrate is not part of the analysed structures
    out = []
    for i, n in res:
        out.append("glucose" if n == "GLC" else f"{NAME.get(n, n.capitalize())}{i + OFFSET}")
    return " ".join(out)


ab = lambda s: s.replace("AutoSite", "AS").replace("DeepPocket", "DP").replace("Fpocket", "Fp") \
    .replace("P2Rank", "P2R").replace("SiteFerret", "SF").replace("SiteHound", "SH").replace("|", ", ")
rows, mdp = [], []
for ph, g in R.groupby("phase"):
    other = R[R.phase != ph]
    for _, r in g.iterrows():
        c = r[["x", "y", "z"]].values.astype(float)
        d = np.linalg.norm(other[["x", "y", "z"]].values - c, axis=1); j = int(d.argmin())
        o = other.iloc[j]
        iou = None if IO is None or (ph, r.region) not in IO.index else IO.loc[(ph, r.region), "iou_occ"]
        mdp.append(dict(id=f"P{ph[-1]}-R{r.region}", n_tools=r.n_tools, lining_residues=lining_model(c)))
        rows.append(dict(region=f"P{ph[-1]}-R{r.region}", tools=ab(r.tools), n_tools=r.n_tools,
                         occupancy=r.occupancy,
                         soft=f"{100 * r.soft_persistence:.1f} ± {100 * r.soft_sd:.1f}",
                         strict=f"{100 * r.strict_persistence:.1f} ± {100 * r.strict_sd:.1f}",
                         iou_occupancy="—" if iou is None or iou != iou else f"{100 * iou:.0f}%",
                         lining=lining(c),
                         other_phase=f"P{o.phase[-1]}-R{o.region} ({d[j]:.1f} Å)" if d[j] <= 6.5 else "—"))
t = pd.DataFrame(rows)
io.banner("11  Table 4")
io.save(t, "11_table4.csv"); io.save(pd.DataFrame(mdp), "11_table4_mdpocket.csv"); print(t.drop(columns=["lining"]).to_string(index=False))
