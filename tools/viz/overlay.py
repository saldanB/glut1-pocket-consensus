"""Build the PDB strings the viewer hands to NGL.

Everything is produced as a string rather than a file: NGL takes a structure from
a blob, so the frame, the pocket pseudo-atoms and the consensus regions can be
loaded as separate components without anything being written to disk.

The deposited frames and the centroids in data/pockets_master.csv are both in the
unsuperposed space of the frame as supplied to the detectors, so they overlay
directly. The consensus regions of run/02 are in the superposed space of the
reference frame instead, so asking for them superposes the frame and its pockets
first; `aligned` carries that choice through every builder here.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from src import io
from src.align import kabsch
from src.consensus import cluster_labels

# one chain and one residue name per detector, so that NGL can select them
CHAIN = {"AutoSite": "A", "Fpocket": "B", "P2Rank": "C",
         "SiteHound": "D", "SiteFerret": "E", "DeepPocket": "F"}
RESNAME = {"AutoSite": "AUT", "Fpocket": "FPO", "P2Rank": "P2R",
           "SiteHound": "SHO", "SiteFerret": "SFE", "DeepPocket": "DPO"}
CONSENSUS_CHAIN, CONSENSUS_RESNAME = "Z", "CON"
PUBLISHED_CHAIN, PUBLISHED_RESNAME = "Y", "PUB"
CLUSTER_RESNAME = "CLU"


def _atom(serial: int, name: str, resname: str, chain: str, resseq: int,
          xyz, bfactor: float = 0.0, occupancy: float = 1.0) -> str:
    # fixed-column PDB: name 13-16, resName 18-20, chainID 22, resSeq 23-26,
    # x/y/z 31-54, occupancy 55-60, tempFactor 61-66, element 77-78
    return ("HETATM%5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s\n"
            % (serial % 100000, name[:4], resname[:3], chain, resseq % 10000,
               xyz[0], xyz[1], xyz[2], occupancy, bfactor, "C"))


def transform(frame: int, reference_frame: int = 1):
    """The superposition of `frame` onto the reference, as (centre, rot, ref_centre)."""
    ca, frames = io.protein_ca()
    index = {int(f): i for i, f in enumerate(frames)}
    for f in (frame, reference_frame):
        if f not in index:
            raise SystemExit("frame %d is not in data/protein_ca.npz" % f)
    ref = ca[index[reference_frame]]
    ref_centre = ref.mean(axis=0)
    this = ca[index[frame]]
    centre = this.mean(axis=0)
    return centre, kabsch(this - centre, ref - ref_centre), ref_centre


def apply(xyz: np.ndarray, tf) -> np.ndarray:
    centre, rot, ref_centre = tf
    return (np.asarray(xyz, dtype=float) - centre) @ rot + ref_centre


def structure_pdb(path: str, tf=None) -> str:
    """The frame as it will be displayed, superposed when `tf` is given."""
    lines = []
    for line in open(path, errors="replace"):
        if not line.startswith(("ATOM", "HETATM")):
            if line.startswith(("TER", "END")):
                lines.append(line)
            continue
        if tf is None:
            lines.append(line)
            continue
        xyz = apply([[float(line[30:38]), float(line[38:46]), float(line[46:54])]], tf)[0]
        lines.append("%s%8.3f%8.3f%8.3f%s" % (line[:30], xyz[0], xyz[1], xyz[2], line[54:]))
    return "".join(lines)


def pockets(frame: int, tools=None, top_k: int | None = None,
            tf=None) -> tuple[str, pd.DataFrame]:
    """One pseudo-atom per predicted pocket of this frame."""
    table = pd.read_csv(os.path.join(io.DATA, "pockets_master.csv"))
    sub = table[table.frame == frame]
    if tools:
        sub = sub[sub.tool.isin(tools)]
    if top_k:
        sub = sub[sub["rank"] <= top_k]
    sub = sub.sort_values(["tool", "rank"]).reset_index(drop=True)
    xyz = sub[["x", "y", "z"]].to_numpy(dtype=float)
    if tf is not None and len(xyz):
        xyz = apply(xyz, tf)
    out = []
    for i, row in sub.iterrows():
        out.append(_atom(i + 1, "PKT", RESNAME[row.tool], CHAIN[row.tool], i + 1,
                         xyz[i], bfactor=float(row.n_points), occupancy=float(row["rank"])))
    moved = sub.copy()
    if len(xyz):
        moved[["x", "y", "z"]] = xyz
    return "".join(out) + "END\n", moved


_CLOUDS: dict = {}


def _cloud_by_pocket(tool: str) -> dict[int, list[tuple[str, np.ndarray]]]:
    """frame -> [(pocket id, its pocket-defining points)], read once per tool."""
    if tool in _CLOUDS:
        return _CLOUDS[tool]
    path = os.path.join(io.DATA, "pointclouds", "pc_GLUT1_%s.npz" % tool)
    if not os.path.exists(path):
        raise SystemExit("no point cloud for %s in %s"
                         % (tool, os.path.join(io.DATA, "pointclouds")))
    archive = np.load(path, allow_pickle=True)
    points, counts = archive["points"], archive["counts"]
    frames, pids = archive["frames"], archive["pids"]
    offsets = np.concatenate([[0], np.cumsum(counts)])
    out: dict[int, list[tuple[str, np.ndarray]]] = {}
    for i in range(len(counts)):
        out.setdefault(int(frames[i]), []).append(
            (str(pids[i]), points[offsets[i]:offsets[i + 1]]))
    _CLOUDS[tool] = out
    return out


def cluster_assignment(phase: str, tools) -> tuple[dict, dict]:
    """(tool, frame, pocket id) -> cluster number for one phase, and the counts.

    The clustering is `src.consensus.cluster_labels` on the superposed centroids,
    the same call `run/02` makes with the same ``min_cluster_size``, so a pocket is
    assigned to the cluster the pipeline assigned it to. Cluster numbers start at
    1 per tool; 0 means HDBSCAN called the pocket noise.
    """
    table = io.pockets(aligned=True)
    sizes = io.parameters()["min_cluster_size"][phase]
    mapping, counts = {}, {}
    for tool in tools:
        sub = table[(table.phase == phase) & (table.tool == tool)]
        xyz = sub[["x", "y", "z"]].to_numpy(dtype=float)
        minimum = sizes[tool]
        labels = (cluster_labels(xyz, minimum) if len(xyz) >= minimum
                  else np.full(len(xyz), -1))
        numbering = {label: i for i, label in enumerate(sorted(set(labels) - {-1}), 1)}
        for frame, pocket, label in zip(sub.frame.astype(int),
                                        sub.pocket_id.astype(str), labels):
            mapping[(tool, int(frame), pocket)] = numbering.get(label, 0)
        counts[tool] = len(numbering)
    return mapping, counts


def clouds(frame: int, tools, tf=None, assignment: dict | None = None) -> tuple[str, int]:
    """The pocket-defining points of this frame's pockets, as (pdb, n_points).

    With `assignment`, each point carries the cluster number of the pocket it
    belongs to in its residue number, so the viewer can colour by cluster; without
    it the residue number is the pocket's own index.
    """
    lines, total = [], 0
    for tool in tools:
        by_frame = _cloud_by_pocket(tool)
        for index, (pocket, points) in enumerate(by_frame.get(int(frame), []), 1):
            if tf is not None and len(points):
                points = apply(points, tf)
            resseq = (assignment.get((tool, int(frame), pocket), 0)
                      if assignment is not None else index)
            for point in points:
                total += 1
                lines.append(_atom(total, "PNT", RESNAME[tool], CHAIN[tool],
                                   resseq, point))
    return "".join(lines) + "END\n", total


def hulls(frame: int, tools, tf=None, assignment: dict | None = None,
          colours: dict | None = None) -> list[dict]:
    """The convex hull of each cluster's points in this frame, as flat-shaded meshes.

    The points of every pocket of this frame that belongs to the same cluster are
    pooled and their convex hull taken, so the surface is a property of the frame
    and moves with it. A cluster needs four points in general position to have a
    hull at all; one that does not -- a single pocket of a detector that reports one
    point per pocket, or a set that happens to be coplanar -- is skipped.

    Each mesh is flat shaded: a triangle carries its own three vertices and its own
    face normal, so the facets of the hull stay visible instead of being rounded
    off. Triangles are wound so that their normals point away from the hull centre.
    """
    from scipy.spatial import ConvexHull
    from scipy.spatial.qhull import QhullError  # type: ignore[import]

    meshes = []
    for tool in tools:
        by_frame = _cloud_by_pocket(tool)
        groups: dict[int, list[np.ndarray]] = {}
        for pocket, points in by_frame.get(int(frame), []):
            number = (assignment.get((tool, int(frame), pocket), 0)
                      if assignment is not None else 0)
            if number == 0:                      # noise has no cluster to enclose
                continue
            groups.setdefault(number, []).append(points)
        for number, parts in sorted(groups.items()):
            points = np.concatenate(parts).astype(float)
            if tf is not None:
                points = apply(points, tf)
            if len(points) < 4:
                continue
            try:
                hull = ConvexHull(points)
            except (QhullError, ValueError):      # coplanar or otherwise degenerate
                continue
            vertices = points[hull.vertices]
            centre = vertices.mean(axis=0)
            remap = {v: i for i, v in enumerate(hull.vertices)}
            position, normal = [], []
            for simplex in hull.simplices:
                a, b, c = (vertices[remap[s]] for s in simplex)
                face = np.cross(b - a, c - a)
                length = np.linalg.norm(face)
                if length == 0:
                    continue
                face = face / length
                if face @ ((a + b + c) / 3.0 - centre) < 0:   # wind it outwards
                    b, c = c, b
                    face = -face
                for vertex in (a, b, c):
                    position.extend(round(float(x), 2) for x in vertex)
                    normal.extend(round(float(x), 3) for x in face)
            if not position:
                continue
            meshes.append(dict(
                tool=tool, cluster=int(number),
                colour=(colours or {}).get((tool, number), "#9aa4b2"),
                n_points=int(len(points)), n_triangles=len(position) // 9,
                volume=round(float(hull.volume), 1),
                position=position, normal=normal,
                index=list(range(len(position) // 3))))
    return meshes


def consensus(phase: str) -> tuple[str, pd.DataFrame]:
    """The consensus regions of run/02, which live in the superposed space."""
    path = os.path.join(io.RESULTS, "02_consensus_regions.csv")
    if not os.path.exists(path):
        raise SystemExit("no %s — run: python run/02_consensus.py" % path)
    table = pd.read_csv(path)
    sub = table[(table.phase == phase) & (table.scenario == "with superposition")]
    sub = sub.sort_values("region").reset_index(drop=True)
    out = []
    for i, row in sub.iterrows():
        out.append(_atom(i + 1, "CEN", CONSENSUS_RESNAME, CONSENSUS_CHAIN, int(row.region),
                         (row.x, row.y, row.z),
                         bfactor=float(row.n_tools), occupancy=float(row.occupancy)))
    return "".join(out) + "END\n", sub


def published(phase: str) -> tuple[str, pd.DataFrame]:
    """The consensus centres reported in the paper, which ship in data/.

    Like the regions of run/02 these are in the superposed space, so they are
    drawn together with a superposed frame.
    """
    table = pd.read_csv(os.path.join(io.DATA, "published_consensus_centres.csv"))
    sub = table[table.phase == phase].sort_values("centre").reset_index(drop=True)
    out = []
    for i, row in sub.iterrows():
        out.append(_atom(i + 1, "CEN", PUBLISHED_RESNAME, PUBLISHED_CHAIN, int(row.centre),
                         (row.x, row.y, row.z), bfactor=float(row.supporting_tools)))
    return "".join(out) + "END\n", sub


def phase_of(frame: int) -> str:
    bounds = io.parameters()["ensemble"]
    lo, hi = bounds["phase1"]
    return "Phase1" if lo <= frame <= hi else "Phase2"
