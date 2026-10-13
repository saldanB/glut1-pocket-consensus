#!/usr/bin/env python3
"""Look at one frame of the ensemble with its predicted pockets projected on it.

The frame is taken from data/frames/ and fetched from the deposit if it is not
there yet, so the first look at a frame costs about 0.6 MB and every later one
nothing. Nothing is written into the repository: the structure and the pockets
are handed to NGL as strings.

    python tools/viz/view_frame.py 42
    python tools/viz/view_frame.py 42 --tools Fpocket P2Rank --top-k 5
    python tools/viz/view_frame.py 1-24 --out ensemble.html      # slider over frames
    python tools/viz/view_frame.py 1-500:25 --out sweep.html     # every 25th frame
    python tools/viz/view_frame.py 42 --rotate 0 0 0            # no initial turn
    python tools/viz/view_frame.py 42 --hulls                    # a surface per cluster
    python tools/viz/view_frame.py 42 --clouds                   # points by cluster
    python tools/viz/view_frame.py 42 --clouds --colour-by tool

In a notebook, `view(42)` returns an NGLWidget instead:

    from tools.viz.view_frame import view
    view(42, consensus=True)

--consensus draws the regions of run/02, which live in the superposed space, so
it superposes the frame and its pockets onto the reference frame first; --align
does that on its own.
"""
from __future__ import annotations

import argparse
import collections
import colorsys
import html
import math
import json
import os
import sys
import webbrowser

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from src import io
from tools.viz import frames as frame_cache
from tools.viz import overlay

COLOUR = {"AutoSite": "#4e79a7", "Fpocket": "#f28e2b", "P2Rank": "#59a14f",
          "SiteHound": "#e15759", "SiteFerret": "#b07aa1", "DeepPocket": "#76b7b2"}
CONSENSUS_COLOUR = "#edc948"
PUBLISHED_COLOUR = "#ff9da7"
NOISE_COLOUR = "#6b7280"

# HDBSCAN is rerun once per phase, not once per frame
_ASSIGNMENT: dict = {}


def assignment_for(phase: str, tools) -> tuple[dict, dict]:
    key = (phase, tuple(tools))
    if key not in _ASSIGNMENT:
        _ASSIGNMENT[key] = overlay.cluster_assignment(phase, tools)
    return _ASSIGNMENT[key]


def cluster_colours(phase: str, tools) -> dict[tuple[str, int], str]:
    """(detector, cluster number) -> colour, one hue per cluster of the phase.

    Cluster numbers restart per detector, so the key is the pair: cluster 7 of
    AutoSite and cluster 7 of Fpocket are different clusters. Hues are spaced by the
    golden angle, so no two clusters of a phase come out close together. The points
    and the hulls read the same table, so a cluster is one colour in both.
    """
    _, counts = assignment_for(phase, tools)
    colours, k = {}, 0
    for tool in tools:
        for number in range(1, counts[tool] + 1):
            hue = (k * 137.508 % 360) / 360.0
            k += 1
            r, g, b = colorsys.hls_to_rgb(hue, 0.56, 0.62)
            colours[(tool, number)] = "#%02x%02x%02x" % (int(r * 255), int(g * 255),
                                                        int(b * 255))
    return colours


def cluster_scheme(phase: str, tools) -> list[list[str]]:
    """The same palette as an NGL selection scheme, for the per-atom point colours.

    Residue number 0 is what HDBSCAN called noise, left grey.
    """
    colours = cluster_colours(phase, tools)
    scheme = [[colour, ":%s and %d" % (overlay.CHAIN[tool], number)]
              for (tool, number), colour in colours.items()]
    scheme += [[NOISE_COLOUR, ":%s and 0" % overlay.CHAIN[tool]] for tool in tools]
    return scheme


def parse_frames(text: str) -> list[int]:
    """A frame specification to a list of frame numbers.

    Accepts a single frame (``42``), a closed interval (``1-24``), an interval
    with a step (``1-500:25``), a comma-separated list (``1,5,9``) and any
    combination of those (``1-10,42,200-400:50``). Reversed intervals are read
    in the order given, and duplicates collapse.
    """
    out: list[int] = []
    for piece in (p.strip() for p in text.split(",")):
        if not piece:
            continue
        step = 1
        if ":" in piece:
            piece, step_text = piece.split(":", 1)
            try:
                step = int(step_text)
            except ValueError:
                raise SystemExit("%r is not a step" % step_text)
            if step < 1:
                raise SystemExit("the step must be 1 or more, not %d" % step)
        try:
            if "-" in piece:
                lo_text, hi_text = piece.split("-", 1)
                lo, hi = int(lo_text), int(hi_text)
                if hi < lo:
                    lo, hi = hi, lo
                out.extend(range(lo, hi + 1, step))
            else:
                out.append(int(piece))
        except ValueError:
            raise SystemExit("%r is not a frame or an interval" % piece)
    return sorted(dict.fromkeys(out))


def build(frame: int, tools, top_k=None, clouds=False, consensus=False,
          align=False, published=False, colour_by="cluster", hulls=False,
          quiet=False) -> dict:
    """The strings and the counts for one frame."""
    by_cluster = clouds and colour_by == "cluster"
    aligned = align or consensus or published or by_cluster or hulls
    tf = overlay.transform(frame) if aligned else None
    path = frame_cache.ensure(frame, quiet=quiet)
    out = {"frame": frame,
           "structure": overlay.structure_pdb(path, tf),
           "aligned": aligned,
           "phase": overlay.phase_of(frame)}
    pockets, table = overlay.pockets(frame, tools, top_k, tf)
    out["pockets"] = pockets
    out["per_tool"] = table.tool.value_counts().to_dict()
    out["clouds"], out["n_points"] = "", 0
    if clouds:
        mapping = (assignment_for(out["phase"], tools)[0] if by_cluster else None)
        out["clouds"], out["n_points"] = overlay.clouds(frame, tools, tf, mapping)
    out["hulls"], out["n_hulls"], out["hull_volume"] = "", 0, 0.0
    if hulls:
        mapping = assignment_for(out["phase"], tools)[0]
        meshes = overlay.hulls(frame, tools, tf, mapping,
                              cluster_colours(out["phase"], tools))
        out["hulls"] = json.dumps(meshes, separators=(",", ":"))
        out["n_hulls"] = len(meshes)
        out["hull_volume"] = round(sum(m["volume"] for m in meshes), 1)
        out["hulls_per_tool"] = collections.Counter(m["tool"] for m in meshes)
    out["consensus"] = ""
    out["n_consensus"] = 0
    if consensus:
        regions, table = overlay.consensus(overlay.phase_of(frame))
        out["consensus"], out["n_consensus"] = regions, len(table)
    out["published"] = ""
    out["n_published"] = 0
    if published:
        centres, table = overlay.published(overlay.phase_of(frame))
        out["published"], out["n_published"] = centres, len(table)
    return out


# --------------------------------------------------------------- notebook ----
def view(frame: int, tools=None, top_k=None, clouds=False, consensus=False,
         align=False, published=False, colour_by="cluster",
         rotate=(90.0, 0.0, 0.0)):
    """An NGLWidget for one frame, for use in a notebook."""
    import nglview

    tools = list(tools) if tools else io.parameters()["tools"]
    data = build(frame, tools, top_k, clouds, consensus, align, published,
                 colour_by, quiet=True)
    widget = nglview.NGLWidget()
    widget.add_component(nglview.TextStructure(data["structure"], ext="pdb"))
    widget.clear_representations(component=0)
    widget.add_representation("cartoon", component=0, color="#b8c4d0", opacity=0.45)

    index = 1
    if data["pockets"].count("HETATM"):
        widget.add_component(nglview.TextStructure(data["pockets"], ext="pdb"))
        widget.clear_representations(component=index)
        for tool in tools:
            widget.add_representation("spacefill", component=index,
                                      selection=":" + overlay.CHAIN[tool],
                                      color=COLOUR[tool], radius=1.6, opacity=0.85)
        index += 1
    if data["clouds"]:
        widget.add_component(nglview.TextStructure(data["clouds"], ext="pdb"))
        widget.clear_representations(component=index)
        for tool in tools:
            widget.add_representation("point", component=index,
                                      selection=":" + overlay.CHAIN[tool],
                                      color=COLOUR[tool])
        index += 1
    if data["consensus"]:
        widget.add_component(nglview.TextStructure(data["consensus"], ext="pdb"))
        widget.clear_representations(component=index)
        widget.add_representation("spacefill", component=index,
                                  color=CONSENSUS_COLOUR, radius=2.6, opacity=0.55)
        index += 1
    if data["published"]:
        widget.add_component(nglview.TextStructure(data["published"], ext="pdb"))
        widget.clear_representations(component=index)
        widget.add_representation("spacefill", component=index,
                                  color=PUBLISHED_COLOUR, radius=2.0, opacity=0.6)
    for axis, degrees in zip(([1, 0, 0], [0, 1, 0], [0, 0, 1]), rotate or ()):
        if degrees:
            widget.control.spin(axis, float(degrees) * math.pi / 180.0)
    return widget


# ------------------------------------------------------------------- html ----
PAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "page.html")


def template() -> str:
    """The page that write_html fills in, read from page.html next to this file.

    Kept as its own file rather than a string in here so that an editor treats it
    as the HTML, CSS and JavaScript it is. The placeholders written __LIKE_THIS__
    are substituted by write_html; nothing else in the page is touched, so the
    file is valid HTML on its own and can be opened to check the markup.
    """
    with open(PAGE) as fh:
        return fh.read()


def write_html(data: list[dict], path: str, tools, rotate=(90.0, 0.0, 0.0),
               colour_by="cluster", scheme=None) -> str:
    payload, entries = [], []
    for item in data:
        ids = {}
        for key in ("structure", "pockets", "clouds", "consensus", "published",
                    "hulls"):
            if not item[key]:
                continue
            element = "f%d-%s" % (item["frame"], key)
            ids[key] = element
            payload.append('<script type="text/plain" id="%s">%s</script>'
                           % (element, html.escape(item[key], quote=False)))
        entries.append(dict(frame=item["frame"], phase=item["phase"],
                            aligned=item["aligned"], per_tool=item["per_tool"],
                            n_consensus=item["n_consensus"],
                            n_published=item["n_published"],
                            n_points=item["n_points"],
                            n_hulls=item["n_hulls"], **ids))
    first, last = data[0]["frame"], data[-1]["frame"]
    title = "GLUT1 frame %d" % first if len(data) == 1 else "GLUT1 frames %d-%d" % (first, last)
    subtitle = "%d detector%s, %s" % (len(tools), "" if len(tools) == 1 else "s",
                                      "superposed" if data[0]["aligned"] else "unsuperposed")
    page = (template()
            .replace("__TITLE__", title)
            .replace("__SUBTITLE__", subtitle)
            .replace("__PAYLOAD__", "\n".join(payload))
            .replace("__FRAMES_JSON__", json.dumps(entries))
            .replace("__COLOUR_JSON__", json.dumps({t: COLOUR[t] for t in tools}))
            .replace("__CHAIN_JSON__", json.dumps({t: overlay.CHAIN[t] for t in tools}))
            .replace("__CONSENSUS_COLOUR__", CONSENSUS_COLOUR)
            .replace("__PUBLISHED_COLOUR__", PUBLISHED_COLOUR)
            .replace("__ROTATE_JSON__", json.dumps([float(d) for d in rotate]))
            .replace("__COLOUR_BY__", colour_by)
            .replace("__CLUSTER_SCHEME_JSON__", json.dumps(scheme or {})))
    with open(path, "w") as fh:
        fh.write(page)
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("frames", help="a frame, a range like 1-24, or a list like 1,5,9")
    ap.add_argument("--tools", nargs="+", metavar="TOOL",
                    help="restrict to these detectors (default: all six)")
    ap.add_argument("--top-k", type=int, metavar="K",
                    help="keep only the K best ranked pockets per tool")
    ap.add_argument("--clouds", action="store_true",
                    help="also draw every pocket-defining point, from data/pointclouds")
    ap.add_argument("--consensus", action="store_true",
                    help="also draw the consensus regions of run/02 (implies --align)")
    ap.add_argument("--hulls", action="store_true",
                    help="draw the convex hull of each cluster's points in this "
                         "frame instead of the points themselves (implies --align)")
    ap.add_argument("--colour-by", "--color-by", dest="colour_by",
                    choices=["cluster", "tool"], default="cluster",
                    help="how to colour the points of --clouds: by the HDBSCAN "
                         "cluster of the pocket they belong to (the default, "
                         "implies --align), or by detector")
    ap.add_argument("--published", action="store_true",
                    help="also draw the consensus centres reported in the paper "
                         "(ships in data/, implies --align)")
    ap.add_argument("--align", action="store_true",
                    help="superpose the frame onto the reference frame")
    ap.add_argument("--rotate", nargs=3, type=float, default=[90.0, 0.0, 0.0],
                    metavar=("X", "Y", "Z"),
                    help="degrees to turn the initial view about each axis "
                         "(default: 90 0 0, a quarter turn about x)")
    ap.add_argument("--out", metavar="FILE", help="write a standalone HTML page here")
    ap.add_argument("--open", action="store_true", dest="open_browser",
                    help="open the page in a browser once written")
    A = ap.parse_args()

    known = io.parameters()["tools"]
    tools = A.tools or known
    unknown = [t for t in tools if t not in known]
    if unknown:
        raise SystemExit("unknown tool(s) %s; choose from %s"
                         % (", ".join(unknown), ", ".join(known)))
    wanted = parse_frames(A.frames)
    if not wanted:
        raise SystemExit("no frame given")

    io.banner("viewer  %d frame%s, %d detector%s"
              % (len(wanted), "" if len(wanted) == 1 else "s",
                 len(tools), "" if len(tools) == 1 else "s"))
    if any(A.rotate):
        print("  initial view turned %g/%g/%g degrees about x/y/z" % tuple(A.rotate))
    data = [build(f, tools, A.top_k, A.clouds, A.consensus, A.align, A.published,
                  A.colour_by, A.hulls) for f in wanted]
    scheme = {}
    if A.hulls:
        missing = sorted({t for t in tools
                          if not any(item["hulls_per_tool"].get(t) for item in data)})
        if missing:
            print("  no hull for %s: a cluster needs four points in general position "
                  "in the frame" % ", ".join(missing))
    if A.clouds and A.colour_by == "cluster":
        for phase in sorted({item["phase"] for item in data}):
            scheme[phase] = cluster_scheme(phase, tools)
            _, counts = assignment_for(phase, tools)
            print("  %s: points coloured by %d clusters (%s)"
                  % (phase, sum(counts.values()),
                     " ".join("%s:%d" % (t[:5], counts[t]) for t in tools)))
    for item in data:
        print("  frame %-4d %s  %3d pockets%s  %s"
              % (item["frame"], item["phase"], item["pockets"].count("HETATM"),
                 ("  %6d points" % item["n_points"] if item["n_points"] else "")
                 + ("  %2d hulls %7.0f A3" % (item["n_hulls"], item["hull_volume"])
                    if item["n_hulls"] else ""),
                 " ".join("%s:%d" % (t[:5], item["per_tool"].get(t, 0)) for t in tools)))

    out = A.out or os.path.join(io.results_dir(), "view_frame_%s.html"
                                % (wanted[0] if len(wanted) == 1
                                   else "%d-%d" % (wanted[0], wanted[-1])))
    write_html(data, out, tools, A.rotate, A.colour_by, scheme)
    print("\n  wrote %s  (%.1f MB)" % (out, os.path.getsize(out) / 1e6))
    if A.open_browser:
        webbrowser.open("file://" + os.path.abspath(out))


if __name__ == "__main__":
    main()
