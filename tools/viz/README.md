# Frame viewer

Look at one frame of the ensemble with the predicted pockets projected onto it.
This is an inspection tool, not part of the reproduction path: nothing in `run/`
depends on it, and it writes nothing into `data/` except its frame cache.

```bash
python tools/viz/view_frame.py 42                      # one frame, all six detectors
python tools/viz/view_frame.py 42 --tools Fpocket P2Rank --top-k 5
python tools/viz/view_frame.py 1-24 --out ensemble.html --open
python tools/viz/view_frame.py 42 --hulls
python tools/viz/view_frame.py 42 --clouds --consensus --published
```

## Which frames

The frame argument takes a single frame or an interval, and more than one of
either:

| argument | frames |
|---|---|
| `42` | just that one |
| `1-24` | the closed interval, 24 frames |
| `1,5,9` | only those |
| `1-500:25` | every 25th frame across the ensemble |
| `1-10,42,200-400:50` | any combination |

Duplicates collapse and a reversed interval such as `10-1` is read as `1-10`.
More than one frame gives a slider with prev/next/play; the camera is kept as you
step, so a feature stays where you put it. Dragging the slider is safe: a frame
whose load has been overtaken by a later one is discarded rather than left on the
stage, and playback waits for each frame instead of firing on a timer. Each frame adds about 0.6 MB to the
page, so a long interval is better sampled with a step than taken whole — and
`--clouds` multiplies that several times over.

## Orientation

The initial view is turned a quarter turn about x, which puts the GLUT1
transmembrane axis across the screen rather than pointing at the camera. Change it
with `--rotate X Y Z` in degrees, applied about x, then y, then z:

```bash
python tools/viz/view_frame.py 42 --rotate 0 0 0     # NGL's own default view
python tools/viz/view_frame.py 42 --rotate 90 0 180  # and turned about z as well
```

Whatever `--rotate` produces becomes the home view, so **Reset view** returns to
it, orientation included.

The default output is `results/view_frame_<frames>.html`, a standalone page that
loads NGL from a CDN; `--open` shows it in a browser. In a notebook you get a
widget instead:

```python
from tools.viz.view_frame import view
view(42, consensus=True)
```

## The frame cache

The frames are deposited as one 60 MB archive
(doi:[10.5281/zenodo.21238931](https://doi.org/10.5281/zenodo.21238931)), which
would be a poor trade for a viewer that looks at a few of them. Zenodo serves
byte ranges, so `frames.py` reads the archive's central directory once, caches it
in `data/frames/.archive_index.json`, and then pulls individual members — about
0.6 MB each. A frame already in `data/frames/` is never fetched again. If the
server stops honouring range requests the whole archive is downloaded once
instead and read locally, so the tool keeps working either way.

`data/frames/` is in `.gitignore`.

## What is drawn

| component | source | representation |
|---|---|---|
| the protein | the deposited frame | cartoon, optional double-sided molecular surface |
| predicted pockets | `data/pockets_master.csv`, one pseudo-atom per pocket | spacefill, one colour per detector |
| convex hull per cluster (`--hulls`) | the same points, hulled per cluster | translucent double-sided mesh, one colour per cluster |
| pocket-defining points (`--clouds`) | `data/pointclouds/`, for this frame only | points, coloured by HDBSCAN cluster or by detector |
| consensus regions (`--consensus`) | `results/02_consensus_regions.csv` | spacefill, gold |
| published centres (`--published`) | `data/published_consensus_centres.csv` | spacefill, pink |

Each detector gets its own chain (`A` to `F`, consensus `Z`) and residue name, so
the NGL selections are plain chain selections and the sidebar can toggle a
detector without rebuilding anything. Pocket rank is carried in the occupancy
column and the number of pocket-defining points in the B-factor column, so both
are available in NGL's tooltip.

## One surface per cluster

The raw points are hard to read once there are ten thousand of them, so `--hulls`
replaces each cluster's points with the convex hull of them:

```bash
python tools/viz/view_frame.py 42 --hulls
python tools/viz/view_frame.py 42 --hulls --clouds   # hull and its points together
```

For the displayed frame, the points of every pocket belonging to one cluster are
pooled and their convex hull taken, so there is one closed surface per cluster per
frame and it changes as the frame changes — on frames 40 to 43 the hulls go
51, 50, 54, 51 and their total volume 60,163, 65,853, 74,572, 65,952 A^3. Each hull
takes the colour its cluster's points would have, and its tooltip gives the
detector, the cluster, the volume and the number of points enclosed.

The meshes are built with `scipy.spatial.ConvexHull` and handed to NGL as a
`Shape`, one per detector so the detector checkboxes still apply. They are flat
shaded: a triangle carries its own three vertices and its own face normal, wound
outwards, so the facets stay visible rather than being rounded off. A frame's hulls
are about 380 KB, less than half what the points cost.

A cluster needs four points in general position to have a hull, so two cases are
skipped and reported rather than drawn: a cluster whose points in this frame are
coplanar, and **SiteHound**, which contributes one point per pocket and therefore
never has enough. Pockets HDBSCAN called noise have no cluster and so no hull.

## Colouring the points by cluster

`--clouds` draws the pocket-defining points of the pockets **of the displayed
frame** — alpha-sphere centres, energy grid points or atomic coordinates,
depending on the detector. By default each point takes the colour of the HDBSCAN
cluster its pocket belongs to, so a site the detector keeps returning to across
the ensemble reads as one colour wherever you are in the trajectory, and a pocket
that HDBSCAN called noise is grey.

```bash
python tools/viz/view_frame.py 42 --clouds                   # by cluster
python tools/viz/view_frame.py 42 --clouds --colour-by tool  # by detector
```

The cluster of a pocket is looked up, not recomputed: `src.consensus.cluster_labels`
runs once per phase on the superposed centroids, the same call `run/02` makes with
the same `min_cluster_size`, and each pocket is matched to its cluster by
`(tool, frame, pocket_id)`. The point archives carry the same `(frame, pocket_id)`
keys as `pockets_master.csv`, so the match is exact. Cluster numbers restart per
detector, so the colours are keyed on the pair: cluster 7 of AutoSite and cluster 7
of Fpocket are different clusters with different colours.

| | Phase 1 | Phase 2 |
|---|---|---|
| clusters over all six detectors | 81 | 95 |
| points in a typical frame | ~10,500 | ~10,000 |
| added per frame | ~0.8 MB | ~0.8 MB |

Colouring by cluster needs the superposed centroids, so it implies `--align` and
therefore `results/pockets_aligned.csv` from `run/01_align.py`; `--colour-by tool`
does not. The colours come from one NGL selection scheme per phase, registered once
rather than per frame, so a cluster keeps its colour as you step the slider.

The detector checkboxes filter both per-detector layers at once — pockets and
points — and the Layers group turns each layer on and off. Every layer is listed
whether or not it was built; one that was not is greyed out and shows the flag
that builds it, so the page says what else it could have shown.

## Superposition

The deposited frames and the centroids in `pockets_master.csv` are both in the
unsuperposed space of the frame as supplied to the detectors, so by default they
are drawn as they are and overlay exactly — the frame's Calpha coordinates agree
with `data/protein_ca.npz` to about 2e-6 A.

The consensus regions of `run/02` are in the superposed space of the reference
frame, so `--consensus` first superposes the frame and its pockets, using the
same Calpha Kabsch fit as `src/align.py`. `--align` does that without drawing the
regions. The sidebar always states which space is on screen, because for GLUT1
the superposition moves a centroid by of the order of 1 A and it would otherwise
be impossible to tell.

`--consensus` needs `results/02_consensus_regions.csv`, so run `run/02_consensus.py`
first; the tool says so if it is missing. `--published` draws the centres reported
in the paper instead, which ship in `data/`, so it works in a fresh clone. Both can
be on at once, and the sidebar toggles either without reloading the frame; their
counts follow the phase, since the two phases have different reference sets.

## The files here

| file | what it is |
|---|---|
| `view_frame.py` | the command line tool and the notebook widget |
| `page.html` | the standalone page `--out` writes, with `__PLACEHOLDERS__` |
| `overlay.py` | the PDB and mesh strings handed to NGL |
| `frames.py` | the frame cache, and the ranged fetch that fills it |

`page.html` is the real page, not a fragment: `write_html` reads it and substitutes
the eleven `__PLACEHOLDER__` tokens, and nothing else. It is therefore valid HTML on
its own, so the markup, the stylesheet and the viewer's JavaScript can be edited and
checked in place with the editor's own tooling. Opening it directly shows the chrome
with no frame loaded.

## Requirements

The HTML output needs only what the analysis environment already has. The
notebook widget additionally needs `nglview`, which is not in `environment.yml`
because it is not needed to reproduce anything:

```bash
conda install -c conda-forge nglview
```
