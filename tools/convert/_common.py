"""Shared helpers for the format converters.

Each converter reads one detector's native output for a frame and writes a
combined PDB (protein + pocket-defining points) using the same file naming and
the same record layout as the deposited files, so that everything downstream
reads them unchanged.
"""
import os


def protein_lines(prepared_dir, frame):
    """ATOM/TER records of the input frame given to the detector."""
    path = os.path.join(prepared_dir, "output.pdb.%d.pdb" % frame)
    with open(path, errors="replace") as fh:
        return [l.rstrip("\n") for l in fh if l.startswith(("ATOM", "TER"))]


def _record(rec, serial, resn, resi, x, y, z):
    return ("%-6s%5d  C   %3s  %4d    %8.3f%8.3f%8.3f  1.00  0.00           C"
            % (rec, serial % 100000, resn[:3], resi % 10000, x, y, z))


def write_pdb(dest, protein, groups, rec="HETATM", resn="STP", remarks=None):
    """groups: list of (pocket_id, [(x, y, z), ...]) in rank order."""
    lines = list(remarks or []) + list(protein)
    serial = len(protein) + 1
    for pid, pts in groups:
        for (x, y, z) in pts:
            lines.append(_record(rec, serial, resn, int(pid), x, y, z))
            serial += 1
    lines.append("END")
    with open(dest, "w") as fh:
        fh.write("\n".join(lines) + "\n")


def points_from_pqr(path):
    """Coordinates from a .pqr file.

    PQR is whitespace-separated rather than column-formatted, and the trailing
    charge and radius columns shift the coordinates, so x, y, z are the fifth,
    fourth and third fields from the end.
    """
    out = []
    with open(path, errors="replace") as fh:
        for l in fh:
            if not l.startswith(("ATOM", "HETATM")):
                continue
            c = l.split()
            try:
                out.append((float(c[-5]), float(c[-4]), float(c[-3])))
            except (ValueError, IndexError):
                pass
    return out


def points_from_pdb(path):
    out = []
    with open(path, errors="replace") as fh:
        for l in fh:
            if l.startswith(("ATOM", "HETATM")):
                try:
                    out.append((float(l[30:38]), float(l[38:46]), float(l[46:54])))
                except ValueError:
                    pass
    return out


def summarise(dest, resn, expected=""):
    """Point count per frame, as a check against the deposited files."""
    import glob
    files = sorted(glob.glob(os.path.join(dest, "*.pdb")))
    if not files:
        print("   no file produced in %s" % dest)
        return
    total = 0
    for p in files[:20]:
        with open(p, errors="replace") as fh:
            total += sum(1 for l in fh
                         if l.startswith(("ATOM", "HETATM")) and l[17:20].strip() == resn)
    print("   %d files; mean %s records per frame over the first 20: %.0f   %s"
          % (len(files), resn, total / min(20, len(files)), expected))


def write_centroids(path, rows, tool, header, fmt, phase1_max=100):
    """Centroid table in the layout of data/pockets_master.csv."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(header + "\n")
        for row in rows:
            fh.write(fmt(row, "Phase1" if row[0] <= phase1_max else "Phase2") + "\n")
    print("   centroids: %d rows -> %s" % (len(rows), path))
