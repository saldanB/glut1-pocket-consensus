"""Build the MDpocket inputs for the two systems of the paper (R3 #9).

GLUT1: the 500 prepared frames, the same ones given to Fpocket.
AR:    59 frames, reconstructed from the ATOM records of the Fpocket outputs
       (each *_out.pdb contains the input protein alongside the alpha spheres).

Writes, per system: topology.pdb (first frame) and traj.pdb (multi-model), plus
traj.xtc if mdtraj or MDAnalysis is available.

    python3 prepare_mdpocket.py glut1 <frames_dir> <out_dir> [--pattern output.pdb.{}.pdb] [--n 500]
    python3 prepare_mdpocket.py ar    <fpocket_AR_dir> <out_dir>
"""
import os, sys, glob, re, argparse

import numpy as np

ap = argparse.ArgumentParser()
ap.add_argument("system", choices=["glut1", "ar"])
ap.add_argument("src")
ap.add_argument("dst")
ap.add_argument("--pattern", default="output.pdb.{}.pdb")
ap.add_argument("--n", type=int, default=500)
ap.add_argument("--no-align", dest="align", action="store_false",
                help="do not superpose the frames (not recommended: MDpocket "
                     "accumulates on a fixed grid)")
A = ap.parse_args()
os.makedirs(A.dst, exist_ok=True)


def kabsch(P, Q):
    """Rotation taking P onto Q; both are assumed already centred."""
    V, S, W = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(V @ W))
    return V @ np.diag([1.0, 1.0, d]) @ W


def coordinates(lines):
    return np.array([[float(l[30:38]), float(l[38:46]), float(l[46:54])] for l in lines])


def rewrite(lines, xyz):
    return ["%s%8.3f%8.3f%8.3f%s" % (l[:30], c[0], c[1], c[2], l[54:])
            for l, c in zip(lines, xyz)]


def ca_indices(lines):
    return [i for i, l in enumerate(lines)
            if l.startswith("ATOM") and l[12:16].strip() == "CA"]


def atoms(path, atom_only=False):
    out = []
    for l in open(path, errors="replace"):
        if l.startswith("ATOM"):
            out.append(l)
        elif l.startswith("HETATM") and not atom_only:
            out.append(l)
    return out


if A.system == "glut1":
    files = [os.path.join(A.src, A.pattern.format(i)) for i in range(1, A.n + 1)]
    files = [f for f in files if os.path.exists(f)]
    atom_only = False
else:
    rx = re.compile(r"AR_frame\.pdb\.(\d+)_out\.pdb$")
    found = []
    for p in glob.glob(os.path.join(A.src, "*_out.pdb")):
        m = rx.search(os.path.basename(p))
        if m:
            found.append((int(m.group(1)), p))
    files = [p for _, p in sorted(found)]
    atom_only = True          # discard Fpocket's STP alpha spheres

if not files:
    sys.exit("no frame found in " + A.src)

counts, frames = set(), []
for p in files:
    a = atoms(p, atom_only)
    counts.add(len(a))
    frames.append(a)

print("%s: %d frames, atoms per frame: %s" % (A.system, len(frames), sorted(counts)))
if len(counts) != 1:
    sys.exit("ERROR: the atom count is not constant across frames; MDpocket "
             "cannot read such a trajectory")

# Superposition. MDpocket accumulates onto a grid fixed in the laboratory frame,
# so an unsuperposed trajectory blurs the frequency map. The same alignment as
# the rest of the pipeline is used: rigid body on the protein CA atoms, with the
# first frame as reference.
if A.align:
    ica = ca_indices(frames[0])
    if not ica:
        sys.exit("no CA atom found: cannot superpose")
    reference = coordinates(frames[0])
    ref_ca = reference[ica]
    ref_centre = ref_ca.mean(0)
    rmsd = []
    for k in range(len(frames)):
        X = coordinates(frames[k])
        ca = X[ica]
        centre = ca.mean(0)
        R = kabsch(ca - centre, ref_ca - ref_centre)
        Y = (X - centre) @ R + ref_centre
        rmsd.append(float(np.sqrt(((Y[ica] - ref_ca) ** 2).sum(1).mean())))
        frames[k] = rewrite(frames[k], Y)
    print("  superposed %d frames on the %d protein CA atoms (reference: frame 1)"
          % (len(frames), len(ica)))
    print("  CA RMSD from the reference after superposition: median %.2f A, max %.2f A"
          % (float(np.median(rmsd)), float(max(rmsd))))
else:
    print("  WARNING: frames NOT superposed, the MDpocket map will be blurred")

topology = os.path.join(A.dst, "topology.pdb")
with open(topology, "w") as fh:
    fh.writelines(frames[0])
    fh.write("END\n")
print("  topology:", topology)

multi = os.path.join(A.dst, "traj.pdb")
with open(multi, "w") as fh:
    for k, a in enumerate(frames, 1):
        fh.write("MODEL %8d\n" % k)
        fh.writelines(a)
        fh.write("ENDMDL\n")
    fh.write("END\n")
print("  multi-model trajectory: %s (%.0f MB)" % (multi, os.path.getsize(multi) / 1e6))

try:
    import mdtraj as md
    t = md.load(multi)
    x = os.path.join(A.dst, "traj.xtc")
    t.save_xtc(x)
    print("  xtc trajectory: %s (%d frames, %d atoms) [mdtraj]" % (x, t.n_frames, t.n_atoms))
except Exception as e1:
    try:
        import MDAnalysis as mda
        u = mda.Universe(multi)
        x = os.path.join(A.dst, "traj.xtc")
        with mda.Writer(x, u.atoms.n_atoms) as w:
            for _ in u.trajectory:
                w.write(u.atoms)
        print("  xtc trajectory: %s [MDAnalysis]" % x)
    except Exception as e2:
        print("  xtc not written (neither mdtraj nor MDAnalysis available):")
        print("    mdtraj:      %s" % e1)
        print("    MDAnalysis:  %s" % e2)
        print("    the multi-model PDB can be used instead, if mdpocket accepts")
        print("    --trajectory_format pdb")
