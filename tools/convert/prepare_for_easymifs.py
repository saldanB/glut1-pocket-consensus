"""Prepare one frame for the EasyMIFs / SiteHound pipeline.

Two edits, both required by pdb2gmx when the structure is typed with the GROMOS
force field:

  * remove D-glucose. In these frames the substrate is written as ATOM records
    with residue name GLC, so prepare_pdb.py -k (which drops HETATM only) does
    not remove it, and the force field has no entry for it:

        Fatal error: Atom O1 in residue GLCA 448 not found in rtp entry

    The SiteHound run is therefore performed on the apo structure, unlike the
    other five detectors. This is stated in the Methods.

  * rename HID to HIS. Delta-protonated histidines are named HID and the force
    field does not define their hydrogens:

        Fatal error: There were 8 missing atoms in molecule Protein

    pdb2gmx -missing would also start, but it would silently alter the charges.

    python3 prepare_for_easymifs.py <input.pdb> <output.pdb>
"""
import sys

src, dst = sys.argv[1], sys.argv[2]
removed = renamed = 0
with open(src) as fh, open(dst, "w") as out:
    for line in fh:
        if line.startswith(("ATOM", "HETATM")):
            resn = line[17:20]
            if resn == "GLC":
                removed += 1
                continue
            if resn == "HID":
                line = line[:17] + "HIS" + line[20:]
                renamed += 1
        out.write(line)

print("[prepare] removed %d glucose atoms, renamed %d HID atoms to HIS"
      % (removed, renamed))
