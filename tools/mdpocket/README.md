# MDpocket cross-check

MDpocket is used in the revised manuscript as an external reference point: it
detects cavities on a trajectory directly, so it answers the same question as
this pipeline by a different route. These two scripts reproduce that comparison.

They are not part of the reproduction path for the main results. They need the
trajectory frames, not the deposited pocket table, and MDpocket itself (shipped
with Fpocket).

```bash
# 1. build the inputs: superposed frames, topology + trajectory
python3 prepare_mdpocket.py glut1 <prepared_frames_dir> mdp_glut1/

# 2. de novo pass: the frequency map over the ensemble
cd mdp_glut1 && mdpocket --trajectory_file traj.xtc \
                         --trajectory_format xtc \
                         -f topology.pdb -o mdpout

# 3. compare the map with the consensus regions
python3 compare_mdpocket.py glut1 mdp_glut1/ mdp_glut1/topology.pdb results/ \
                            --table4 <Table4.csv>
```

For AR, `prepare_mdpocket.py ar <fpocket_AR_dir> mdp_ar/` rebuilds the 59 frames
from the ATOM records of the Fpocket outputs, and `compare_mdpocket.py ar` scores
the recovery of the four residues of the specificity pocket instead of the
Jaccard overlap.

Two points about step 1 that matter for the result. MDpocket accumulates
occupancy on a grid fixed in the laboratory frame, so the frames must be
superposed first or the frequency map is blurred by the protein's own rotation;
`prepare_mdpocket.py` applies the same rigid-body CA alignment used everywhere
else in the pipeline, and prints the resulting RMSD. And the atom count must be
identical in every frame, which the script checks and refuses to proceed
without.

The comparison in step 3 is made on lining residues at 5 A, the convention of
Table 4, rather than on centre-to-centre distances; the reason is given in the
header of `compare_mdpocket.py`.
