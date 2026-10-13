"""IoU trajectory clusters and their correspondence with the HDBSCAN consensus regions.

1. Every pocket is rasterised on a 1 A grid; a voxel is occupied when within R = 5 A
   of a pocket-defining element (same representation as the pairwise IoU analysis).
2. In each frame, two pockets from different tools are linked when their IoU exceeds
   the threshold (0.5); connected components of linked pockets are frame-level
   volumetric consensus sites (tool set, centroid = mean of the superposed centroids
   of the member pockets).
3. Frame-level sites are linked across frames by agglomerative clustering of their
   centroids (average linkage, cut at 6.5 A): each cluster is an IoU trajectory.
4. IoU occupancy = fraction of the phase frames in which the trajectory has a site.
   A tool supports a trajectory when it takes part in it in at least 20% of the
   frames in which the trajectory is present. Trajectories with >= 4 tools and
   occupancy >= 20% are retained.
5. Each retained trajectory is matched to the nearest HDBSCAN consensus region
   (centre distance <= 6.5 A).
Usage: python run/07_iou_trajectories.py [iou_threshold] [radius]   (default 0.5 5.0; the paper also reports 0.3 and 0.7)
Writes results/07_iou_trajectories_<tag>.csv, 07_iou_vs_hdbscan_<tag>.csv,
       07_iou_occupancy_by_simultaneous_tools_<tag>.csv and, for 0.5 / 5 A, 07_iou_vs_hdbscan.csv
"""
import sys, os, itertools
import numpy as np, pandas as pd
from multiprocessing import Pool
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.sparse.csgraph import connected_components
from scipy.sparse import coo_matrix

THR = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
R = float(sys.argv[2]) if len(sys.argv) > 2 else 5.0
VOX = 1.0
TOOLS = ['AutoSite', 'DeepPocket', 'Fpocket', 'P2Rank', 'SiteFerret', 'SiteHound']
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
r = int(np.ceil(R / VOX))
g = np.mgrid[-r:r + 1, -r:r + 1, -r:r + 1].reshape(3, -1).T
OFF = g[(np.linalg.norm(g * VOX, axis=1) <= R)]

def keys(pts):
    v = np.unique(np.floor(pts / VOX).astype(np.int64), axis=0)
    d = (v[:, None, :] + OFF[None]).reshape(-1, 3) + 1000
    return np.unique(d[:, 0] * 4_000_000 + d[:, 1] * 2000 + d[:, 2])

CL = {}
def load():
    for t in TOOLS:
        d = np.load(f'{ROOT}/data/pointclouds/pc_GLUT1_{t}.npz', allow_pickle=True)
        off = np.concatenate([[0], np.cumsum(d['counts'])])
        CL[t] = (d['points'], off, d['frames'], d['pids'])

def frame_sites(f):
    pk = []
    for t in TOOLS:
        P, off, fr, pid = CL[t]
        for i in np.where(fr == f)[0]:
            pts = P[off[i]:off[i + 1]]
            if len(pts): pk.append((t, int(pid[i]), pts.mean(0), keys(pts)))
    n = len(pk); rows, cols = [], []
    for a, b in itertools.combinations(range(n), 2):
        if pk[a][0] == pk[b][0] or np.linalg.norm(pk[a][2] - pk[b][2]) > 2 * R + 15: continue
        inter = np.intersect1d(pk[a][3], pk[b][3], assume_unique=True).size
        if inter == 0: continue
        iou = inter / (pk[a][3].size + pk[b][3].size - inter)
        if iou > THR: rows.append(a); cols.append(b)
    out = []
    if not rows: return out
    m = coo_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
    nc, lab = connected_components(m, directed=False)
    for c in range(nc):
        mem = np.where(lab == c)[0]
        if len(mem) < 2: continue
        out.append(dict(frame=f, members=[(pk[i][0], pk[i][1]) for i in mem]))
    return out

if __name__ == '__main__':
    print(f'\n07  IoU trajectories, threshold {THR}, R {R} A\n' + '-' * 50)
    load()
    frames = list(range(1, 501))
    with Pool(2, initializer=load) as pool:
        res = pool.map(frame_sites, frames, chunksize=5)
    al = pd.read_csv(f'{ROOT}/results/pockets_aligned.csv').set_index(['tool', 'frame', 'pocket_id'])
    sites = []
    for lst in res:
        for s in lst:
            xyz = np.array([al.loc[(t, s['frame'], p), ['x', 'y', 'z']].values.astype(float) for t, p in s['members']])
            sites.append(dict(frame=s['frame'], tools='|'.join(sorted({t for t, _ in s['members']})), x=xyz[:, 0].mean(), y=xyz[:, 1].mean(), z=xyz[:, 2].mean()))
    S = pd.DataFrame(sites)
    S['phase'] = np.where(S.frame <= 100, 'Phase1', 'Phase2')
    reg = pd.read_csv(f'{ROOT}/results/02_consensus_regions.csv'); reg = reg[reg.scenario == 'with superposition']
    tag = f't{THR}_R{R}'
    out = []
    for ph, nfr in (('Phase1', 100), ('Phase2', 400)):
        s = S[S.phase == ph].reset_index(drop=True)
        lab = fcluster(linkage(s[['x', 'y', 'z']].values, 'average'), 6.5, 'distance')
        s['traj'] = lab
        rg = reg[reg.phase == ph]
        for k, gg in s.groupby('traj'):
            fr = gg.frame.unique(); occ = len(fr) / nfr
            cnt = {t: gg[gg.tools.str.contains(t)].frame.nunique() for t in TOOLS}
            sup = [t for t in TOOLS if cnt[t] >= 0.2 * len(fr)]
            c = gg[['x', 'y', 'z']].mean().values
            d = np.linalg.norm(rg[['x', 'y', 'z']].values - c, axis=1); j = int(d.argmin())
            out.append(dict(phase=ph, traj=k, occupancy=round(occ, 3), n_tools=len(sup), tools='|'.join(sup),
                            x=c[0], y=c[1], z=c[2], n_sites=len(gg),
                            region=int(rg.iloc[j].region) if d[j] <= 6.5 else None, distance=round(float(d[j]), 2),
                            region_tools=rg.iloc[j].tools if d[j] <= 6.5 else None))
    T = pd.DataFrame(out)
    T['retained'] = (T.n_tools >= 4) & (T.occupancy >= 0.2)
    T.sort_values(['phase', 'occupancy'], ascending=[True, False]).to_csv(f'{ROOT}/results/07_iou_trajectories_{tag}.csv', index=False)
    K = T[T.retained]
    print(K.sort_values(['phase', 'occupancy'], ascending=[True, False]).to_string(index=False))

    # correspondence with the HDBSCAN regions and occupancy by number of simultaneous tools
    rows, rows2 = [], []
    for ph, nf in (('Phase1', 100), ('Phase2', 400)):
        s = S[S.phase == ph].reset_index(drop=True)
        s['traj'] = fcluster(linkage(s[['x', 'y', 'z']].values, 'average'), 6.5, 'distance')
        s['nt'] = s.tools.str.count(r'\|') + 1
        for _, r in reg[reg.phase == ph].iterrows():
            mm = K[(K.phase == ph) & (K.region == r.region)].sort_values('distance')
            b = None if mm.empty else mm.iloc[0]
            rows.append(dict(phase=ph, region=r.region, ntools=r.n_tools, soft=round(100 * r.soft_persistence, 1),
                             strict=round(100 * r.strict_persistence, 1),
                             iou_occ=None if b is None else b.occupancy, dist=None if b is None else b.distance,
                             iou_ntools=None if b is None else b.n_tools, iou_tools=None if b is None else b.tools))
            row = dict(region=f"P{ph[-1]}-R{r.region}", soft_persistence=round(100 * r.soft_persistence, 1))
            for k in (2, 3, 4, 5, 6):
                row[f'iou_occ_min{k}tools'] = None if b is None else round(100 * s[(s.traj == b.traj) & (s.nt >= k)].frame.nunique() / nf)
            rows2.append(row)
    pd.DataFrame(rows).to_csv(f'{ROOT}/results/07_iou_vs_hdbscan_{tag}.csv', index=False)
    pd.DataFrame(rows2).to_csv(f'{ROOT}/results/07_iou_occupancy_by_simultaneous_tools_{tag}.csv', index=False)
    if THR == 0.5 and R == 5.0:   # the manuscript setting, read by run/11
        pd.DataFrame(rows).to_csv(f'{ROOT}/results/07_iou_vs_hdbscan.csv', index=False)
