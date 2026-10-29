import sys, itertools, numpy as np, MDAnalysis as mda
from MDAnalysis.lib.distances import capped_distance
top, trj = sys.argv[1], sys.argv[2]
stride = int(sys.argv[3]) if len(sys.argv) > 3 else 1
u = mda.Universe(top, trj)
atoms = u.select_atoms('not name H*')
shifts = [s for s in itertools.product((-1, 0, 1), repeat=3) if s != (0, 0, 0)]
overall = 1e9
for ts in u.trajectory[::stride]:
    P = atoms.positions
    B = ts.triclinic_dimensions
    best = 1e9
    for s in shifts:
        _, d = capped_distance(P, P + np.dot(s, B), 40.0, return_distances=True)
        if len(d):
            best = min(best, d.min())
    overall = min(overall, best)
    print(f'frame {ts.frame:5d}  min distance to image {best:6.2f} A')
print(f'MINIMUM over trajectory: {overall:.2f} A')
