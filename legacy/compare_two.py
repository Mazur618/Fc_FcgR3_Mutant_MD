#!/usr/bin/env python3
"""Compare MM-GBSA per-frame energies of two runs (e.g. WT vs a mutant).
usage: python compare_two.py WT_perframe_s301.csv D265A/rep1/perframe.csv [--nblocks 5] [--labels WT D265A]
Both runs must use the same startframe/interval. Prints mean per term, ddG = B - A,
and an uncertainty from block means. Saves compare_two.png."""
import argparse
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('a'); ap.add_argument('b')
ap.add_argument('--nblocks', type=int, default=5)
ap.add_argument('--labels', nargs=2, default=['A', 'B'])
args = ap.parse_args()

def load(path):
    L = open(path).read().splitlines()
    s = next(i for i, l in enumerate(L) if 'Energy Terms' in l and l.strip().lower().startswith('delta'))
    hdr = [h.strip() for h in L[s + 1].split(',')]
    rows = []
    for l in L[s + 2:]:
        try: rows.append([float(x) for x in l.split(',')])
        except ValueError: break
    return hdr, np.array(rows)

def mean_se(x, nb):
    m = np.mean(x)
    bm = np.array([b.mean() for b in np.array_split(x, nb)])
    return m, bm.std(ddof=1) / np.sqrt(nb)

ha, A = load(args.a); hb, B = load(args.b)
print(f'{args.labels[0]}: {len(A)} frames   {args.labels[1]}: {len(B)} frames')
if len(A) != len(B):
    print('WARNING: different frame counts, so the windows probably differ. Rerun with the same startframe/interval.')
terms = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'DELTA G gas', 'DELTA G solv', 'DELTA TOTAL']
print(f'\n{"term":14s} {args.labels[0]:>10s} {args.labels[1]:>10s} {"diff(B-A)":>10s} {"+-":>6s}')
for t in terms:
    ma, sa = mean_se(A[:, ha.index(t)], args.nblocks)
    mb, sb = mean_se(B[:, hb.index(t)], args.nblocks)
    print(f'{t:14s} {ma:10.2f} {mb:10.2f} {mb - ma:10.2f} {np.hypot(sa, sb):6.2f}')
print('\n+- is the block-mean SEM of each run combined in quadrature. It reflects within-run '
      'noise only; with one replicate each it understates the real uncertainty.')

ta, tb = A[:, ha.index('DELTA TOTAL')], B[:, hb.index('DELTA TOTAL')]
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(ta, lw=0.7, label=args.labels[0]); ax[0].plot(tb, lw=0.7, label=args.labels[1])
ax[0].set_xlabel('MM-GBSA frame index'); ax[0].set_ylabel('DELTA TOTAL (kcal/mol)'); ax[0].legend()
bins = np.linspace(min(ta.min(), tb.min()), max(ta.max(), tb.max()), 25)
ax[1].hist(ta, bins, alpha=0.5, label=args.labels[0]); ax[1].hist(tb, bins, alpha=0.5, label=args.labels[1])
ax[1].set_xlabel('DELTA TOTAL (kcal/mol)'); ax[1].legend()
plt.tight_layout(); plt.savefig('compare_two.png', dpi=150)
print('saved compare_two.png')
