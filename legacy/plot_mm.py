#!/usr/bin/env python3
"""Minimal MM-GBSA plots for one replicate.
usage: python plot_mm.py [rep_folder] [--start 301] [--step 5] [--dt 0.01] [--csv perframe.csv]
Writes mm.png in the folder and prints mean / SD / block means of DELTA TOTAL."""
import argparse, os
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('folder', nargs='?', default='.')
ap.add_argument('--start', type=int, default=301, help='MM-GBSA startframe (1-based)')
ap.add_argument('--step', type=int, default=5, help='MM-GBSA interval')
ap.add_argument('--dt', type=float, default=0.01, help='ns per trajectory frame')
ap.add_argument('--csv', default='perframe.csv')
a = ap.parse_args()

L = open(os.path.join(a.folder, a.csv)).read().splitlines()
s = next(i for i, l in enumerate(L) if 'Energy Terms' in l and l.strip().lower().startswith('delta'))
hdr = [h.strip() for h in L[s + 1].split(',')]
rows = []
for l in L[s + 2:]:
    try:
        rows.append([float(x) for x in l.split(',')])
    except ValueError:
        break
d = np.array(rows)
col = {h: j for j, h in enumerate(hdr)}
t = (a.start - 1 + a.step * np.arange(len(d))) * a.dt
tot = d[:, col['DELTA TOTAL']]

w = max(5, len(tot) // 10)
run = np.convolve(np.pad(tot, (w // 2, w - 1 - w // 2), mode='edge'), np.ones(w) / w, mode='valid')

fig, ax = plt.subplots(1, 2, figsize=(11, 4))
ax[0].plot(t, tot, lw=0.6, alpha=0.4)
ax[0].plot(t, run, lw=1.6)
ax[0].set_xlabel('time (ns)'); ax[0].set_ylabel('kcal/mol'); ax[0].set_title('MM-GBSA DELTA TOTAL')
for k in ('VDWAALS', 'EEL', 'EGB', 'ESURF'):
    ax[1].plot(t, d[:, col[k]], lw=0.8, label=k)
ax[1].plot(t, d[:, col['EEL']] + d[:, col['EGB']], lw=1.2, color='k', label='EEL+EGB')
ax[1].set_xlabel('time (ns)'); ax[1].set_title('Components'); ax[1].legend(fontsize=7, ncol=2)
plt.tight_layout()
plt.savefig(os.path.join(a.folder, 'mm.png'), dpi=150)

blocks = [b.mean() for b in np.array_split(tot, 4)]
print(f'{len(tot)} frames  mean {tot.mean():.2f}  sd {tot.std(ddof=1):.2f}')
print('4-block means: ' + ' '.join(f'{b:.2f}' for b in blocks))
