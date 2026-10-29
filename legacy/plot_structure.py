#!/usr/bin/env python3
"""Heatmap and per-metric plots from structure_results/replicates.csv.

usage (prep folder, after collect_structure.py):
    python plot_structure.py

Writes structure_results/heatmap_vs_WT.png and strip_by_variant.png.
Heatmap colour = (variant mean - WT mean) / pooled within-system SD of that
metric, so values near 0 are within normal replicate scatter. Cell text is the
raw difference from WT.
"""
import argparse
import csv
import math

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--indir', default='structure_results')
a = ap.parse_args()

METRICS = ['CH3_fit', 'Fc_CH2', 'CD16a_rmsd', 'gly_Fc', 'gly_stub',
           'CD16a_self', 'Fc_self', 'COM', 'Rg', 'BSA', 'dG']
COLORS = {'wt': 'black', 'reduced': 'tab:red', 'wtlike': 'tab:gray',
          'improved': 'tab:green'}

rows = list(csv.DictReader(open(a.indir + '/replicates.csv')))
systems = []
cls = {}
for r in rows:
    if r['system'] not in systems:
        systems.append(r['system'])
        cls[r['system']] = r['class']
order = ['WT'] + sorted([s for s in systems if s != 'WT'],
                        key=lambda s: ({'reduced': 0, 'wtlike': 1,
                                        'improved': 2}.get(cls[s], 3), s))


def vals(system, m):
    return np.array([float(r[m]) for r in rows if r['system'] == system])


mean = {s: {m: vals(s, m).mean() for m in METRICS} for s in order}
pooled = {}
for m in METRICS:
    var = [vals(s, m).var(ddof=1) for s in order if len(vals(s, m)) > 1]
    pooled[m] = math.sqrt(np.mean(var)) if var else float('nan')

others = [s for s in order if s != 'WT']
z = np.array([[(mean[s][m] - mean['WT'][m]) / pooled[m] for m in METRICS]
              for s in others])
diff = np.array([[mean[s][m] - mean['WT'][m] for m in METRICS]
                 for s in others])

fig, ax = plt.subplots(figsize=(1.0 * len(METRICS) + 3, 0.6 * len(others) + 2.5))
lim = max(2.0, np.nanmax(np.abs(z)))
im = ax.imshow(z, cmap='RdBu_r', vmin=-lim, vmax=lim, aspect='auto')
for i in range(len(others)):
    for j in range(len(METRICS)):
        tc = 'white' if abs(z[i, j]) > 0.6 * lim else 'black'
        ax.text(j, i, '{:+.1f}'.format(diff[i, j]), ha='center', va='center',
                fontsize=7, color=tc)
ax.set_xticks(range(len(METRICS)))
ax.set_xticklabels(METRICS, rotation=45, ha='right')
ax.set_yticks(range(len(others)))
ax.set_yticklabels(['{} [{}]'.format(s, cls[s]) for s in others])
fig.colorbar(im, ax=ax, label='difference from WT in pooled-SD units')
ax.set_title('Variant minus WT (text = raw difference; colour = in SD units)',
             fontsize=9)
plt.tight_layout()
plt.savefig(a.indir + '/heatmap_vs_WT.png', dpi=150)
plt.close()

ncol = 4
nrow = math.ceil(len(METRICS) / ncol)
fig, axs = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.2 * nrow))
rng = np.random.default_rng(1)
for ax, m in zip(axs.ravel(), METRICS):
    for i, s in enumerate(order):
        v = vals(s, m)
        c = COLORS.get(cls[s], 'b')
        ax.scatter(np.full(len(v), i) + rng.uniform(-0.08, 0.08, len(v)), v,
                   color=c, s=22, zorder=3)
        ax.hlines(v.mean(), i - 0.3, i + 0.3, color=c, lw=1.5)
    ax.axhline(mean['WT'][m], color='k', lw=0.6, ls='--')
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=60, fontsize=7)
    ax.set_title(m, fontsize=9)
for ax in axs.ravel()[len(METRICS):]:
    ax.axis('off')
plt.tight_layout()
plt.savefig(a.indir + '/strip_by_variant.png', dpi=150)
plt.close()
print('wrote {}/heatmap_vs_WT.png and strip_by_variant.png'.format(a.indir))
