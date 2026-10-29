#!/usr/bin/env python3
"""Collect structural metrics from every NAME/repN folder and compare with WT.

usage (prep folder, after queue_analysis.sh finished):
    python collect_structure.py [--from-frame 301]

Reads per replicate: rmsd.dat rmsd_self.dat comdist.dat rg.dat sasa.dat
perframe.csv summary.txt. Metrics are averaged over trajectory frames
>= --from-frame (default 301, the same start as the MM-GBSA window).

Writes structure_results/: replicates.csv variants.csv correlations.csv
scatter_dG_vs_structure.png and prints the tables.
"""
import argparse
import csv
import glob
import math
import os
import re

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('--variants', default='variants.txt')
ap.add_argument('--from-frame', type=int, default=301)
ap.add_argument('--outdir', default='structure_results')
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

METRICS = ['CH3_fit', 'Fc_CH2', 'CD16a_rmsd', 'gly_Fc', 'gly_stub',
           'CD16a_self', 'Fc_self', 'COM', 'Rg', 'BSA']
notes = []


def read_dat(path):
    with open(path) as fh:
        hdr = fh.readline().replace('#', '').split()
    return hdr, np.loadtxt(path)


def window_mean(path, col, start):
    hdr, d = read_dat(path)
    keep = d[:, 0] >= start
    return float(d[keep, hdr.index(col)].mean())


def bsa_mean(path, start):
    hdr, d = read_dat(path)
    f, cplx, fc, rc = d[:, 0], d[:, 1], d[:, 2], d[:, 3]
    bsa = (fc + rc - cplx) / 2
    good = (bsa > 0) & (bsa <= np.minimum(fc, rc))
    keep = (f >= start) & good
    dropped = int(((f >= start) & ~good).sum())
    return float(bsa[keep].mean()), dropped


def mm_dG(path):
    lines = open(path).read().splitlines()
    s = None
    for i, line in enumerate(lines):
        low = line.strip().lower()
        if 'energy terms' in low and low.startswith('delta'):
            s = i
            break
    hdr = [h.strip() for h in lines[s + 1].split(',')]
    rows = []
    for line in lines[s + 2:]:
        try:
            rows.append([float(x) for x in line.split(',')])
        except ValueError:
            break
    d = np.array(rows)
    return float(d[:, hdr.index('DELTA TOTAL')].mean())


def read_summary(path):
    out = {'min_image_A': float('nan'), 'offchair_frames': float('nan')}
    if not os.path.exists(path):
        return out
    off = 0
    found = False
    for line in open(path):
        m = re.search(r'MINIMUM over trajectory:\s+([\d.]+)', line)
        if m:
            out['min_image_A'] = float(m.group(1))
        m = re.match(r'^P\d+\s+\S+\s+\S+\s+\S+\s+(\d+)\s+\(', line)
        if m:
            off += int(m.group(1))
            found = True
    if found:
        out['offchair_frames'] = off
    return out


def variants_list():
    out = [('WT', 'wt')]
    for line in open(a.variants):
        s = line.strip()
        if not s:
            continue
        tok = s.lstrip('#').split()
        if len(tok) >= 3 and tok[2] in ('reduced', 'wtlike', 'improved'):
            if not s.startswith('#') or os.path.isdir(tok[0]):
                out.append((tok[0], tok[2]))
    return out


rows = []
for name, cls in variants_list():
    for d in sorted(glob.glob('{}/rep*/'.format(name))):
        rep = d.rstrip('/').split('rep')[-1]
        need = ['rmsd.dat', 'rmsd_self.dat', 'comdist.dat', 'rg.dat',
                'sasa.dat', 'perframe.csv']
        missing = [n for n in need if not os.path.exists(d + n)]
        if missing:
            notes.append('{} rep{}: missing {}'.format(name, rep, missing))
            continue
        st = a.from_frame
        bsa, dropped = bsa_mean(d + 'sasa.dat', st)
        row = {'system': name, 'class': cls, 'rep': rep,
               'CH3_fit': window_mean(d + 'rmsd.dat', 'CH3_fit', st),
               'Fc_CH2': window_mean(d + 'rmsd.dat', 'Fc_CH2', st),
               'CD16a_rmsd': window_mean(d + 'rmsd.dat', 'CD16a', st),
               'gly_Fc': window_mean(d + 'rmsd.dat', 'gly_Fc', st),
               'gly_stub': window_mean(d + 'rmsd.dat', 'gly_stub', st),
               'CD16a_self': window_mean(d + 'rmsd_self.dat', 'CD16a_self', st),
               'Fc_self': window_mean(d + 'rmsd_self.dat', 'Fc_self', st),
               'COM': window_mean(d + 'comdist.dat', 'COM_Fc_CD16a', st),
               'Rg': window_mean(d + 'rg.dat', 'Rg', st),
               'BSA': bsa, 'BSA_dropped_frames': dropped,
               'dG': mm_dG(d + 'perframe.csv')}
        row.update(read_summary(d + 'summary.txt'))
        rows.append(row)

if not rows:
    raise SystemExit('no complete replicate folders found')

cols = (['system', 'class', 'rep'] + METRICS + ['dG', 'BSA_dropped_frames',
        'min_image_A', 'offchair_frames'])
with open(a.outdir + '/replicates.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols)
    w.writeheader()
    w.writerows(rows)

# per-variant mean and SD between replicates, and difference from WT
systems = []
for r in rows:
    if r['system'] not in systems:
        systems.append(r['system'])
stats = {}
for s in systems:
    sub = [r for r in rows if r['system'] == s]
    stats[s] = {'n': len(sub), 'class': sub[0]['class']}
    for m in METRICS + ['dG']:
        v = np.array([r[m] for r in sub])
        stats[s][m] = (float(v.mean()),
                       float(v.std(ddof=1)) if len(v) > 1 else float('nan'))

with open(a.outdir + '/variants.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['system', 'class', 'n_reps'] +
               [m + suffix for m in METRICS + ['dG'] for suffix in ('_mean', '_sd')])
    for s in systems:
        line = [s, stats[s]['class'], stats[s]['n']]
        for m in METRICS + ['dG']:
            line += [round(stats[s][m][0], 3), round(stats[s][m][1], 3)]
        w.writerow(line)


def fmt(x, p=2):
    return 'nan' if (x is None or math.isnan(x)) else '{:.{}f}'.format(x, p)


print('\nPer-system mean (SD between replicates), frames >= {}:'.format(a.from_frame))
head = '{:8s} {:>3s}'.format('system', 'n') + ''.join(
    '{:>16s}'.format(m) for m in ['CD16a_rmsd', 'CD16a_self', 'COM', 'BSA', 'dG'])
print(head)
for s in systems:
    line = '{:8s} {:3d}'.format(s, stats[s]['n'])
    for m in ['CD16a_rmsd', 'CD16a_self', 'COM', 'BSA', 'dG']:
        mu, sd = stats[s][m]
        line += '{:>16s}'.format('{} ({})'.format(fmt(mu, 1), fmt(sd, 1)))
    print(line)

# does the replicate-to-replicate scatter in dG track a structural metric?


def rank(v):
    order = np.argsort(v)
    r = np.empty(len(v))
    r[order] = np.arange(len(v))
    return r


corr = []
dg = np.array([r['dG'] for r in rows])
for m in METRICS:
    x = np.array([r[m] for r in rows])
    if len(rows) > 3 and x.std() > 0 and dg.std() > 0:
        rho = float(np.corrcoef(rank(x), rank(dg))[0, 1])
    else:
        rho = float('nan')
    corr.append((m, rho))
with open(a.outdir + '/correlations.csv', 'w', newline='') as fh:
    w = csv.writer(fh)
    w.writerow(['metric', 'spearman_with_dG_across_all_replicates', 'n'])
    for m, rho in corr:
        w.writerow([m, round(rho, 3), len(rows)])
print('\nSpearman correlation of each metric with dG across all {} replicates '
      '(negative = larger metric goes with more favourable dG):'.format(len(rows)))
for m, rho in corr:
    print('  {:12s} {}'.format(m, fmt(rho)))

fig, axs = plt.subplots(2, 5, figsize=(18, 6.5), sharey=True)
colors = {'wt': 'black', 'reduced': 'tab:red', 'wtlike': 'tab:gray',
          'improved': 'tab:green'}
for ax, m in zip(axs.ravel(), METRICS):
    for r in rows:
        ax.scatter(r[m], r['dG'], color=colors.get(r['class'], 'b'), s=22)
    ax.set_xlabel(m)
for ax in axs[:, 0]:
    ax.set_ylabel('dG (kcal/mol)')
fig.suptitle('Each dot = one replicate. black WT, red reduced, grey WT-like, '
             'green improved', fontsize=10)
plt.tight_layout()
plt.savefig(a.outdir + '/scatter_dG_vs_structure.png', dpi=150)
plt.close()

if notes:
    print('\nNotes:')
    for n in notes:
        print('  - ' + n)
print('\nwrote {}/ (replicates.csv, variants.csv, correlations.csv, '
      'scatter_dG_vs_structure.png)'.format(a.outdir))
