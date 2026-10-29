#!/usr/bin/env python3
"""Overlay stability metrics and MM-GBSA per-frame energies for one or more runs.

Each run folder should contain the .dat files from the cpptraj analysis
(rmsd.dat, rmsd_self.dat, comdist.dat, rg.dat, sasa.dat) and the per-frame
MM-GBSA csv (default WT_perframe.csv).

Usage:
    python stability_plot.py                 # current folder
    python stability_plot.py run1 run2       # overlay two runs
    python stability_plot.py run1 run2 --start 157 --step 5 --dt 0.01
"""
import argparse, os
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ap = argparse.ArgumentParser()
ap.add_argument('dirs', nargs='*', default=['.'], help='run folders to overlay')
ap.add_argument('--dt', type=float, default=0.01, help='ns per trajectory frame (10 ps = 0.01)')
ap.add_argument('--start', type=int, default=157, help='MM-GBSA startframe (1-based)')
ap.add_argument('--step', type=int, default=5, help='MM-GBSA interval')
ap.add_argument('--mm', default='WT_perframe.csv', help='per-frame MM-GBSA csv name')
ap.add_argument('--win', type=int, default=25, help='running-mean window (points)')
ap.add_argument('-o', '--out', default='stability_all.png')
a = ap.parse_args()


def load(path):
    if not os.path.exists(path):
        print(f'  missing: {path}')
        return None, None
    names = open(path).readline().replace('#', '').split()
    return names, np.loadtxt(path)


def smooth(y, w):
    if len(y) < w:
        return y
    p = np.pad(y, (w // 2, w - 1 - w // 2), mode='edge')
    return np.convolve(p, np.ones(w) / w, mode='valid')


def parse_delta(path):
    if not os.path.exists(path):
        print(f'  missing: {path}')
        return None, None
    L = open(path).read().splitlines()
    s = next((i for i, l in enumerate(L)
              if 'Energy Terms' in l and l.strip().lower().startswith('delta')), None)
    if s is None:
        return None, None
    hdr = [h.strip() for h in L[s + 1].split(',')]
    rows = []
    for l in L[s + 2:]:
        try:
            rows.append([float(x) for x in l.split(',')])
        except ValueError:
            break
    return hdr, np.array(rows)


panels = [
    ('RMSD after fit on Fc CH3', 'A'),
    ('CD16a / glycan RMSD (fit on CH3)', 'A'),
    ('Self-fit RMSD (internal stability)', 'A'),
    ('Fc-CD16a COM distance', 'A'),
    ('Radius of gyration', 'A'),
    ('Buried surface area per partner', 'A^2'),
    ('MM-GBSA DELTA TOTAL', 'kcal/mol'),
    ('EEL, EGB and EEL+EGB', 'kcal/mol'),
]
fig, axs = plt.subplots(4, 2, figsize=(13, 15))
axs = axs.ravel()
for ax_, (title, ylab) in zip(axs, panels):
    ax_.set_title(title, fontsize=10)
    ax_.set_ylabel(ylab)
    ax_.set_xlabel('time (ns)')


def put(ax, t, y, color, label=None, ls='-'):
    ax.plot(t, y, color=color, alpha=0.25, lw=0.6)
    ax.plot(t, smooth(y, a.win), color=color, lw=1.4, label=label, ls=ls)


def blocks(tag, name, y):
    m = [b.mean() for b in np.array_split(y, 4)]
    print(f'  {tag:12s} {name:16s} blocks: ' + ' '.join(f'{x:9.2f}' for x in m))


colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
for k, d in enumerate(a.dirs):
    c = colors[k % len(colors)]
    tag = os.path.basename(os.path.abspath(d))
    print(f'== {tag}')

    n, r = load(f'{d}/rmsd.dat')
    if r is not None:
        t = (r[:, 0] - 1) * a.dt
        for nm, pi in (('CH3_fit', 0), ('Fc_CH2', 0), ('CD16a', 1), ('gly_Fc', 1), ('gly_stub', 1)):
            j = n.index(nm)
            ls = '-' if nm in ('CH3_fit', 'CD16a') else '--' if nm in ('Fc_CH2', 'gly_Fc') else ':'
            put(axs[pi], t, r[:, j], c, f'{tag} {nm}', ls)
            blocks(tag, nm, r[:, j])

    n, r = load(f'{d}/rmsd_self.dat')
    if r is not None:
        t = (r[:, 0] - 1) * a.dt
        for nm, ls in (('CD16a_self', '-'), ('Fc_self', '--')):
            j = n.index(nm)
            put(axs[2], t, r[:, j], c, f'{tag} {nm}', ls)
            blocks(tag, nm, r[:, j])

    n, r = load(f'{d}/comdist.dat')
    if r is not None:
        put(axs[3], (r[:, 0] - 1) * a.dt, r[:, 1], c, tag)
        blocks(tag, 'COM', r[:, 1])

    n, r = load(f'{d}/rg.dat')
    if r is not None:
        put(axs[4], (r[:, 0] - 1) * a.dt, r[:, 1], c, tag)
        blocks(tag, 'Rg', r[:, 1])

    n, r = load(f'{d}/sasa.dat')
    if r is not None:
        f, cx, fc, rc = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
        bsa = (fc + rc - cx) / 2
        ok = (bsa > 0) & (bsa <= np.minimum(fc, rc))  # drop physically impossible frames only
        print(f'  {tag:12s} BSA: dropped {int((~ok).sum())} impossible frame(s)')
        put(axs[5], (f[ok] - 1) * a.dt, bsa[ok], c, tag)
        blocks(tag, 'BSA', bsa[ok])

    h, m = parse_delta(f'{d}/{a.mm}')
    if m is not None:
        t = (a.start - 1 + a.step * np.arange(len(m))) * a.dt
        tot = m[:, h.index('DELTA TOTAL')]
        put(axs[6], t, tot, c, tag)
        blocks(tag, 'DELTA TOTAL', tot)
        eel, egb = m[:, h.index('EEL')], m[:, h.index('EGB')]
        put(axs[7], t, eel, c, f'{tag} EEL', ':')
        put(axs[7], t, egb, c, f'{tag} EGB', '--')
        put(axs[7], t, eel + egb, c, f'{tag} EEL+EGB', '-')
        blocks(tag, 'EEL+EGB', eel + egb)

for ax_ in axs:
    if ax_.get_legend_handles_labels()[0]:
        ax_.legend(fontsize=6, ncol=2)
plt.tight_layout()
plt.savefig(a.out, dpi=150)
print('saved', a.out)
