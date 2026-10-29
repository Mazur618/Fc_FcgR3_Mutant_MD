#!/usr/bin/env python3
"""usage: python plot_compare.py --wt path/FINAL_WT_s301.dat [more WT FINAL files]
Dot plot of DELTA TOTAL per replicate for WT and every variant in variants.txt,
coloured by class; saves compare_dG.png."""
import argparse, glob, re
import numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ap = argparse.ArgumentParser()
ap.add_argument('--wt', nargs='+', required=True)
a = ap.parse_args()
pat = re.compile(r'DELTA TOTAL\s+(-?\d+\.\d+)')
def dg(p):
    m = pat.findall(open(p).read()); return float(m[-1]) if m else None
rows = [('WT', 'wt', [dg(f) for f in a.wt])]
for l in open('variants.txt'):
    if not l.strip() or l.startswith('#'): continue
    n, s, c = l.split()[:3]
    rows.append((n, c, [dg(f) for f in sorted(glob.glob(f'{n}/rep*/FINAL.dat'))]))
order = {'wt': 0, 'reduced': 1, 'wtlike': 2, 'improved': 3}
rows.sort(key=lambda r: (order.get(r[1], 9), r[0]))
col = {'wt': 'k', 'reduced': 'tab:red', 'wtlike': 'tab:gray', 'improved': 'tab:green'}
fig, ax = plt.subplots(figsize=(9, 4.5))
for i, (n, c, v) in enumerate(rows):
    v = np.array([x for x in v if x is not None])
    if len(v) == 0: continue
    ax.scatter([i] * len(v), v, color=col.get(c, 'b'), zorder=3)
    ax.hlines(v.mean(), i - 0.3, i + 0.3, color=col.get(c, 'b'))
    print(f'{n:8s} {c:9s} n={len(v)} mean {v.mean():7.2f}' + (f' sd {v.std(ddof=1):5.2f}' if len(v) > 1 else ''))
ax.set_xticks(range(len(rows))); ax.set_xticklabels([r[0] for r in rows])
ax.set_ylabel('MM-GBSA DELTA TOTAL (kcal/mol)')
ax.set_title('One dot per replicate; bar = mean. Black WT, red reduced, grey WT-like, green improved')
plt.tight_layout(); plt.savefig('compare_dG.png', dpi=150)
print('saved compare_dG.png')
