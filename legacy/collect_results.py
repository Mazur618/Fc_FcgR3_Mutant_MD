#!/usr/bin/env python3
"""usage: python collect_results.py --wt path/FINAL_WT.dat [more WT files ...]
Reads NAME/rep*/FINAL.dat for each variant in variants.txt and prints
per-replicate DELTA TOTAL, the replicate mean/SD, and ddG versus WT."""
import argparse, glob, re
import numpy as np
ap = argparse.ArgumentParser()
ap.add_argument('--wt', nargs='+', required=True, help='FINAL dat files of WT replicates')
a = ap.parse_args()
pat = re.compile(r'DELTA TOTAL\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)')

def dg(path):
    m = pat.findall(open(path).read())
    return float(m[-1][0]) if m else None

def summ(vals):
    v = np.array([x for x in vals if x is not None])
    return v, (v.mean() if len(v) else np.nan), (v.std(ddof=1) if len(v) > 1 else np.nan)

wt, wm, ws = summ([dg(f) for f in a.wt])
print(f'WT: reps {np.round(wt, 1).tolist()}  mean {wm:.2f}  sd {ws:.2f}  (n={len(wt)})')
print(f'{"variant":8s} {"class":9s} {"n":>2s} {"mean":>8s} {"sd":>6s} {"ddG vs WT":>10s}  replicate values')
for line in open('variants.txt'):
    if not line.strip() or line.startswith('#'):
        continue
    name, spec, cls = line.split()[:3]
    v, m, s = summ([dg(f) for f in sorted(glob.glob(f'{name}/rep*/FINAL.dat'))])
    print(f'{name:8s} {cls:9s} {len(v):2d} {m:8.2f} {s:6.2f} {m - wm:10.2f}  {np.round(v, 1).tolist()}')
print('\nA ddG is only meaningful if it is larger than the replicate SD (WT and variant).')
