#!/usr/bin/env python3
"""MM-GBSA analysis across variants and replicates, compared with WT.

usage (run from the prep folder):
    python analyze_mm.py --wt-extra /path/to/WT_perframe_s301.csv
    python analyze_mm.py --wt WT/rep1/perframe.csv WT/rep2/perframe.csv

Reads NAME/rep*/perframe.csv for each variant in variants.txt (commented-out
lines are still used when their folder exists). WT replicates come from
WT/rep*/perframe.csv plus any --wt-extra files, or from an explicit --wt list.

Writes mm_results/: replicates.csv variants.csv components.csv results.json
report.md dG_dotplot.png ddG_forest.png ddG_components.png timeseries_grid.png

The replicate is the unit of replication. Frames inside one run are correlated
and are used only for within-run diagnostics.
"""
import argparse
import csv
import glob
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

TERMS = ['VDWAALS', 'EEL', 'EGB', 'ESURF',
         'DELTA G gas', 'DELTA G solv', 'DELTA TOTAL']
COMP = ['VDWAALS', 'EEL', 'EGB', 'ESURF']
CLASS_SCORE = {'improved': -1, 'wtlike': 0, 'reduced': 1}
CLASS_COLOR = {'wt': 'black', 'reduced': 'tab:red',
               'wtlike': 'tab:gray', 'improved': 'tab:green'}

ap = argparse.ArgumentParser()
ap.add_argument('--variants', default='variants.txt')
ap.add_argument('--wt', nargs='*', default=None,
                help='explicit WT perframe.csv files')
ap.add_argument('--wt-extra', nargs='*', default=[],
                help='extra WT perframe.csv files')
ap.add_argument('--csv', default='perframe.csv')
ap.add_argument('--outdir', default='mm_results')
ap.add_argument('--start', type=int, default=301)
ap.add_argument('--step', type=int, default=5)
ap.add_argument('--dt', type=float, default=0.01,
                help='ns per trajectory frame')
ap.add_argument('--temp', type=float, default=310.15)
ap.add_argument('--nboot', type=int, default=5000)
a = ap.parse_args()

os.makedirs(a.outdir, exist_ok=True)
RT = 0.0019872041 * a.temp
rng = np.random.default_rng(12345)
warnings = []


# ------------------------------------------------------------------
# statistics helpers
# ------------------------------------------------------------------
def _tiny(x):
    return 1e-300 if abs(x) < 1e-300 else x


def _betacf(x, a_, b_, maxit=300, eps=3e-14):
    qab = a_ + b_
    qap = a_ + 1
    qam = a_ - 1
    c = 1.0
    d = _tiny(1 - qab * x / qap)
    d = 1 / d
    h = d
    for m in range(1, maxit + 1):
        m2 = 2 * m
        aa = m * (b_ - m) * x / ((qam + m2) * (a_ + m2))
        d = 1 / _tiny(1 + aa * d)
        c = _tiny(1 + aa / c)
        h *= d * c
        aa = -(a_ + m) * (qab + m) * x / ((a_ + m2) * (qap + m2))
        d = 1 / _tiny(1 + aa * d)
        c = _tiny(1 + aa / c)
        de = d * c
        h *= de
        if abs(de - 1) < eps:
            break
    return h


def betainc(a_, b_, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    lg = math.lgamma(a_ + b_) - math.lgamma(a_) - math.lgamma(b_)
    bt = math.exp(lg + a_ * math.log(x) + b_ * math.log(1 - x))
    if x < (a_ + 1) / (a_ + b_ + 2):
        return bt * _betacf(x, a_, b_) / a_
    return 1 - bt * _betacf(1 - x, b_, a_) / b_


def t_cdf(t, df):
    p = 0.5 * betainc(df / 2, 0.5, df / (df + t * t))
    return 1 - p if t > 0 else p


def t_ppf(q, df):
    """Upper-tail Student-t quantile by bisection."""
    lo, hi = 0.0, 2000.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if t_cdf(mid, df) < q:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def block_sem(x, nb=5):
    bm = np.array([b.mean() for b in np.array_split(x, nb)])
    return bm.std(ddof=1) / math.sqrt(nb)


def n_eff(x):
    x = x - x.mean()
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    acf = np.fft.irfft(f * np.conj(f))[:n]
    acf = acf / acf[0]
    tau = 0.5
    for k in range(1, n):
        if acf[k] <= 0:
            break
        tau += acf[k]
    return max(1.0, min(float(n), n / (2 * tau)))


def rankdata(v):
    v = np.asarray(v, float)
    order = np.argsort(v)
    r = np.empty(len(v))
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        r[order[i:j + 1]] = 0.5 * (i + j) + 1
        i = j + 1
    return r


def spearman(x, y):
    rx, ry = rankdata(x), rankdata(y)
    if rx.std() == 0 or ry.std() == 0:
        return float('nan')
    return float(np.corrcoef(rx, ry)[0, 1])


def holm(p):
    p = np.array(p, float)
    out = np.full(len(p), np.nan)
    ok = np.where(~np.isnan(p))[0]
    order = ok[np.argsort(p[ok])]
    m = len(order)
    run = 0.0
    for i, idx in enumerate(order):
        run = max(run, (m - i) * p[idx])
        out[idx] = min(1.0, run)
    return out


def isnan(x):
    return x is None or (isinstance(x, float) and math.isnan(x))


def f(x, p=2):
    return 'nan' if isnan(x) else '{:.{}f}'.format(x, p)


def g(x):
    return 'nan' if isnan(x) else '{:.3g}'.format(x)


# ------------------------------------------------------------------
# data loading
# ------------------------------------------------------------------
def load_delta(path):
    lines = open(path).read().splitlines()
    start = None
    for i, line in enumerate(lines):
        low = line.strip().lower()
        if 'energy terms' in low and low.startswith('delta'):
            start = i
            break
    if start is None:
        raise ValueError('no Delta Energy Terms section')
    hdr = [h.strip() for h in lines[start + 1].split(',')]
    rows = []
    for line in lines[start + 2:]:
        try:
            rows.append([float(x) for x in line.split(',')])
        except ValueError:
            break
    d = np.array(rows)
    return {t: d[:, hdr.index(t)] for t in TERMS if t in hdr}


def read_variants(path):
    if not os.path.exists(path):
        sys.exit(path + ' not found (run from the prep folder)')
    out = []
    for line in open(path):
        s = line.strip()
        if not s:
            continue
        commented = s.startswith('#')
        tok = s.lstrip('#').split()
        if len(tok) < 3 or tok[2] not in CLASS_SCORE:
            continue
        if commented and not os.path.isdir(tok[0]):
            continue
        out.append((tok[0], tok[2]))
    return out


systems = {}


def add_rep(name, cls, path):
    try:
        d = load_delta(path)
    except Exception as e:
        warnings.append('could not parse {}: {}'.format(path, e))
        return
    entry = systems.setdefault(name, {'class': cls, 'reps': []})
    entry['reps'].append({'path': path, 'data': d})


if a.wt is not None:
    wt_files = list(a.wt)
else:
    wt_files = sorted(glob.glob('WT/rep*/' + a.csv))
wt_files = [w for w in wt_files + list(a.wt_extra) if w]
if not wt_files:
    sys.exit('no WT files found: use --wt or --wt-extra')
for p in wt_files:
    add_rep('WT', 'wt', p)
for name, cls in read_variants(a.variants):
    files = sorted(glob.glob('{}/rep*/{}'.format(name, a.csv)))
    if not files:
        warnings.append(name + ': no ' + a.csv + ' found (not finished?)')
        continue
    for p in files:
        add_rep(name, cls, p)

if 'WT' not in systems:
    sys.exit('no WT replicate could be read')

nfr = set()
for s in systems.values():
    for r in s['reps']:
        nfr.add(len(r['data']['DELTA TOTAL']))
if len(nfr) > 1:
    warnings.append('replicates have different frame counts {}: the '
                    'MM-GBSA windows differ'.format(sorted(nfr)))

# ------------------------------------------------------------------
# per-replicate statistics
# ------------------------------------------------------------------
rep_rows = []
for name, s in systems.items():
    for i, r in enumerate(s['reps'], 1):
        x = r['data']['DELTA TOTAL']
        z = (x - x.mean()) / x.std(ddof=1)
        r['mean'] = {t: float(r['data'][t].mean()) for t in r['data']}
        r['bsem'] = float(block_sem(x))
        r['neff'] = float(n_eff(x))
        row = {'system': name, 'class': s['class'], 'rep': i,
               'file': r['path'], 'n_frames': len(x),
               'dG_mean': r['mean']['DELTA TOTAL'],
               'dG_sd_frames': float(x.std(ddof=1)),
               'dG_blockSEM': r['bsem'], 'n_eff': r['neff'],
               'n_outlier_frames_4sd': int((abs(z) > 4).sum())}
        for t in TERMS:
            if t in r['mean']:
                row[t + '_mean'] = r['mean'][t]
        rep_rows.append(row)

# ------------------------------------------------------------------
# per-system statistics
# ------------------------------------------------------------------
for name, s in systems.items():
    n = len(s['reps'])
    s['n'] = n
    s['m'] = {}
    s['sd'] = {}
    for t in TERMS:
        v = np.array([r['mean'][t] for r in s['reps'] if t in r['mean']])
        s['m'][t] = float(v.mean())
        s['sd'][t] = float(v.std(ddof=1)) if n > 1 else float('nan')
    s['rep_means'] = np.array([r['mean']['DELTA TOTAL'] for r in s['reps']])
    if n > 1:
        s['sem'] = s['sd']['DELTA TOTAL'] / math.sqrt(n)
        s['ci95'] = t_ppf(0.975, n - 1) * s['sem']
        s['sem_kind'] = 'between-replicate'
    else:
        s['sem'] = float(np.mean([r['bsem'] for r in s['reps']]))
        s['ci95'] = float('nan')
        s['sem_kind'] = 'within-run only (1 replicate)'
    if n < 3:
        warnings.append('{}: only {} replicate(s); weak statistics'.format(name, n))

wt = systems['WT']
names = [n for n in systems if n != 'WT']


def sem_term(sy, t):
    v = np.array([r['mean'][t] for r in sy['reps']])
    return v.std(ddof=1) / math.sqrt(len(v)) if len(v) > 1 else float('nan')


def ddg(name, t='DELTA TOTAL'):
    s = systems[name]
    d = s['m'][t] - wt['m'][t]
    if t == 'DELTA TOTAL':
        se = math.hypot(s['sem'], wt['sem'])
    else:
        se = math.hypot(sem_term(s, t), sem_term(wt, t))
    return d, se


res = {}
for name in names:
    s = systems[name]
    d, se = ddg(name)
    r = {'class': s['class'], 'n_reps': s['n'],
         'dG_mean': s['m']['DELTA TOTAL'],
         'dG_sd_between_reps': s['sd']['DELTA TOTAL'],
         'dG_sem': s['sem'], 'sem_kind': s['sem_kind'],
         'dG_ci95_halfwidth': s['ci95'], 'ddG': d, 'ddG_se': se,
         'p_welch': float('nan'), 'ddG_ci_lo': float('nan'),
         'ddG_ci_hi': float('nan')}
    n1, n2 = s['n'], wt['n']
    if n1 > 1 and n2 > 1:
        v1 = s['sd']['DELTA TOTAL'] ** 2 / n1
        v2 = wt['sd']['DELTA TOTAL'] ** 2 / n2
        if v1 + v2 > 0:
            tt = d / math.sqrt(v1 + v2)
            df = (v1 + v2) ** 2 / (v1 ** 2 / (n1 - 1) + v2 ** 2 / (n2 - 1))
            r['welch_t'] = tt
            r['welch_df'] = df
            r['p_welch'] = 2 * (1 - t_cdf(abs(tt), df))
            h = t_ppf(0.975, df) * math.sqrt(v1 + v2)
            r['ddG_ci_lo'] = d - h
            r['ddG_ci_hi'] = d + h
    wt_dg = wt['m']['DELTA TOTAL']
    if abs(wt_dg) > 1:
        r['ratio_dG'] = s['m']['DELTA TOTAL'] / wt_dg
        r['pct_change_dG'] = (r['ratio_dG'] - 1) * 100
    else:
        r['ratio_dG'] = float('nan')
        r['pct_change_dG'] = float('nan')
    r['implied_Kd_fold'] = math.exp(d / RT)
    if not math.isnan(r['ddG_ci_lo']):
        r['implied_Kd_fold_lo'] = math.exp(r['ddG_ci_lo'] / RT)
        r['implied_Kd_fold_hi'] = math.exp(r['ddG_ci_hi'] / RT)
    for t in COMP + ['DELTA G gas', 'DELTA G solv']:
        r['dd_' + t], r['dd_' + t + '_se'] = ddg(name, t)
    r['dd_polar_EEL+EGB'] = r['dd_EEL'] + r['dd_EGB']
    r['dd_nonpolar_VDW+ESURF'] = r['dd_VDWAALS'] + r['dd_ESURF']
    res[name] = r

holm_p = holm([res[n]['p_welch'] for n in names])
for name, p in zip(names, holm_p):
    res[name]['p_holm'] = float(p)

# ------------------------------------------------------------------
# ranking (1 = most favourable) with a bootstrap over replicate means
# ------------------------------------------------------------------
allnames = ['WT'] + names
means = np.array([systems[n]['m']['DELTA TOTAL'] for n in allnames])
rank = rankdata(means)
boot = np.zeros((a.nboot, len(allnames)))
for b in range(a.nboot):
    mb = []
    for n in allnames:
        rm = systems[n]['rep_means']
        mb.append(rng.choice(rm, size=len(rm)).mean())
    boot[b] = rankdata(mb)
for i, n in enumerate(allnames):
    info = {'rank': float(rank[i]),
            'rank_median_boot': float(np.median(boot[:, i])),
            'rank_lo90': float(np.percentile(boot[:, i], 5)),
            'rank_hi90': float(np.percentile(boot[:, i], 95))}
    if n == 'WT':
        wt.update(info)
    else:
        better = (boot[:, i] < boot[:, 0]).mean()
        tie = (boot[:, i] == boot[:, 0]).mean()
        info['P_more_favourable_than_WT'] = float(better + 0.5 * tie)
        res[n].update(info)

# ------------------------------------------------------------------
# agreement with the expected class order
# ------------------------------------------------------------------
cls_names = [n for n in names if res[n]['class'] in CLASS_SCORE]
agree = {}
if len(cls_names) > 2:
    sc = [CLASS_SCORE[res[n]['class']] for n in cls_names]
    dd = [res[n]['ddG'] for n in cls_names]
    agree['spearman_ddG_vs_class'] = spearman(dd, sc)
    pairs = []
    for i in range(len(cls_names)):
        for j in range(len(cls_names)):
            if sc[i] < sc[j]:
                pairs.append((i, j))
    ok = sum(1 for i, j in pairs if dd[i] < dd[j])
    agree['cross_class_pairs'] = len(pairs)
    agree['cross_class_concordant'] = int(ok)
    agree['concordance_fraction'] = ok / len(pairs) if pairs else float('nan')
    checks = {}
    for n in cls_names:
        c, d = res[n]['class'], res[n]['ddG']
        if c == 'wtlike':
            checks[n] = 'n/a (wt-like)'
        elif (c == 'improved' and d < 0) or (c == 'reduced' and d > 0):
            checks[n] = 'as expected'
        else:
            checks[n] = 'OPPOSITE'
    agree['sign_check'] = checks

order = sorted(names, key=lambda n: res[n]['dG_mean'])

# ------------------------------------------------------------------
# tables
# ------------------------------------------------------------------
with open(a.outdir + '/replicates.csv', 'w', newline='') as fh:
    keys = []
    for row in rep_rows:
        for k in row:
            if k not in keys:
                keys.append(k)
    w = csv.DictWriter(fh, fieldnames=keys)
    w.writeheader()
    w.writerows(rep_rows)

cols = ['system', 'class', 'n_reps', 'dG_mean', 'dG_sd_between_reps',
        'dG_sem', 'sem_kind', 'dG_ci95_halfwidth', 'ddG', 'ddG_se',
        'ddG_ci_lo', 'ddG_ci_hi', 'p_welch', 'p_holm', 'ratio_dG',
        'pct_change_dG', 'implied_Kd_fold', 'implied_Kd_fold_lo',
        'implied_Kd_fold_hi', 'rank', 'rank_median_boot', 'rank_lo90',
        'rank_hi90', 'P_more_favourable_than_WT']
wt_row = {'system': 'WT', 'class': 'wt', 'n_reps': wt['n'],
          'dG_mean': wt['m']['DELTA TOTAL'],
          'dG_sd_between_reps': wt['sd']['DELTA TOTAL'],
          'dG_sem': wt['sem'], 'sem_kind': wt['sem_kind'],
          'dG_ci95_halfwidth': wt['ci95'], 'ddG': 0.0, 'ratio_dG': 1.0,
          'pct_change_dG': 0.0, 'implied_Kd_fold': 1.0,
          'rank': wt['rank'], 'rank_median_boot': wt['rank_median_boot'],
          'rank_lo90': wt['rank_lo90'], 'rank_hi90': wt['rank_hi90']}
with open(a.outdir + '/variants.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, extrasaction='ignore')
    w.writeheader()
    w.writerow(wt_row)
    for n in order:
        row = dict(res[n])
        row['system'] = n
        w.writerow(row)

ccols = (['system', 'class']
         + ['dd_' + t for t in COMP + ['DELTA G gas', 'DELTA G solv']]
         + ['dd_' + t + '_se' for t in COMP]
         + ['dd_polar_EEL+EGB', 'dd_nonpolar_VDW+ESURF'])
with open(a.outdir + '/components.csv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=ccols, extrasaction='ignore')
    w.writeheader()
    for n in order:
        row = dict(res[n])
        row['system'] = n
        w.writerow(row)

wt_json = {k: wt[k] for k in ('n', 'm', 'sd', 'sem', 'ci95', 'rank')}
out_json = {'wt': wt_json, 'variants': res,
            'agreement_with_expected_class_order': agree,
            'warnings': warnings}
with open(a.outdir + '/results.json', 'w') as fh:
    json.dump(out_json, fh, indent=1, default=str)


# ------------------------------------------------------------------
# plots
# ------------------------------------------------------------------
def style(ax):
    ax.grid(axis='y', alpha=0.25)


# 1. replicate means per system, ordered by mean dG
srt = ['WT'] + order
fig, ax = plt.subplots(figsize=(1.1 * len(allnames) + 3, 4.8))
wm, ws = wt['m']['DELTA TOTAL'], wt['sem']
ax.axhspan(wm - ws, wm + ws, color='k', alpha=0.1)
ax.axhline(wm, color='k', lw=0.8, ls='--')
for i, n in enumerate(srt):
    s = systems[n]
    c = CLASS_COLOR.get(s['class'], 'b')
    k = len(s['rep_means'])
    jitter = rng.uniform(-0.08, 0.08, k)
    ax.scatter(np.full(k, i) + jitter, s['rep_means'], color=c,
               zorder=3, s=28)
    ax.errorbar(i, s['m']['DELTA TOTAL'], yerr=s['sem'], color=c,
                capsize=5, lw=1.6, zorder=2)
ax.set_xticks(range(len(srt)))
ax.set_xticklabels(['{}\n(n={})'.format(n, systems[n]['n']) for n in srt])
ax.set_ylabel('MM-GBSA DELTA TOTAL (kcal/mol)')
style(ax)
ax.set_title('Replicate means (dots), mean +/- SEM (bar); shaded = WT +/- SEM\n'
             'black WT, green improved, grey WT-like, red reduced', fontsize=9)
plt.tight_layout()
plt.savefig(a.outdir + '/dG_dotplot.png', dpi=150)
plt.close()

# 2. forest plot of ddG
fig, ax = plt.subplots(figsize=(7, 0.55 * len(order) + 1.8))
rev = order[::-1]
for i, n in enumerate(rev):
    r = res[n]
    c = CLASS_COLOR.get(r['class'], 'b')
    ax.scatter(r['ddG'], i, color=c, zorder=3)
    if not math.isnan(r['ddG_ci_lo']):
        ax.hlines(i, r['ddG_ci_lo'], r['ddG_ci_hi'], color=c, lw=2)
    else:
        ax.hlines(i, r['ddG'] - r['ddG_se'], r['ddG'] + r['ddG_se'],
                  color=c, lw=1, ls=':')
ax.axvline(0, color='k', lw=0.8)
ax.set_yticks(range(len(rev)))
ax.set_yticklabels(['{}  [{}]'.format(n, res[n]['class']) for n in rev])
ax.set_xlabel('ddG vs WT (kcal/mol); negative = more favourable than WT')
ax.grid(axis='x', alpha=0.25)
ax.set_title('ddG with 95% CI (solid), or +/- SE if n=1 (dotted)', fontsize=9)
plt.tight_layout()
plt.savefig(a.outdir + '/ddG_forest.png', dpi=150)
plt.close()

# 3. energy terms behind each ddG
fig, ax = plt.subplots(figsize=(1.0 * len(order) + 3, 4.5))
wdt = 0.2
for k, t in enumerate(COMP):
    xs = np.arange(len(order)) + (k - 1.5) * wdt
    ax.bar(xs, [res[n]['dd_' + t] for n in order], wdt, label=t)
ax.scatter(range(len(order)), [res[n]['ddG'] for n in order],
           color='k', zorder=3, marker='D', label='DELTA TOTAL')
ax.axhline(0, color='k', lw=0.8)
ax.set_xticks(range(len(order)))
ax.set_xticklabels(order)
ax.set_ylabel('difference from WT (kcal/mol)')
ax.legend(fontsize=7, ncol=3)
style(ax)
ax.set_title('Energy terms behind each ddG (EEL and EGB largely cancel)',
             fontsize=9)
plt.tight_layout()
plt.savefig(a.outdir + '/ddG_components.png', dpi=150)
plt.close()

# 4. per-frame time series, one panel per system
nsys = len(allnames)
nc = 3
nr = math.ceil(nsys / nc)
fig, axs = plt.subplots(nr, nc, figsize=(4.2 * nc, 2.6 * nr),
                        sharey=True, squeeze=False)
for ax, n in zip(axs.ravel(), allnames):
    for i, r in enumerate(systems[n]['reps'], 1):
        x = r['data']['DELTA TOTAL']
        t = (a.start - 1 + a.step * np.arange(len(x))) * a.dt
        w_ = max(5, len(x) // 10)
        pad = np.pad(x, (w_ // 2, w_ - 1 - w_ // 2), mode='edge')
        run = np.convolve(pad, np.ones(w_) / w_, mode='valid')
        col = 'C{}'.format(i - 1)
        ax.plot(t, x, lw=0.4, alpha=0.3, color=col)
        ax.plot(t, run, lw=1.4, color=col, label='rep{}'.format(i))
    ax.set_title('{} [{}]'.format(n, systems[n]['class']), fontsize=9)
    ax.legend(fontsize=6)
    ax.set_xlabel('time (ns)')
for ax in axs.ravel()[nsys:]:
    ax.axis('off')
axs[0][0].set_ylabel('DELTA TOTAL (kcal/mol)')
plt.tight_layout()
plt.savefig(a.outdir + '/timeseries_grid.png', dpi=150)
plt.close()

# ------------------------------------------------------------------
# readable report
# ------------------------------------------------------------------
L = ['# MM-GBSA summary', '']
L.append('WT: {} replicate(s), DELTA TOTAL {} kcal/mol, SEM {} ({}).'.format(
    wt['n'], f(wt['m']['DELTA TOTAL']), f(wt['sem']), wt['sem_kind']))
L.append('Window: startframe {}, interval {}; {} frames per replicate. '
         'T = {} K.'.format(a.start, a.step, sorted(nfr), a.temp))
L += ['', '## Ranking by DELTA TOTAL (1 = most favourable)', '']
L.append('| rank | system | class | n | dG (kcal/mol) | SD between reps | '
         'ddG vs WT | 95% CI | p (Welch) | p (Holm) | rank 90% range | '
         'P(better than WT) |')
L.append('|---|---|---|---|---|---|---|---|---|---|---|---|')
L.append('| {} | WT | wt | {} | {} | {} | 0 | | | | {}-{} | |'.format(
    f(wt['rank'], 0), wt['n'], f(wt['m']['DELTA TOTAL']),
    f(wt['sd']['DELTA TOTAL']), f(wt['rank_lo90'], 0), f(wt['rank_hi90'], 0)))
for n in order:
    r = res[n]
    ci = '{} to {}'.format(f(r['ddG_ci_lo']), f(r['ddG_ci_hi']))
    rk = '{}-{}'.format(f(r['rank_lo90'], 0), f(r['rank_hi90'], 0))
    L.append('| {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} | {} |'.format(
        f(r['rank'], 0), n, r['class'], r['n_reps'], f(r['dG_mean']),
        f(r['dG_sd_between_reps']), f(r['ddG']), ci, f(r['p_welch'], 3),
        f(r['p_holm'], 3), rk, f(r['P_more_favourable_than_WT'])))

L += ['', '## Fold-style comparison with WT', '']
L.append('| variant | dG ratio (mut/WT) | change in dG (%) | '
         'implied Kd fold (mut/WT) | range from 95% CI |')
L.append('|---|---|---|---|---|')
for n in order:
    r = res[n]
    if 'implied_Kd_fold_lo' in r:
        rg = '{} to {}'.format(g(r['implied_Kd_fold_lo']),
                               g(r['implied_Kd_fold_hi']))
    else:
        rg = 'n/a'
    L.append('| {} | {} | {} | {} | {} |'.format(
        n, f(r['ratio_dG'], 3), f(r['pct_change_dG'], 1),
        g(r['implied_Kd_fold']), rg))
L += ['',
      'ratio > 1 means a more negative (more favourable) dG than WT. '
      '"Implied Kd fold" = exp(ddG/RT); >1 means weaker binding. It is NOT '
      'a quantitative affinity prediction: entropy is omitted, MM-GBSA '
      'energies are inflated in magnitude, and the exponential amplifies '
      'noise. Use the sign, the ranking and the CI, not the size of the fold.',
      '', '## Energy terms behind ddG (kcal/mol, difference from WT)', '']
L.append('| variant | VDWAALS | EEL | EGB | ESURF | polar (EEL+EGB) | '
         'nonpolar (VDW+ESURF) | total |')
L.append('|---|---|---|---|---|---|---|---|')
for n in order:
    r = res[n]
    L.append('| {} | {} | {} | {} | {} | {} | {} | {} |'.format(
        n, f(r['dd_VDWAALS']), f(r['dd_EEL']), f(r['dd_EGB']),
        f(r['dd_ESURF']), f(r['dd_polar_EEL+EGB']),
        f(r['dd_nonpolar_VDW+ESURF']), f(r['ddG'])))
if agree:
    L += ['', '## Agreement with the expected class order '
          '(improved < WT-like < reduced in dG)', '']
    L.append('- Spearman correlation of ddG with class: {}'.format(
        f(agree['spearman_ddG_vs_class'])))
    L.append('- Cross-class pairs ordered as expected: {} of {} ({}%)'.format(
        agree['cross_class_concordant'], agree['cross_class_pairs'],
        f(100 * agree['concordance_fraction'], 0)))
    sc_txt = ', '.join('{} {}'.format(n, v)
                       for n, v in agree['sign_check'].items())
    L.append('- Sign check: ' + sc_txt)
L += ['', '## Notes', '',
      '- The replicate is the unit of replication. SD and SEM are across '
      'replicate means; frames inside a replicate are correlated (see n_eff '
      'in replicates.csv).',
      '- With 3 replicates per system the Welch test has little power, so a '
      'large p does not mean no effect. Use the CI and the consistency '
      'across replicates (timeseries_grid.png). Several variants are tested, '
      'so read the Holm column rather than raw p.',
      '- The bootstrap rank range resamples replicate means and is crude '
      'with n=3.',
      '- Entropy is omitted, receptor glycans are single-NAG stubs, and '
      'GLYCAM under OPC3 is not validated.']
if warnings:
    L += ['', '## Warnings', ''] + ['- ' + w for w in warnings]
with open(a.outdir + '/report.md', 'w') as fh:
    fh.write('\n'.join(L) + '\n')
print('\n'.join(L))
print('\nwrote {}/ (report.md, variants.csv, components.csv, replicates.csv, '
      'results.json, 4 png)'.format(a.outdir))
