#!/usr/bin/env python3
"""Plot existing whole-partner MM/GBSA alongside completed structural analysis.

No energy calculations are performed. Rows in perframe.csv are matched in order
to the configured snapshots, assuming the approved original MM/GBSA window.
Frame counts and ordinal IDs are checked; CSV alone cannot prove window identity.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

from ..io import write_csv
from ..structure.plot_structure import load_runs


TERMS = ['VDWAALS', 'EEL', 'EGB', 'ESURF', 'DELTA TOTAL']


def load_energy(path):
    with Path(path).open() as fh:
        for line in fh:
            if line.strip().lower().startswith('delta') and 'energy terms' in line.lower():
                break
        else:
            raise ValueError(f'no Delta Energy Terms section: {path}')
        header = [v.strip() for v in next(csv.reader(fh))]
        rows = []
        for row in csv.reader(fh):
            if not row or not any(v.strip() for v in row):
                break
            values = dict(zip(header, map(float, row)))
            if len(row) != len(header):
                raise ValueError(f'ragged energy row: {path}')
            rows.append(values)
    if not rows:
        raise ValueError(f'empty Delta Energy Terms section: {path}')
    for i, row in enumerate(rows):
        if row['Frame #'] != i:
            raise ValueError(f'nonconsecutive snapshot ordinal: {path}')
        for term in TERMS:
            if term not in row:
                raise ValueError(f'missing {term}: {path}')
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    cfg, runs = load_runs(args.runs)
    import numpy as np
    from scipy.stats import t
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    systems = cfg['systems']
    rows, sources = [], []
    for run in runs:
        source = Path(run['meta']['input']['mmgbsa_csv'])
        energy = load_energy(source)
        if len(energy) != run['meta']['n_frames']:
            raise ValueError(f'energy and structural windows have different counts: {source}')
        values = np.array([[r[k] for k in TERMS] for r in energy])
        if not np.all(np.isfinite(values)):
            raise ValueError(f'nonfinite energies: {source}')
        if np.max(np.abs(values[:, :4].sum(axis=1) - values[:, 4])) > 0.02:
            raise ValueError(f'four-component sum does not match total: {source}')
        row = dict(system=run['meta']['input']['system'], replicate=run['meta']['input']['replicate'])
        row.update({k: float(values[:, i].mean()) for i, k in enumerate(TERMS)})
        row['EEL+EGB'] = row['EEL'] + row['EGB']
        for metric in ('contact_residue_pairs', 'salt_bridge_residue_pairs'):
            row[metric] = float(np.mean([f[metric] for f in run['frames']]))
        rows.append(row)
        sources.append(dict(path=str(source), sha256=hashlib.sha256(source.read_bytes()).hexdigest()))
    def vals(system, term):
        return np.array([r[term] for r in rows if r['system'] == system])
    def mean(system, term):
        return float(vals(system, term).mean())
    stats = []
    wt = vals('WT', 'DELTA TOTAL')
    for s in systems:
        v = vals(s, 'DELTA TOTAL')
        row = dict(system=s, n_replicates=len(v), mean_kcal_mol=float(v.mean()),
                   sd_between_replicates=float(v.std(ddof=1)), ddG_vs_WT=float(v.mean() - wt.mean()))
        if s != 'WT':
            if min(len(v), len(wt)) < 2:
                raise ValueError('Welch intervals require at least two replicates per system')
            a, b = v.var(ddof=1) / len(v), wt.var(ddof=1) / len(wt)
            if a + b == 0:
                raise ValueError('Welch interval is undefined with zero variance in both systems')
            df = (a + b) ** 2 / (a * a / (len(v) - 1) + b * b / (len(wt) - 1))
            half = float(t.ppf(0.975, df) * np.sqrt(a + b))
            row.update(ci95_lo=row['ddG_vs_WT'] - half, ci95_hi=row['ddG_vs_WT'] + half,
                       p_welch=float(2 * t.sf(abs(row['ddG_vs_WT']) / np.sqrt(a + b), df)))
        stats.append(row)
    comparisons = [r for r in stats if r['system'] != 'WT']
    adjusted = 0.
    for i, r in enumerate(sorted(comparisons, key=lambda r: r['p_welch'])):
        adjusted = max(adjusted, (len(comparisons) - i) * r['p_welch'])
        r['p_holm'] = min(1., adjusted)
    args.out.mkdir(parents=True, exist_ok=False)
    write_csv(args.out / 'replicates.csv', list(rows[0]), rows)
    write_csv(args.out / 'comparisons.csv', ['system', 'n_replicates', 'mean_kcal_mol',
              'sd_between_replicates', 'ddG_vs_WT', 'ci95_lo', 'ci95_hi', 'p_welch', 'p_holm'], stats)
    def save(fig, name):
        for ext in ('png', 'pdf'):
            fig.savefig(args.out / f'{name}.{ext}', dpi=200, bbox_inches='tight')
        plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), gridspec_kw={'width_ratios': [1.3, 1]})
    for i, s in enumerate(systems):
        v = vals(s, 'DELTA TOTAL')
        axes[0].scatter(i + np.linspace(-0.12, 0.12, len(v)), v, color='0.3', s=28)
        axes[0].errorbar(i, v.mean(), yerr=v.std(ddof=1), color='black', capsize=4, marker='_')
    axes[0].set_xticks(range(len(systems)), systems, rotation=55, ha='right')
    axes[0].set_ylabel('Whole-partner MM/GBSA (kcal/mol)')
    axes[0].set_title('Replicate means; bars = between-replicate SD')
    for i, r in enumerate(comparisons):
        axes[1].plot([r['ci95_lo'], r['ci95_hi']], [i, i], color='0.35', lw=2)
        axes[1].scatter(r['ddG_vs_WT'], i, color='black')
    axes[1].axvline(0, color='0.5', ls='--', lw=1)
    axes[1].set_yticks(range(len(comparisons)), [r['system'] for r in comparisons])
    axes[1].invert_yaxis()
    axes[1].set_xlabel('Variant − WT (kcal/mol)')
    axes[1].set_title('Difference with unadjusted 95% Welch CI')
    fig.tight_layout()
    save(fig, 'binding_comparison')
    variants = [s for s in systems if s != 'WT']
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for ax, terms, title in [(axes[0], ['EEL', 'EGB'], 'Direct electrostatics and GB solvation'),
                             (axes[1], ['EEL+EGB', 'VDWAALS', 'ESURF'], 'Net polar and nonpolar components')]:
        width = 0.75 / len(terms)
        for i, term in enumerate(terms):
            delta = [mean(s, term) - mean('WT', term) for s in variants]
            ax.bar(np.arange(len(variants)) + (i - (len(terms) - 1) / 2) * width,
                   delta, width, label=term)
        ax.axhline(0, color='0.4', lw=0.8)
        ax.set_xticks(range(len(variants)), variants, rotation=55, ha='right')
        ax.set_ylabel('Variant − WT (kcal/mol)')
        ax.set_title(title)
        ax.legend(fontsize=8)
    fig.tight_layout()
    save(fig, 'energy_components')
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    colors = plt.get_cmap('tab10')
    for ax, metric in zip(axes, ['contact_residue_pairs', 'salt_bridge_residue_pairs']):
        for i, s in enumerate(systems):
            selected = [r for r in rows if r['system'] == s]
            ax.scatter([r[metric] for r in selected], [r['DELTA TOTAL'] for r in selected],
                       label=s, color=colors(i), s=32)
        ax.set_xlabel('Mean ' + metric.replace('_', ' '))
    axes[0].set_ylabel('Whole-partner MM/GBSA (kcal/mol)')
    axes[1].legend(bbox_to_anchor=(1.03, 1), loc='upper left', fontsize=8)
    fig.suptitle('Each point is one replicate; exploratory structural–energetic comparison')
    fig.tight_layout()
    save(fig, 'contacts_vs_binding')
    with (args.out / 'provenance.json').open('x') as fh:
        json.dump(dict(protocol=cfg, energy_sources=sources, structural_runs=[r['path'] for r in runs],
                       caveat='Whole-partner energies; entropy omitted. CSV row order assumes original approved window. '
                              'Confidence intervals are unadjusted; p_holm adjusts comparisons. Component bars are descriptive.'), fh, indent=2)
    print('system mean SD ddG [95% CI] p_Holm')
    for r in stats:
        print(r)
    print(f'Wrote MM/GBSA figures and tables -> {args.out}')


if __name__ == '__main__':
    main()
