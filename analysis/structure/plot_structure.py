#!/usr/bin/env python3
"""Compare completed analysis runs; never read trajectories or replace output.

All supplied runs must share a protocol/reference and include every configured
replicate. Interface plots use a union of pairs passing the display threshold in
any system. Absent pairs count as zero only in successfully completed runs.
"""
import argparse
import csv
import json
import math
from pathlib import Path

from ..io import write_csv


METRIC_LABELS = {
    'CH3_common_reference_rmsd_A': ('Fc CH3 after CH3 alignment', 'Backbone RMSD (Å)'),
    'CH2_common_reference_rmsd_A': ('Fc CH2 after CH3 alignment', 'Backbone RMSD (Å)'),
    'receptor_common_reference_rmsd_A': ('Receptor after CH3 alignment', 'Backbone RMSD (Å)'),
    'Fc_self_fit_rmsd_A': ('Fc internal deformation', 'Self-fit backbone RMSD (Å)'),
    'receptor_self_fit_rmsd_A': ('Receptor internal deformation', 'Self-fit backbone RMSD (Å)'),
    'partner_COM_distance_A': ('Glycosylated partner separation', 'COM distance (Å)'),
    'complex_Rg_A': ('Complex size', 'Radius of gyration (Å)'),
    'contact_residue_pairs': ('Interface heavy-atom contacts', 'Mean residue-pair count'),
    'salt_bridge_residue_pairs': ('Interface salt-bridge candidates', 'Mean residue-pair count'),
}


def load_runs(paths):
    runs, seen = [], set()
    for path in paths:
        path = Path(path)
        meta = json.loads((path / 'complete.json').read_text())
        key = (meta['input']['system'], meta['input']['replicate'])
        source = meta['input']['trajectory']
        if key in seen or source in seen:
            raise ValueError(f'duplicate run or trajectory: {path}')
        seen.update((key, source))
        if runs and any(meta[k] != runs[0]['meta'][k] for k in
                        ('protocol_sha256', 'reference_sha256', 'n_frames')):
            raise ValueError(f'incompatible analysis protocol, reference, or window: {path}')
        with (path / 'pairs.csv').open() as fh:
            pairs = {(r['kind'], r['Fc_label'], r['receptor_label']): float(r['occupancy'])
                     for r in csv.DictReader(fh)}
        with (path / 'frames.csv').open() as fh:
            frames = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(fh)]
        if len(frames) != meta['n_frames']:
            raise ValueError(f'incomplete frames table: {path}')
        runs.append(dict(meta=meta, pairs=pairs, frames=frames, path=str(path)))
    if not runs:
        raise ValueError('no completed runs supplied')
    cfg = runs[0]['meta']['protocol']
    expected = {(s, r) for s in cfg['systems'] for r in cfg['replicates']}
    actual = {(r['meta']['input']['system'], r['meta']['input']['replicate']) for r in runs}
    if actual != expected:
        raise ValueError(f'input set differs from protocol; missing={sorted(expected - actual)}')
    return cfg, runs


def system_runs(runs, system):
    return sorted([r for r in runs if r['meta']['input']['system'] == system],
                  key=lambda r: r['meta']['input']['replicate'])


def save_figure(fig, out, name):
    import matplotlib.pyplot as plt
    for ext in ('png', 'pdf'):
        fig.savefig(out / f'{name}.{ext}', dpi=200, bbox_inches='tight')
    plt.close(fig)


def summarize_replicates(cfg, runs):
    import numpy as np
    systems = cfg['systems']
    metrics = [k for k in runs[0]['frames'][0] if k not in ('frame', 'time_ns')]
    summary = []
    for s in systems:
        for r in system_runs(runs, s):
            row = dict(system=s, replicate=r['meta']['input']['replicate'])
            row.update({m: float(np.mean([f[m] for f in r['frames']])) for m in metrics})
            summary.append(row)
    return metrics, summary


def compare_residue_pairs(cfg, runs):
    import numpy as np
    systems = cfg['systems']
    pair_union = sorted(set().union(*(r['pairs'] for r in runs)))
    pair_rows = []
    for key in pair_union:
        for s in systems:
            values = [r['pairs'].get(key, 0.0) for r in system_runs(runs, s)]
            wt = [r['pairs'].get(key, 0.0) for r in system_runs(runs, 'WT')]
            pair_rows.append(dict(kind=key[0], Fc_label=key[1], receptor_label=key[2], system=s,
                                  mean_occupancy=float(np.mean(values)),
                                  sd_between_replicates=float(np.std(values, ddof=1)) if len(values) > 1 else '',
                                  delta_occupancy_vs_WT=float(np.mean(values) - np.mean(wt))))
    return pair_union, pair_rows


def plot_structural_summary(cfg, metrics, summary, out):
    import numpy as np
    import matplotlib.pyplot as plt
    systems = cfg['systems']
    nrow = math.ceil(len(metrics) / 3)
    fig, axes = plt.subplots(nrow, 3, figsize=(14, 3.4 * nrow), squeeze=False)
    for ax, metric in zip(axes.flat, metrics):
        for i, s in enumerate(systems):
            values = [r[metric] for r in summary if r['system'] == s]
            offsets = np.linspace(-0.12, 0.12, len(values))
            ax.scatter(i + offsets, values, color='0.3', s=20)
            ax.errorbar(i, np.mean(values), yerr=np.std(values, ddof=1) if len(values) > 1 else None,
                        color='black', marker='_', capsize=3)
        ax.set_xticks(range(len(systems)), systems, rotation=55, ha='right')
        title, ylabel = METRIC_LABELS[metric]
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(ylabel)
    for ax in list(axes.flat)[len(metrics):]:
        ax.set_visible(False)
    fig.suptitle('Structure and interface: dots = replicate means; bars = mean ± between-replicate SD')
    fig.tight_layout()
    save_figure(fig, out, 'structure_and_interface')


def plot_timecourses(cfg, runs, metrics, out):
    import numpy as np
    import matplotlib.pyplot as plt
    systems = cfg['systems']
    # Separate time-course panels preserve differences among replicates.
    for metric in metrics:
        nrow = math.ceil(len(systems) / 4)
        fig, axes = plt.subplots(nrow, 4, figsize=(14, 3 * nrow), sharex=True, sharey=True,
                                 squeeze=False)
        for ax, s in zip(axes.flat, systems):
            for r in system_runs(runs, s):
                ax.plot([f['time_ns'] for f in r['frames']], [f[metric] for f in r['frames']],
                        lw=0.8, label=f"rep{r['meta']['input']['replicate']}")
            ax.set_title(s)
            ax.set_xlabel('Stored trajectory time (ns)')
            ax.legend(fontsize=7)
        for ax in list(axes.flat)[len(systems):]:
            ax.set_visible(False)
        title, ylabel = METRIC_LABELS[metric]
        fig.supylabel(ylabel)
        fig.suptitle(title + '; common WT reference where applicable')
        fig.tight_layout()
        save_figure(fig, out, metric + '_timecourses')


def plot_contact_changes(cfg, runs, pair_union, out):
    import numpy as np
    import matplotlib.pyplot as plt
    systems = cfg['systems']
    for kind in ('contact', 'salt_bridge'):
        keys = [k for k in pair_union if k[0] == kind and max(
            np.mean([r['pairs'].get(k, 0.0) for r in system_runs(runs, s)]) for s in systems
        ) >= cfg['display_min_mean_occupancy']]
        keys.sort(key=lambda k: max(abs(np.mean([r['pairs'].get(k, 0.0) for r in system_runs(runs, s)]) -
                                       np.mean([r['pairs'].get(k, 0.0) for r in system_runs(runs, 'WT')]))
                                   for s in systems), reverse=True)
        # Paginate instead of silently dropping lower-ranked interface pairs.
        for page, start in enumerate(range(0, len(keys), 35), 1):
            selected = keys[start:start + 35]
            delta = np.array([[np.mean([r['pairs'].get(k, 0.0) for r in system_runs(runs, s)]) -
                               np.mean([r['pairs'].get(k, 0.0) for r in system_runs(runs, 'WT')])
                               for s in systems if s != 'WT'] for k in selected])
            fig, ax = plt.subplots(figsize=(9, max(4, len(selected) * 0.25)))
            im = ax.imshow(delta * 100, cmap='RdBu_r', vmin=-100, vmax=100, aspect='auto')
            ax.set_yticks(range(len(selected)), [f'{k[1]} — {k[2]}' for k in selected])
            ax.set_xticks(range(len(systems) - 1), [s for s in systems if s != 'WT'], rotation=45)
            ax.set_title(f'{kind.replace("_", " ")} occupancy changes; page {page}')
            fig.colorbar(im, ax=ax, orientation='horizontal', pad=0.3, shrink=0.75,
                         label='Variant − WT occupancy (percentage points)')
            fig.tight_layout()
            save_figure(fig, out, f'{kind}_changes_{page}')


def run_comparison(paths, out):
    """Load completed runs, summarize replicates, and render structural comparisons."""
    out = Path(out)
    if out.exists():
        raise FileExistsError(f'refusing existing output directory: {out}')
    cfg, runs = load_runs(paths)
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    out.mkdir(parents=True, exist_ok=False)
    metrics, summary = summarize_replicates(cfg, runs)
    write_csv(out / 'replicate_means.csv', ['system', 'replicate'] + metrics, summary)
    plot_structural_summary(cfg, metrics, summary, out)
    plot_timecourses(cfg, runs, metrics, out)
    pair_union, pair_rows = compare_residue_pairs(cfg, runs)
    fields = ['kind', 'Fc_label', 'receptor_label', 'system', 'mean_occupancy',
              'sd_between_replicates', 'delta_occupancy_vs_WT']
    write_csv(out / 'pair_comparisons.csv', fields, pair_rows)
    plot_contact_changes(cfg, runs, pair_union, out)
    with (out / 'figure_manifest.json').open('x') as fh:
        json.dump(dict(protocol=cfg, sources=[r['path'] for r in runs],
                       uncertainty='between-replicate SD; descriptive, not confidence intervals',
                       interface='union of residue pairs; display threshold on maximum system mean'), fh, indent=2)
    print(f'Wrote tables and PNG/PDF figures -> {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    run_comparison(args.runs, args.out)


if __name__ == '__main__':
    main()
