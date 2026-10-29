#!/usr/bin/env python3
"""One-pass contact sensitivity plus window checks on existing structure/energies.

No simulations or new energy calculations. Paths and reference identity are
validated against completed primary runs. Outputs always require a new directory.
"""
import argparse
from collections import Counter
import csv
import itertools
import json
from pathlib import Path

from defined_analysis import digest, write_csv
from plot_binding import load_energy
from plot_defined import load_runs, METRIC_LABELS


STARTS = (301, 501, 701)
STRIDES = (1, 5)
CUTOFFS = {'contact': (4.0, 4.5, 5.0), 'salt_bridge': (3.5, 4.0, 4.5)}
PRIMARY_CUTOFF = {'contact': 4.5, 'salt_bridge': 4.0}


def conditions():
    return [(kind, start, stride, cutoff) for kind, cutoffs in CUTOFFS.items()
            for start, stride, cutoff in itertools.product(STARTS, STRIDES, cutoffs)]


def frame_selected(frame, start, stride, end=1000):
    return start <= frame <= end and (frame - start) % stride == 0


def minimum_pairs(atom_pairs, distances, left_residues, right_residues, nres):
    """Encode residue pairs and reduce multiple atom distances to their minimum."""
    import numpy as np
    if not len(atom_pairs):
        return np.array([], dtype=int), np.array([], dtype=float)
    encoded = left_residues[atom_pairs[:, 0]] * nres + right_residues[atom_pairs[:, 1]]
    keys, inverse = np.unique(encoded, return_inverse=True)
    minimum = np.full(len(keys), np.inf)
    np.minimum.at(minimum, inverse, distances)
    return keys, minimum


def analyze_run(run, out):
    import numpy as np
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance
    if out.exists():
        raise FileExistsError(out)
    meta = run['meta']
    source = meta['input']
    if digest(source['topology']) != meta['topology_sha256']:
        raise ValueError('topology changed since primary analysis')
    u = mda.Universe(source['topology'], source['trajectory'])
    with (Path(run['path']) / 'residue_map.csv').open() as fh:
        mapping = list(csv.DictReader(fh))
    if len(u.residues) != len(mapping) or len(u.trajectory) < 1000:
        raise ValueError('unexpected residue count or incomplete trajectory')
    if any(str(r.resname) != row['name'] for r, row in zip(u.residues, mapping)):
        raise ValueError('topology identities differ from primary residue map')
    nres = len(mapping)
    heavy = u.atoms.select_atoms('not name H*')
    fc = heavy[np.array([mapping[i]['partner'] == 'Fc' for i in heavy.resindices])]
    rc = heavy[np.array([mapping[i]['partner'] == 'FcGR3a' for i in heavy.resindices])]
    positive = '(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ) or (resname HIP and name ND1 NE2)'
    negative = '(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)'
    groups = {'contact': [(fc, rc)],
              'salt_bridge': [(fc.select_atoms(positive), rc.select_atoms(negative)),
                              (fc.select_atoms(negative), rc.select_atoms(positive))]}
    counters = {c: Counter() for c in conditions()}
    nframes = Counter()
    clock = []
    for ts in u.trajectory[300:1000]:
        frame = ts.frame + 1
        if not np.all(np.isfinite(u.atoms.positions)) or not np.isfinite(ts.time):
            raise ValueError(f'nonfinite coordinates/time at frame {frame}')
        if clock and ts.time <= clock[-1]:
            raise ValueError('stored trajectory times are not increasing')
        clock.append(float(ts.time))
        for kind, cutoffs in CUTOFFS.items():
            hits = {cutoff: set() for cutoff in cutoffs}
            for left, right in groups[kind]:
                if not len(left) or not len(right):
                    continue
                atoms, distance = capped_distance(left.positions, right.positions,
                                                 max(cutoffs), return_distances=True)
                keys, minimum = minimum_pairs(atoms, distance, left.resindices, right.resindices, nres)
                for cutoff in cutoffs:
                    hits[cutoff].update(map(int, keys[minimum <= cutoff]))
            for start, stride in itertools.product(STARTS, STRIDES):
                if frame_selected(frame, start, stride):
                    for cutoff in cutoffs:
                        condition = (kind, start, stride, cutoff)
                        counters[condition].update(hits[cutoff])
                        nframes[condition] += 1
    # Exact regression against primary contact counts, not just matching means.
    index = {row['label']: i for i, row in enumerate(mapping)}
    expected = {kind: Counter() for kind in CUTOFFS}
    with (Path(run['path']) / 'pairs.csv').open() as fh:
        for row in csv.DictReader(fh):
            key = index[row['Fc_label']] * nres + index[row['receptor_label']]
            expected[row['kind']][key] = int(row['observed_frames'])
    for kind in CUTOFFS:
        if counters[(kind, 301, 5, PRIMARY_CUTOFF[kind])] != expected[kind]:
            raise ValueError(f'sensitivity baseline differs from primary: {source}, {kind}')
    rows, counts = [], []
    for condition, counter in counters.items():
        kind, start, stride, cutoff = condition
        count = nframes[condition]
        info = dict(system=source['system'], replicate=source['replicate'], kind=kind,
                    start_frame=start, stride=stride, cutoff_A=cutoff, analyzed_frames=count)
        counts.append(dict(info, mean_pairs=sum(counter.values()) / count))
        for key, observed in sorted(counter.items()):
            i, j = divmod(key, nres)
            rows.append(dict(info, Fc_label=mapping[i]['label'], receptor_label=mapping[j]['label'],
                             occupancy=observed / count))
    out.mkdir(parents=True, exist_ok=False)
    fields = ['system', 'replicate', 'kind', 'start_frame', 'stride', 'cutoff_A', 'analyzed_frames']
    write_csv(out / 'counts.csv', fields + ['mean_pairs'], counts)
    write_csv(out / 'pairs.csv', fields + ['Fc_label', 'receptor_label', 'occupancy'], rows)
    with (out / 'complete.json').open('x') as fh:
        json.dump(dict(primary_run=run['path'], input=source, protocol_sha256=meta['protocol_sha256'],
                       script_sha256=digest(__file__), baseline_counts_match=True,
                       starts=STARTS, strides=STRIDES, cutoffs=CUTOFFS,
                       dense_frames=len(clock), first_time_ns=clock[0] / 1000,
                       last_time_ns=clock[-1] / 1000, MDAnalysis=mda.__version__, numpy=np.__version__), fh, indent=2)
    print(f"Sensitivity completed {source['system']} rep{source['replicate']}", flush=True)


def summarize(cfg, runs, root):
    import numpy as np
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'pdf.fonttype': 42})
    folder = root / 'summary'
    folder.mkdir(exist_ok=False)
    systems = cfg['systems']
    variants = [s for s in systems if s != 'WT']
    pairs, count_rows = {}, []
    for run in runs:
        source = run['meta']['input']
        path = root / source['system'] / f"rep{source['replicate']}"
        if not (path / 'complete.json').is_file():
            raise ValueError(f'incomplete sensitivity run: {path}')
        with (path / 'counts.csv').open() as fh:
            for r in csv.DictReader(fh):
                count_rows.append(dict(system=r['system'], replicate=int(r['replicate']), kind=r['kind'],
                                       start_frame=int(r['start_frame']), stride=int(r['stride']),
                                       cutoff_A=float(r['cutoff_A']), analyzed_frames=int(r['analyzed_frames']),
                                       mean_pairs=float(r['mean_pairs'])))
        with (path / 'pairs.csv').open() as fh:
            for r in csv.DictReader(fh):
                key = (r['system'], int(r['replicate']), r['kind'], int(r['start_frame']),
                       int(r['stride']), float(r['cutoff_A']), r['Fc_label'], r['receptor_label'])
                pairs[key] = float(r['occupancy'])
    write_csv(folder / 'replicate_counts.csv', list(count_rows[0]), count_rows)
    def occupancy(s, condition, pair):
        return float(np.mean([pairs.get((s, rep, *condition, *pair), 0.) for rep in cfg['replicates']]))
    def delta(s, condition, pair):
        return occupancy(s, condition, pair) - occupancy('WT', condition, pair)
    def count(s, kind, start, stride, cutoff):
        return float(np.mean([r['mean_pairs'] for r in count_rows if
                              (r['system'], r['kind'], r['start_frame'], r['stride'], r['cutoff_A']) ==
                              (s, kind, start, stride, cutoff)]))
    union = {kind: sorted({(k[-2], k[-1]) for k in pairs if k[2] == kind}) for kind in CUTOFFS}
    pair_rows = []
    for kind in CUTOFFS:
        for pair in union[kind]:
            for s in systems:
                for condition in [c for c in conditions() if c[0] == kind]:
                    pair_rows.append(dict(system=s, kind=kind, Fc_label=pair[0], receptor_label=pair[1],
                                          start_frame=condition[1], stride=condition[2], cutoff_A=condition[3],
                                          mean_occupancy=occupancy(s, condition, pair),
                                          delta_vs_WT=delta(s, condition, pair)))
    write_csv(folder / 'pair_sensitivity.csv', list(pair_rows[0]), pair_rows)
    colors = {s: plt.get_cmap('tab10')(i) for i, s in enumerate(variants)}
    def save(fig, name):
        for ext in ('png', 'pdf'):
            fig.savefig(folder / f'{name}.{ext}', dpi=180, bbox_inches='tight')
        plt.close(fig)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    for i, kind in enumerate(CUTOFFS):
        for s in variants:
            axes[i, 0].plot(CUTOFFS[kind], [count(s, kind, 301, 5, c) - count('WT', kind, 301, 5, c)
                                          for c in CUTOFFS[kind]], marker='o', color=colors[s], label=s)
            axes[i, 1].plot(STARTS, [count(s, kind, start, 5, PRIMARY_CUTOFF[kind]) -
                                   count('WT', kind, start, 5, PRIMARY_CUTOFF[kind]) for start in STARTS],
                            marker='o', color=colors[s])
        for ax in axes[i]:
            ax.axhline(0, color='0.4', ls='--', lw=0.8)
            ax.set_ylabel('Variant − WT mean residue-pair count')
        axes[i, 0].set_title(kind.replace('_', ' ') + ': cutoff sensitivity')
        axes[i, 0].set_xlabel('Distance cutoff (Å)')
        axes[i, 0].set_xticks(CUTOFFS[kind])
        axes[i, 1].set_title(kind.replace('_', ' ') + ': window sensitivity')
        axes[i, 1].set_xlabel('First trajectory frame (all end at frame 1000)')
        axes[i, 1].set_xticks(STARTS)
    axes[0, 0].legend(ncol=4, fontsize=8)
    fig.suptitle('Sensitivity of mean interface counts; 50 ps sampling, three replicate means per system')
    fig.tight_layout()
    save(fig, 'contact_count_sensitivity')
    sampling = []
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, kind in zip(axes, CUTOFFS):
        all_errors = []
        for s in variants:
            coarse = np.array([delta(s, (kind, 301, 5, PRIMARY_CUTOFF[kind]), p) for p in union[kind]])
            dense = np.array([delta(s, (kind, 301, 1, PRIMARY_CUTOFF[kind]), p) for p in union[kind]])
            ax.scatter(coarse * 100, dense * 100, s=9, alpha=0.5, color=colors[s], label=s)
            error = (dense - coarse) * 100
            all_errors.extend(error)
            sampling.append(dict(kind=kind, system=s, max_abs_difference_pp=float(np.max(np.abs(error))),
                                 rms_difference_pp=float(np.sqrt(np.mean(error ** 2)))))
        ax.plot([-100, 100], [-100, 100], ls='--', color='0.4')
        ax.set_title(kind.replace('_', ' ') + f'; max deviation {max(map(abs, all_errors)):.1f} pp')
        ax.set_xlabel('Variant − WT occupancy at 50 ps (%)')
        ax.set_ylabel('Variant − WT occupancy at 10 ps (%)')
        ax.set_aspect('equal')
    axes[1].legend(fontsize=7, loc='upper left')
    fig.tight_layout()
    save(fig, 'sampling_sensitivity')
    write_csv(folder / 'sampling_comparison.csv', list(sampling[0]), sampling)
    findings = []
    for kind, number in [('contact', 12), ('salt_bridge', 3)]:
        candidates = [(abs(delta(s, (kind, 301, 5, PRIMARY_CUTOFF[kind]), p)), s, p)
                      for s in variants for p in union[kind]]
        for _, s, pair in sorted(candidates, reverse=True)[:number]:
            values = [delta(s, c, pair) for c in conditions() if c[0] == kind]
            findings.append(dict(kind=kind, system=s, Fc_label=pair[0], receptor_label=pair[1],
                                 primary_delta_pp=100 * delta(s, (kind, 301, 5, PRIMARY_CUTOFF[kind]), pair),
                                 minimum_delta_pp=100 * min(values), maximum_delta_pp=100 * max(values),
                                 direction_retained=min(values) >= 0 or max(values) <= 0))
    write_csv(folder / 'finding_ranges.csv', list(findings[0]), findings)
    fig, ax = plt.subplots(figsize=(13, 7))
    matrix = [[100 * delta(r['system'], (r['kind'], start, stride, cutoff),
                           (r['Fc_label'], r['receptor_label']))
               for start, stride, cutoff_index in itertools.product(STARTS, STRIDES, range(3))
               for cutoff in [CUTOFFS[r['kind']][cutoff_index]]] for r in findings]
    im = ax.imshow(matrix, vmin=-100, vmax=100, cmap='RdBu_r', aspect='auto')
    labels = [f"{r['system']} {r['kind']}: {r['Fc_label']}–{r['receptor_label']}" for r in findings]
    ax.set_yticks(range(len(labels)), labels, fontsize=8)
    ax.set_xticks(range(18), [f'{start}/{stride}/{ci + 1}' for start, stride, ci in
                            itertools.product(STARTS, STRIDES, range(3))], rotation=60)
    ax.set_xlabel('Start frame / stride / cutoff index (contact: 4, 4.5, 5 Å; salt: 3.5, 4, 4.5 Å)')
    ax.set_title('Largest primary occupancy changes across 18 analysis choices per row')
    fig.colorbar(im, ax=ax, label='Variant − WT occupancy (percentage points)')
    fig.tight_layout()
    save(fig, 'finding_robustness')
    structural, binding = [], []
    for run in runs:
        source = run['meta']['input']
        energy = load_energy(source['mmgbsa_csv'])
        if len(energy) != len(run['frames']):
            raise ValueError('energy/structure frame counts differ')
        for start in STARTS:
            chosen = [i for i, f in enumerate(run['frames']) if f['frame'] >= start]
            common = dict(system=source['system'], replicate=source['replicate'], start_frame=start, n_frames=len(chosen))
            structural.append(dict(common, **{k: float(np.mean([run['frames'][i][k] for i in chosen]))
                                              for k in run['frames'][0] if k not in ('frame', 'time_ns')}))
            binding.append(dict(common, energy_mean=float(np.mean([energy[i]['DELTA TOTAL'] for i in chosen]))))
    write_csv(folder / 'structural_windows.csv', list(structural[0]), structural)
    write_csv(folder / 'binding_windows.csv', list(binding[0]), binding)
    def window_mean(table, s, start, key):
        return float(np.mean([r[key] for r in table if r['system'] == s and r['start_frame'] == start]))
    binding_summary = []
    fig, ax = plt.subplots(figsize=(10, 5))
    for j, start in enumerate(STARTS):
        differences = [window_mean(binding, s, start, 'energy_mean') -
                       window_mean(binding, 'WT', start, 'energy_mean') for s in variants]
        ax.plot(range(len(variants)), differences, marker='o', label=f'First frame {start}')
        for s, value in zip(variants, differences):
            binding_summary.append(dict(system=s, start_frame=start, ddG_vs_WT=value))
    ax.set_xticks(range(len(variants)), variants)
    ax.axhline(0, color='0.4', ls='--')
    ax.set_ylabel('Whole-partner MM/GBSA variant − WT (kcal/mol)')
    ax.set_title('Window sensitivity of mean energy differences; existing 50 ps snapshots')
    ax.legend()
    fig.tight_layout()
    save(fig, 'binding_window_sensitivity')
    write_csv(folder / 'binding_window_comparisons.csv', list(binding_summary[0]), binding_summary)
    keys = ['CH2_common_reference_rmsd_A', 'receptor_common_reference_rmsd_A',
            'Fc_self_fit_rmsd_A', 'receptor_self_fit_rmsd_A', 'partner_COM_distance_A', 'complex_Rg_A']
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    for ax, key in zip(axes.flat, keys):
        for start in STARTS:
            differences = [window_mean(structural, s, start, key) - window_mean(structural, 'WT', start, key)
                           for s in variants]
            ax.plot(range(len(variants)), differences, marker='o', label=f'First frame {start}')
        ax.axhline(0, color='0.4', ls='--', lw=0.8)
        title, ylabel = METRIC_LABELS[key]
        ax.set_title(title)
        ax.set_ylabel('Variant − WT ' + ylabel)
        ax.set_xticks(range(len(variants)), variants, rotation=45)
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('Window sensitivity of structural means; common WT reference, 50 ps snapshots')
    fig.tight_layout()
    save(fig, 'structure_window_sensitivity')
    with (folder / 'complete.json').open('x') as fh:
        json.dump(dict(protocol=cfg, starts=STARTS, strides=STRIDES, cutoffs=CUTOFFS,
                       analyzed_runs=len(runs), baselines_match=True,
                       limitations='Descriptive sensitivity checks on nested windows of the same trajectories; '
                                   'not independent replicates, equivalence tests, or proof of convergence.'), fh, indent=2)
    print(f'Sensitivity summary -> {folder}', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('runs', nargs='+', type=Path)
    parser.add_argument('--out', required=True, type=Path)
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError(args.out)
    cfg, runs = load_runs(args.runs)
    if (cfg['start_frame'], cfg['end_frame'], cfg['stride']) != (301, 1000, 5):
        raise ValueError('this sensitivity grid requires the approved primary 301/1000/5 protocol')
    args.out.mkdir(parents=True, exist_ok=False)
    for run in runs:
        source = run['meta']['input']
        analyze_run(run, args.out / source['system'] / f"rep{source['replicate']}")
    summarize(cfg, runs, args.out)


if __name__ == '__main__':
    main()
