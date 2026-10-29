#!/usr/bin/env python3
"""Explicit input inventory and fixed-window structure/contact analysis.

Inventory uses only the standard library. The run command requires numpy and
MDAnalysis and reads a single trajectory. Follow local analysis compute policies.
"""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform

from ..io import digest, write_csv


def read_protocol(path):
    path = Path(path).resolve()
    cfg = json.loads(path.read_text())
    if 'WT' not in cfg['systems'] or len(set(cfg['systems'])) != len(cfg['systems']):
        raise ValueError('systems must contain WT exactly once and have no duplicates')
    if not cfg['replicates'] or any(type(r) is not int or r < 1 for r in cfg['replicates']):
        raise ValueError('replicates must be positive integers')
    if len(set(cfg['replicates'])) != len(cfg['replicates']):
        raise ValueError('replicates must be unique')
    for key in ('start_frame', 'end_frame', 'stride'):
        if type(cfg[key]) is not int or cfg[key] < 1:
            raise ValueError(f'{key} must be a positive integer')
    if cfg['end_frame'] < cfg['start_frame']:
        raise ValueError('end_frame precedes start_frame')
    for key in ('heavy_atom_contact_A', 'salt_bridge_A'):
        if cfg[key] <= 0:
            raise ValueError(f'{key} must be positive')
    if not 0 <= cfg['display_min_mean_occupancy'] <= 1:
        raise ValueError('display occupancy must be between zero and one')
    return path, cfg


def reference_residues(path):
    """Read protein backbone and original numbering without arithmetic offsets."""
    residues = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM  ') or line[21] not in 'ABC':
                continue
            if line[16] not in (' ', 'A'):
                continue
            key = (line[21], int(line[22:26]), line[26].strip())
            r = residues.setdefault(key, {'chain': key[0], 'number': key[1],
                                         'icode': key[2], 'name': line[17:20].strip(),
                                         'backbone': {}})
            atom = line[12:16].strip()
            if atom in ('N', 'CA', 'C'):
                r['backbone'][atom] = [float(line[i:i + 8]) for i in (30, 38, 46)]
    return [r for chain in 'ABC' for key, r in residues.items() if key[0] == chain]


def residue_map(reference, cfg):
    rows = reference_residues(reference)
    for chain, expected in cfg['expected_protein_residue_counts'].items():
        if sum(r['chain'] == chain for r in rows) != expected:
            raise ValueError(f'reference chain {chain}: unexpected residue count')
    for i, row in enumerate(rows, 1):
        if set(row['backbone']) != {'N', 'CA', 'C'}:
            raise ValueError(f'incomplete backbone: {row}')
        row.update(index=i, label=f"{row['chain']}:{row['number']}{row['icode']}",
                   partner='Fc' if row['chain'] in 'AB' else 'FcGR3a',
                   component=f"protein_{row['chain']}")
    for k, name in enumerate(cfg['expected_glycan_names']):
        if k < 16:
            chain, site, ring = ('A' if k < 8 else 'B'), 297, k % 8 + 1
        else:
            chain, site, ring = 'C', (162 if k == 16 else 45), 1
        rows.append(dict(index=len(rows) + 1, label=f'{chain}:N{site}:glycan{ring}',
                         name=name, chain=chain, number=site, icode='',
                         partner='Fc' if chain in 'AB' else 'FcGR3a',
                         component=f'glycan_{chain}_N{site}', backbone={}))
    return rows


def input_records(protocol, cfg):
    root = (protocol.parent / cfg['prep_root']).resolve()
    records = []
    trajectories = set()
    for system in cfg['systems']:
        for rep in cfg['replicates']:
            folder = root / system / f'rep{rep}'
            record = dict(system=system, replicate=rep)
            for key, name in [('topology', 'dry.prmtop'), ('trajectory', 'prod_dry.nc'),
                              ('mmgbsa_csv', 'perframe.csv')]:
                p = (folder / name).resolve()
                record[key] = str(p)
                record[key + '_exists'] = p.is_file()
            if record['trajectory'] in trajectories:
                raise ValueError('duplicate resolved trajectory: ' + record['trajectory'])
            trajectories.add(record['trajectory'])
            records.append(record)
    return records


def pair_counts(atom_pairs, left_residues, right_residues):
    """One residue-pair observation per frame, regardless of atom-pair count."""
    return {(int(left_residues[i]), int(right_residues[j])) for i, j in atom_pairs}


def kabsch_fit(x, y):
    """Return a proper rotation and centers for row-vector coordinates."""
    import numpy as np
    xc, yc = x.mean(axis=0), y.mean(axis=0)
    left, _, right = np.linalg.svd((x - xc).T @ (y - yc))
    flip = np.eye(3)
    flip[-1, -1] = np.linalg.det(left @ right)
    return left @ flip @ right, xc, yc


def load_and_validate(protocol, cfg, system, replicate):
    """Locate inputs and verify topology, mutation identity, masses, and window."""
    import numpy as np
    import MDAnalysis as mda

    records = input_records(protocol, cfg)
    rec = next((r for r in records if r['system'] == system and
                r['replicate'] == replicate), None)
    if rec is None:
        raise ValueError('system/replicate is not in protocol')
    if not rec['topology_exists'] or not rec['trajectory_exists']:
        raise ValueError('missing topology or trajectory')
    reference = (protocol.parent / cfg['reference_pdb']).resolve()
    mapping = residue_map(reference, cfg)
    universe = mda.Universe(rec['topology'], rec['trajectory'])
    if len(universe.residues) != len(mapping):
        raise ValueError('wrong residue count; use original dry topology, not MM reordered topology')
    if not np.all(np.isfinite(universe.atoms.masses)) or np.any(universe.atoms.masses <= 0):
        raise ValueError('positive finite atom masses are required for COM and Rg')
    mutation = None if system == 'WT' else (int(system[1:-1]), system[-1])
    aa3 = dict(A='ALA', S='SER')
    aliases = {'NLN': 'ASN', 'CYX': 'CYS', 'HID': 'HIS', 'HIE': 'HIS', 'HIP': 'HIS'}
    for residue, row in zip(universe.residues, mapping):
        expected = row['name']
        if mutation and row['component'].startswith('protein_') and row['chain'] in 'AB' and row['number'] == mutation[0]:
            expected = aa3[mutation[1]]
        actual = str(residue.resname)
        if row['component'].startswith('protein_'):
            actual, expected = aliases.get(actual, actual), aliases.get(expected, expected)
        if actual != expected:
            raise ValueError(f"residue layout mismatch at {row['index']}: {actual} != {expected}")
        row['name'] = str(residue.resname)
    if len(universe.trajectory) < cfg['end_frame']:
        raise ValueError('trajectory shorter than specified window; refusing partial analysis')
    return universe, mapping, rec, reference


def prepare_selections(universe, mapping, cfg):
    """Build fixed atom groups and backbone/reference indices once per run."""
    import numpy as np

    backbone, target = [], []
    for residue, row in zip(universe.residues[:596], mapping[:596]):
        for name in ('N', 'CA', 'C'):
            matches = residue.atoms.select_atoms('name ' + name)
            if len(matches) != 1:
                raise ValueError(f"missing/ambiguous {name} at {row['label']}")
            backbone.append(matches[0].index)
            target.append(row['backbone'][name])
    target = np.asarray(target)
    backbone = np.asarray(backbone)
    def group(ranges):
        return np.array([3 * (i - 1) + j for lo, hi in ranges
                         for i in range(lo, hi + 1) for j in range(3)])
    ch3 = group(cfg['ch3_resids'])
    ch2 = group(cfg['ch2_resids'])
    receptor = group([[426, 596]])
    fc_bb = group([[1, 425]])
    heavy = universe.atoms.select_atoms('not name H*')
    fc_heavy = heavy[np.array([mapping[i]['partner'] == 'Fc' for i in heavy.resindices])]
    receptor_heavy = heavy[np.array([mapping[i]['partner'] == 'FcGR3a' for i in heavy.resindices])]
    charged = '(resname ARG and name NE NH1 NH2) or (resname LYS and name NZ) or (resname HIP and name ND1 NE2)'
    acidic = '(resname ASP and name OD1 OD2) or (resname GLU and name OE1 OE2)'
    selections = [(fc_heavy, receptor_heavy)]
    salt_selections = [(fc_heavy.select_atoms(charged), receptor_heavy.select_atoms(acidic)),
                       (fc_heavy.select_atoms(acidic), receptor_heavy.select_atoms(charged))]
    fc_all = universe.atoms[np.array([mapping[i]['partner'] == 'Fc' for i in universe.atoms.resindices])]
    receptor_all = universe.atoms[np.array([mapping[i]['partner'] == 'FcGR3a' for i in universe.atoms.resindices])]
    return dict(backbone=backbone, target=target, ch3=ch3,
                domains=[('CH3', ch3), ('CH2', ch2), ('receptor', receptor)],
                self_fit_domains=[('Fc', fc_bb), ('receptor', receptor)],
                fc_all=fc_all, receptor_all=receptor_all,
                contacts=[('contact', selections, cfg['heavy_atom_contact_A']),
                          ('salt_bridge', salt_selections, cfg['salt_bridge_A'])])


def rms_distance(x, y):
    """Root-mean-square Cartesian displacement between matching coordinates."""
    import numpy as np
    return float(np.sqrt(np.mean(np.sum((x - y) ** 2, axis=1))))


def validate_frame(universe, timestep, previous_time):
    """Reject invalid coordinates and non-increasing sampled times."""
    import numpy as np
    if not np.all(np.isfinite(universe.atoms.positions)) or not np.isfinite(timestep.time):
        raise ValueError(f'nonfinite coordinates or time at frame {timestep.frame + 1}')
    if previous_time is not None and timestep.time / 1000 <= previous_time:
        raise ValueError('trajectory time must increase between sampled frames')


def calculate_structural_metrics(universe, selections):
    """Calculate common-reference/self-fit RMSDs, partner distance, and size."""
    import numpy as np

    backbone, target, ch3 = (selections[key] for key in ('backbone', 'target', 'ch3'))
    x = universe.atoms.positions[backbone].astype(float)
    rot, center, refcenter = kabsch_fit(x[ch3], target[ch3])
    aligned = (x - center) @ rot + refcenter
    metrics = {}
    for name, ids in selections['domains']:
        metrics[name + '_common_reference_rmsd_A'] = rms_distance(aligned[ids], target[ids])
    for name, ids in selections['self_fit_domains']:
        rot2, cx, cy = kabsch_fit(x[ids], target[ids])
        metrics[name + '_self_fit_rmsd_A'] = rms_distance((x[ids] - cx) @ rot2 + cy, target[ids])
    fc_center = selections['fc_all'].center_of_mass()
    receptor_center = selections['receptor_all'].center_of_mass()
    metrics['partner_COM_distance_A'] = float(np.linalg.norm(fc_center - receptor_center))
    metrics['complex_Rg_A'] = float(universe.atoms.radius_of_gyration())
    return metrics


def calculate_contact_pairs(contact_groups):
    """Return unique residue pairs per contact type for the current frame."""
    from MDAnalysis.lib.distances import capped_distance

    observations = {}
    for kind, groups, cutoff in contact_groups:
        pairs = set()
        for left, right in groups:
            if len(left) and len(right):
                # Deliberately no minimum-image distances: use the imaged complex.
                atom_pairs = capped_distance(left.positions, right.positions, cutoff,
                                             return_distances=False)
                pairs.update(pair_counts(atom_pairs, left.resindices, right.resindices))
        observations[kind] = pairs
    return observations


def analyze_frames(universe, selections, cfg):
    """Analyze the configured frame window in one pass through the trajectory."""
    counts = {'contact': Counter(), 'salt_bridge': Counter()}
    frames = []
    window = slice(cfg['start_frame'] - 1, cfg['end_frame'], cfg['stride'])
    for timestep in universe.trajectory[window]:
        previous_time = frames[-1]['time_ns'] if frames else None
        validate_frame(universe, timestep, previous_time)
        metrics = calculate_structural_metrics(universe, selections)
        contacts = calculate_contact_pairs(selections['contacts'])
        for kind, pairs in contacts.items():
            counts[kind].update(pairs)
        observations = {kind + '_residue_pairs': len(pairs)
                        for kind, pairs in contacts.items()}
        frames.append(dict(frame=timestep.frame + 1, time_ns=float(timestep.time) / 1000,
                           **metrics, **observations))
    return frames, counts


def write_results(out, frames, counts, mapping, protocol, cfg, rec, reference):
    """Write the existing CSV schema and mark completion only after all outputs."""
    import numpy as np
    import MDAnalysis as mda

    out.mkdir(parents=True, exist_ok=False)
    write_csv(out / 'frames.csv', list(frames[0]), frames)
    map_fields = ['index', 'label', 'name', 'partner', 'component', 'chain', 'number', 'icode']
    write_csv(out / 'residue_map.csv', map_fields, mapping)
    fields = ['kind', 'Fc_label', 'Fc_component', 'receptor_label', 'receptor_component',
              'observed_frames', 'analyzed_frames', 'occupancy']
    rows = []
    for kind, counter in counts.items():
        for (i, j), n in sorted(counter.items()):
            rows.append(dict(kind=kind, Fc_label=mapping[i]['label'], Fc_component=mapping[i]['component'],
                             receptor_label=mapping[j]['label'], receptor_component=mapping[j]['component'],
                             observed_frames=n, analyzed_frames=len(frames), occupancy=n / len(frames)))
    write_csv(out / 'pairs.csv', fields, rows)
    provenance = dict(protocol=cfg, protocol_sha256=digest(protocol), input=rec,
                      reference_sha256=digest(reference), topology_sha256=digest(rec['topology']),
                      python=platform.python_version(), numpy=np.__version__, MDAnalysis=mda.__version__,
                      n_frames=len(frames), first_time_ns=frames[0]['time_ns'], last_time_ns=frames[-1]['time_ns'])
    # Written last: downstream tools accept only completed runs.
    with (out / 'complete.json').open('x') as fh:
        json.dump(provenance, fh, indent=2)


def run(protocol, cfg, system, replicate, out):
    """Validate inputs, prepare selections, analyze frames, and write results."""
    out = Path(out)
    if out.exists():
        raise FileExistsError(f'refusing to replace existing output directory: {out}')
    universe, mapping, record, reference = load_and_validate(protocol, cfg, system, replicate)
    selections = prepare_selections(universe, mapping, cfg)
    frames, counts = analyze_frames(universe, selections, cfg)
    write_results(out, frames, counts, mapping, protocol, cfg, record, reference)
    print(f'Completed {system} rep{replicate}: {len(frames)} frames -> {out}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--protocol', type=Path, default=Path(__file__).resolve().parents[2] / 'config/analysis_protocol.json')
    sub = parser.add_subparsers(dest='command', required=True)
    inventory = sub.add_parser('inventory', help='lightweight input audit; no trajectory reads')
    inventory.add_argument('--out', type=Path, required=True)
    analyze = sub.add_parser('run', help='trajectory analysis; follow local compute policies')
    analyze.add_argument('--system', required=True)
    analyze.add_argument('--replicate', required=True, type=int)
    analyze.add_argument('--out', type=Path, required=True)
    args = parser.parse_args()
    protocol, cfg = read_protocol(args.protocol)
    if args.command == 'inventory':
        records = input_records(protocol, cfg)
        mapping = residue_map((protocol.parent / cfg['reference_pdb']).resolve(), cfg)
        with args.out.open('x') as fh:
            json.dump(dict(protocol=cfg, inputs=records, residue_map=mapping), fh, indent=2)
        print(f'{len(records)} runs, {len(mapping)} residues; inventory -> {args.out}')
    else:
        run(protocol, cfg, args.system, args.replicate, args.out)


if __name__ == '__main__':
    main()
