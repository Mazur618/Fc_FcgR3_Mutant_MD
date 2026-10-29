"""Synthetic checks for import safety, preprocessing, statistics, and plotting stages."""
import csv
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import unittest
from types import SimpleNamespace
from unittest.mock import patch

from analysis import sensitivity
from analysis.io import write_csv
from analysis.mmgbsa import plot_binding
from analysis.structure import plot_structure
from preprocessing import build_complex
from slurm import mutate_protein
from test_analysis import SharedFilesystemTemporaryDirectory


HERE = Path(__file__).resolve().parent


def make_summary_runs(root):
    """Six tiny completed runs with known means and nonzero replicate variance."""
    cfg = dict(systems=['WT', 'D265A'], replicates=[1, 2, 3], start_frame=301,
               end_frame=1000, stride=5, display_min_mean_occupancy=.2)
    paths = []
    for system in cfg['systems']:
        for rep in cfg['replicates']:
            folder = root / system / f'rep{rep}'
            folder.mkdir(parents=True)
            energy = folder / 'energy.csv'
            rows = ['DELTA Energy Terms', 'Frame #,VDWAALS,EEL,EGB,ESURF,DELTA TOTAL']
            for ordinal in range(3):
                vdw = -10 + rep + ordinal - (10 if system == 'D265A' else 0)
                rows.append(f'{ordinal},{vdw},-20,5,1,{vdw-14}')
            energy.write_text('\n'.join(rows) + '\n\n')
            frames = [dict(frame=frame, time_ns=frame / 100,
                           **{key: float(rep + i + 1) for i, key in enumerate(plot_structure.METRIC_LABELS)})
                      for frame in (301, 501, 701)]
            write_csv(folder / 'frames.csv', list(frames[0]), frames)
            pairs = [dict(kind=kind, Fc_label='A:265', receptor_label='C:120',
                          occupancy=(.8 if system == 'WT' else .2))
                     for kind in ('contact', 'salt_bridge')]
            write_csv(folder / 'pairs.csv', list(pairs[0]), pairs)
            meta = dict(input=dict(system=system, replicate=rep, mmgbsa_csv=str(energy),
                                   trajectory=str(folder / 'synthetic.nc')),
                        protocol=cfg, protocol_sha256='synthetic protocol',
                        reference_sha256='synthetic reference', n_frames=3)
            (folder / 'complete.json').write_text(json.dumps(meta))
            paths.append(folder)
    return cfg, paths


class PreprocessingStageTests(unittest.TestCase):
    def test_imports_have_no_scientific_initialization_or_file_writes(self):
        command = '''import sys
from preprocessing import build_complex, extract_protein
from slurm import mutate_protein
assert 'pyrosetta' not in sys.modules
assert 'gemmi' not in sys.modules
'''
        result = subprocess.run([sys.executable, '-B', '-c', command], cwd=HERE.parent,
                                capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_preparation_help_does_not_require_scientific_dependencies(self):
        for name in ('preprocessing/build_complex.py', 'slurm/mutate_protein.py'):
            result = subprocess.run([sys.executable, '-B', name, '--help'], cwd=HERE.parent,
                                    capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn('usage:', result.stdout)

    def test_mutation_specification_and_protected_sites(self):
        args = build_complex.parse_args(['--name', 'D265A', '--mut', 'Fc:265:ALA'])
        mutations, histidines = build_complex.validate_sites(args)
        self.assertEqual(mutations, [('A', 265, 'ALA'), ('B', 265, 'ALA')])
        self.assertEqual(histidines, {})
        for site in ('A:261:ALA', 'A:297:ALA'):
            args.mut = [site]
            with self.assertRaises(SystemExit):
                build_complex.validate_sites(args)

    def test_extraction_keeps_the_original_cleanup_order(self):
        from unittest.mock import Mock, call
        from preprocessing.extract_protein import extract_protein
        structure = Mock()
        gemmi = SimpleNamespace(read_structure=Mock(return_value=structure))
        with patch.dict(sys.modules, gemmi=gemmi):
            extract_protein('source.cif', 'protein.pdb')
        gemmi.read_structure.assert_called_once_with('source.cif')
        self.assertEqual(structure.mock_calls, [call.remove_hydrogens(),
                         call.remove_alternative_conformations(),
                         call.remove_ligands_and_waters(), call.write_pdb('protein.pdb')])

    def test_standalone_mutation_resolves_pdb_numbering(self):
        pose = SimpleNamespace(pdb_info=lambda: SimpleNamespace(pdb2pose=lambda chain, num: 579))
        args = mutate_protein.parse_args(['reference.pdb', '--chain', 'C', '--num', '158', '--aa', 'VAL'])
        self.assertEqual(mutate_protein.resolve_residue(pose, args), 579)


@unittest.skipUnless(all(importlib.util.find_spec(name) for name in ('numpy', 'scipy', 'matplotlib')),
                     'NumPy, SciPy, and Matplotlib are required for comparison tests')
class ComparisonStageTests(unittest.TestCase):
    def setUp(self):
        self.temp = SharedFilesystemTemporaryDirectory(dir=HERE)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        cache = patch.dict(os.environ, MPLCONFIGDIR=str(self.root / 'matplotlib-cache'))
        cache.start()
        self.addCleanup(cache.stop)
        import matplotlib
        matplotlib.use('Agg')
        self.cfg, self.paths = make_summary_runs(self.root)
        _, self.runs = plot_structure.load_runs(self.paths)

    def test_binding_means_differences_and_variability(self):
        rows, sources = plot_binding.collect_replicate_energies(self.runs)
        stats = plot_binding.compare_binding_energies(self.cfg['systems'], rows)
        self.assertEqual(len(sources), 6)
        wt, mutant = stats
        self.assertAlmostEqual(wt['mean_kcal_mol'], -21.)
        self.assertAlmostEqual(mutant['mean_kcal_mol'], -31.)
        self.assertAlmostEqual(mutant['ddG_vs_WT'], -10.)
        self.assertAlmostEqual(mutant['sd_between_replicates'], 1.)
        self.assertLess(mutant['ci95_lo'], -10.)
        self.assertGreater(mutant['ci95_hi'], -10.)
        self.assertAlmostEqual(mutant['p_holm'], mutant['p_welch'])

    def test_missing_contact_uses_zero_without_omitting_replicate(self):
        self.runs[-1]['pairs'] = {}
        _, rows = plot_structure.compare_residue_pairs(self.cfg, self.runs)
        row = next(r for r in rows if r['system'] == 'D265A' and r['kind'] == 'contact')
        self.assertAlmostEqual(row['mean_occupancy'], (.2 + .2 + 0.) / 3)
        self.assertAlmostEqual(row['delta_occupancy_vs_WT'], (.2 + .2) / 3 - .8)

    def test_comparison_workflows_reach_every_figure_family(self):
        from matplotlib.figure import Figure
        for module, name, expected in ((plot_structure, 'structure', 12), (plot_binding, 'binding', 3)):
            with self.subTest(workflow=name), patch.object(Figure, 'savefig', autospec=True) as save:
                out = self.root / name
                module.run_comparison(self.paths, out)
                names = {Path(call.args[1]).name for call in save.call_args_list}
                self.assertEqual(len(names), expected * 2)
                self.assertTrue(all(path.suffix in ('.png', '.pdf') for path in map(Path, names)))
                with self.assertRaises(FileExistsError):
                    module.run_comparison(self.paths, out)
                self.assertTrue(any(out.glob('*.csv')))

    def test_sensitivity_summary_tables_and_figures(self):
        from matplotlib.figure import Figure
        root = self.root / 'sensitivity'
        for run in self.runs:
            source = run['meta']['input']
            out = root / source['system'] / f"rep{source['replicate']}"
            out.mkdir(parents=True)
            count_rows, pairs = [], []
            for kind, start, stride, cutoff in sensitivity.conditions():
                info = dict(system=source['system'], replicate=source['replicate'], kind=kind,
                            start_frame=start, stride=stride, cutoff_A=cutoff, analyzed_frames=3)
                occupancy = .8 if source['system'] == 'WT' else .2
                count_rows.append(dict(info, mean_pairs=occupancy))
                pairs.append(dict(info, Fc_label='A:265', receptor_label='C:120', occupancy=occupancy))
            write_csv(out / 'counts.csv', list(count_rows[0]), count_rows)
            write_csv(out / 'pairs.csv', list(pairs[0]), pairs)
            (out / 'complete.json').write_text('{}')
        with patch.object(Figure, 'savefig', autospec=True) as save:
            sensitivity.summarize(self.cfg, self.runs, root)
        self.assertEqual(len(save.call_args_list), 10)
        with (root / 'summary/binding_window_comparisons.csv').open() as fh:
            rows = list(csv.DictReader(fh))
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(abs(float(row['ddG_vs_WT']) + 10) < 1e-12 for row in rows))
        self.assertTrue((root / 'summary/complete.json').is_file())


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('MDAnalysis'),
                     'NumPy and MDAnalysis are required for dense contact tests')
class SensitivityScanTests(unittest.TestCase):
    def test_dense_scan_matches_primary_counts_and_occupancy(self):
        import numpy as np
        import MDAnalysis as mda
        from MDAnalysis.coordinates.memory import MemoryReader
        from analysis.io import digest
        with SharedFilesystemTemporaryDirectory(dir=HERE) as folder:
            root = Path(folder)
            topology = root / 'topology'
            topology.write_text('synthetic topology')
            mapping = [dict(label='A:265', name='ASP', partner='Fc'),
                       dict(label='C:120', name='LYS', partner='FcGR3a')]
            write_csv(root / 'residue_map.csv', list(mapping[0]), mapping)
            primary = [dict(kind=kind, Fc_label='A:265', receptor_label='C:120', observed_frames=140)
                       for kind in ('contact', 'salt_bridge')]
            write_csv(root / 'pairs.csv', list(primary[0]), primary)
            universe = mda.Universe.empty(4, n_residues=2, atom_resindex=[0, 0, 1, 1])
            universe.add_TopologyAttr('names', ['OD1', 'OD2', 'NZ', 'H1'])
            universe.add_TopologyAttr('resnames', ['ASP', 'LYS'])
            coordinates = np.repeat(np.array([[[0., 0., 0.], [0., .5, 0.],
                                                 [3., 0., 0.], [0., 0., 0.]]], dtype=np.float32), 1000, axis=0)
            universe.load_new(coordinates, format=MemoryReader, dt=10.)
            source = dict(system='WT', replicate=1, topology=str(topology), trajectory='synthetic')
            run = dict(path=str(root), meta=dict(input=source, topology_sha256=digest(topology),
                                               protocol_sha256='synthetic'))
            out = root / 'result'
            with patch('MDAnalysis.Universe', return_value=universe):
                sensitivity.analyze_run(run, out)
            with (out / 'counts.csv').open() as fh:
                counts = list(csv.DictReader(fh))
            self.assertEqual(len(counts), 36)
            self.assertTrue(all(float(row['mean_pairs']) == 1. for row in counts))
            primary_count = next(row for row in counts if row['kind'] == 'contact' and
                                 row['start_frame'] == '301' and row['stride'] == '5' and row['cutoff_A'] == '4.5')
            self.assertEqual(int(primary_count['analyzed_frames']), 140)
            metadata = json.loads((out / 'complete.json').read_text())
            self.assertEqual(metadata['dense_frames'], 700)
            self.assertTrue(metadata['baseline_counts_match'])
            with self.assertRaisesRegex(ValueError, 'baseline differs'):
                from collections import Counter
                counters = {condition: Counter() for condition in sensitivity.conditions()}
                sensitivity.validate_primary_counts(run, mapping, counters)


if __name__ == '__main__':
    unittest.main()
