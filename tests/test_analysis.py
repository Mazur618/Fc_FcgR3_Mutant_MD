"""Small workflow/numerical tests; no simulation or trajectory reads."""
import json
import importlib.util
import subprocess
import sys
import errno
from pathlib import Path
import tempfile
import time
import unittest

from analysis.trajectory import input_records, kabsch_fit, pair_counts, read_protocol, residue_map, write_csv
from analysis.plot_structure import load_runs
from analysis.plot_binding import load_energy
from analysis.sensitivity import conditions, frame_selected, minimum_pairs


HERE = Path(__file__).resolve().parent


class SharedFilesystemTemporaryDirectory(tempfile.TemporaryDirectory):
    """Retry transient ENOTEMPTY when shared storage delays directory updates."""
    def cleanup(self):
        for attempt in range(6):
            try:
                return super().cleanup()
            except OSError as error:
                if error.errno != errno.ENOTEMPTY or attempt == 5:
                    raise
                time.sleep(0.2 * (attempt + 1))


class ProtocolTests(unittest.TestCase):
    def setUp(self):
        self.protocol, self.cfg = read_protocol(HERE.parent / 'config/analysis_protocol.json')

    def test_mapping_uses_original_numbering_and_separates_glycans(self):
        # Synthetic reference: exercises numbering without distributing study data.
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            reference = Path(d) / 'reference.pdb'
            lines = []
            serial = 0
            for chain, first, count in [('A', 233, 212), ('B', 232, 213), ('C', 1, 171)]:
                for number in range(first, first + count):
                    name = 'VAL' if (chain, number) == ('C', 158) else 'ALA'
                    for atom in ('N', 'CA', 'C'):
                        serial += 1
                        lines.append(f'ATOM  {serial:5d} {atom:^4s} {name} {chain}{number:4d}    '
                                     f'{0.:8.3f}{0.:8.3f}{0.:8.3f}  1.00  0.00           C\n')
            reference.write_text(''.join(lines))
            rows = residue_map(reference, self.cfg)
        self.assertEqual(len(rows), 614)
        self.assertEqual(rows[64]['label'], 'A:297')
        self.assertEqual(rows[277]['label'], 'B:297')
        self.assertEqual(rows[596]['label'], 'A:N297:glycan1')
        self.assertEqual(rows[604]['label'], 'B:N297:glycan1')
        self.assertEqual(rows[612]['label'], 'C:N162:glycan1')
        self.assertEqual(rows[613]['label'], 'C:N45:glycan1')
        val158 = next(r for r in rows if r['label'] == 'C:158')
        self.assertEqual(val158['name'], 'VAL')

    def test_window_is_same_140_snapshots_as_mmgbsa(self):
        frames = list(range(self.cfg['start_frame'], self.cfg['end_frame'] + 1,
                            self.cfg['stride']))
        self.assertEqual((len(frames), frames[0], frames[-1]), (140, 301, 996))

    def test_sensitivity_windows_and_cadences(self):
        self.assertEqual(len(conditions()), 36)
        for start, dense, coarse in [(301, 700, 140), (501, 500, 100), (701, 300, 60)]:
            self.assertEqual(sum(frame_selected(f, start, 1) for f in range(1, 1001)), dense)
            self.assertEqual(sum(frame_selected(f, start, 5) for f in range(1, 1001)), coarse)
        self.assertFalse(frame_selected(1001, 301, 5))

    def test_atom_contacts_are_counted_once_per_residue_pair(self):
        pairs = [(0, 0), (1, 0), (1, 1), (2, 2)]
        self.assertEqual(pair_counts(pairs, [64, 64, 277], [549, 549, 580]),
                         {(64, 549), (277, 580)})
        self.assertEqual(pair_counts([], [], []), set())

    def test_inventory_resolves_wt_symlink_and_does_not_duplicate_runs(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            root = Path(d)
            cfg = dict(self.cfg, prep_root='.')
            folder = root / 'WT' / 'rep1'
            folder.mkdir(parents=True)
            target = root / 'example.nc'
            target.touch()
            (folder / 'prod_dry.nc').symlink_to('../../example.nc')
            records = input_records(root / 'protocol.json', cfg)
            self.assertEqual(len(records), 24)
            self.assertEqual(len({r['trajectory'] for r in records}), 24)
            wt = next(r for r in records if r['system'] == 'WT' and r['replicate'] == 1)
            self.assertEqual(Path(wt['trajectory']), target.resolve())
            self.assertTrue(wt['trajectory_exists'])

    def test_duplicate_trajectory_is_rejected(self):
        self.cfg['replicates'] = [1, 1]
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            input_records(self.protocol, self.cfg)

    def test_invalid_window_is_rejected(self):
        cfg = dict(self.cfg, stride=0)
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            path = Path(d) / 'invalid.json'
            path.write_text(json.dumps(cfg))
            with self.assertRaisesRegex(ValueError, 'stride'):
                read_protocol(path)

    def test_existing_output_is_protected(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            path = Path(d) / 'table.csv'
            write_csv(path, ['value'], [dict(value=1)])
            with self.assertRaises(FileExistsError):
                write_csv(path, ['value'], [dict(value=2)])
            self.assertIn('1', path.read_text())

    def test_failed_run_cannot_enter_comparison(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            with self.assertRaises(FileNotFoundError):
                load_runs([d])

    def test_complete_run_with_missing_replicates_is_rejected(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            folder = Path(d)
            meta = dict(input=dict(system='WT', replicate=1, trajectory='example.nc'),
                        protocol=self.cfg, protocol_sha256='example',
                        reference_sha256='example', n_frames=1)
            (folder / 'complete.json').write_text(json.dumps(meta))
            write_csv(folder / 'pairs.csv', ['kind', 'Fc_label', 'receptor_label', 'occupancy'], [])
            write_csv(folder / 'frames.csv', ['frame', 'time_ns'], [dict(frame=301, time_ns=3.01)])
            with self.assertRaisesRegex(ValueError, 'missing='):
                load_runs([folder])

    def test_energy_parser_selects_delta_instead_of_complex(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            path = Path(d) / 'energy.csv'
            path.write_text('Complex Energy Terms\nFrame #,TOTAL\n0,-999\n\n'
                            'DELTA Energy Terms\nFrame #,VDWAALS,EEL,EGB,ESURF,DELTA TOTAL\n'
                            '0,-100,-300,370,-15,-45\n1,-101,-301,371,-15,-46\n')
            rows = load_energy(path)
            self.assertEqual([r['DELTA TOTAL'] for r in rows], [-45, -46])

    def test_energy_parser_rejects_missing_snapshot_ordinals(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            path = Path(d) / 'energy.csv'
            path.write_text('DELTA Energy Terms\nFrame #,VDWAALS,EEL,EGB,ESURF,DELTA TOTAL\n'
                            '0,-100,-300,370,-15,-45\n2,-101,-301,371,-15,-46\n')
            with self.assertRaisesRegex(ValueError, 'nonconsecutive'):
                load_energy(path)


class PortabilityTests(unittest.TestCase):
    def test_package_commands_show_help_without_analysis_dependencies(self):
        for module in ('trajectory', 'plot_structure', 'plot_binding',
                       'sensitivity', 'build_report'):
            with self.subTest(module=module):
                result = subprocess.run(
                    [sys.executable, '-B', '-m', 'analysis.' + module, '--help'],
                    cwd=HERE.parent, capture_output=True, text=True)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)

    def test_protocol_paths_survive_moving_the_dataset_tree(self):
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            root = Path(d)
            source = root / 'original tree'
            source.mkdir()
            (source / 'data').mkdir()
            (source / 'reference.pdb').touch()
            protocol = source / 'configured.json'
            command = [sys.executable, str(HERE.parent / 'slurm/configure_protocol.py'),
                       '--data-root', str(source / 'data'),
                       '--reference', str(source / 'reference.pdb'), '--out', str(protocol)]
            subprocess.run(command, check=True, capture_output=True, text=True)
            moved = root / 'moved tree'
            source.rename(moved)
            path, cfg = read_protocol(moved / 'configured.json')
            self.assertEqual((path.parent / cfg['prep_root']).resolve(), moved / 'data')
            self.assertTrue((path.parent / cfg['reference_pdb']).is_file())
            command[-1] = str(moved / 'configured.json')
            refused = subprocess.run(command, capture_output=True, text=True)
            self.assertNotEqual(refused.returncode, 0)
            self.assertIn('FileExistsError', refused.stderr)

    def test_dataset_export_materializes_files_and_hashes_trajectories(self):
        spec = importlib.util.spec_from_file_location('export_inputs', HERE.parent / 'slurm/export_inputs.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with SharedFilesystemTemporaryDirectory(dir=HERE) as d:
            root = Path(d)
            source = root / 'source'
            source.mkdir()
            (source / 'example.nc').write_bytes(b'synthetic trajectory placeholder')
            manifest = dict(entries=[dict(path='WT/rep1/prod_dry.nc', source='example.nc',
                                         group='analysis', size_bytes=32, sha256=None)])
            manifest['entries'][0]['size_bytes'] = (source / 'example.nc').stat().st_size
            out = root / 'dataset'
            copied = module.export(manifest, source, out, {'analysis'})
            self.assertEqual((out / 'WT/rep1/prod_dry.nc').read_bytes(), (source / 'example.nc').read_bytes())
            self.assertFalse((out / 'WT/rep1/prod_dry.nc').is_symlink())
            self.assertEqual(len(copied[0]['sha256']), 64)
            self.assertTrue((out / 'dataset-manifest.json').is_file())
            with self.assertRaises(FileExistsError):
                module.export(manifest, source, out, {'analysis'})
            manifest['entries'][0]['sha256'] = 'incorrect'
            with self.assertRaisesRegex(ValueError, 'checksum mismatch'):
                module.export(manifest, source, root / 'bad-hash', {'analysis'})
            manifest['entries'][0]['path'] = '../escape.nc'
            with self.assertRaisesRegex(ValueError, 'unsafe'):
                module.selected_inputs(manifest, source, {'analysis'})

    def test_report_paths_do_not_assume_a_prep_parent(self):
        from analysis.build_report import figure_catalog
        root = Path('renamed-repository')
        figures = figure_catalog(root)
        self.assertEqual(len(figures), 23)
        for figure in figures:
            relative = figure['source'].relative_to(root)
            self.assertIn(relative.parts[0], ('defined_figures_review',
                                             'defined_binding_figures',
                                             'sensitivity_results'))


class NumericalTests(unittest.TestCase):
    def setUp(self):
        try:
            import numpy as np
        except ImportError:
            self.skipTest('NumPy is required for numerical tests')
        self.np = np

    def test_fit_recovers_rotation_and_translation(self):
        np = self.np
        x = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        rotation = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]])
        y = x @ rotation + [10., -5., 2.]
        fitted, xc, yc = kabsch_fit(x, y)
        np.testing.assert_allclose((x - xc) @ fitted + yc, y, atol=1e-12)

    def test_fit_does_not_reflect_chirality(self):
        np = self.np
        x = np.array([[0., 0., 0.], [1., 0., 0.], [0., 2., 0.], [0., 0., 3.]])
        y = x * [-1., 1., 1.]
        rotation, xc, yc = kabsch_fit(x, y)
        self.assertAlmostEqual(np.linalg.det(rotation), 1.)
        self.assertGreater(float(np.linalg.norm((x - xc) @ rotation + yc - y)), 0.1)

    def test_capped_contacts_match_brute_force_distances(self):
        np = self.np
        try:
            from MDAnalysis.lib.distances import capped_distance
        except ImportError:
            self.skipTest('MDAnalysis is required for contact tests')
        left = np.array([[0., 0., 0.], [1., 0., 0.], [10., 0., 0.]])
        right = np.array([[0., 2., 0.], [10., 1., 0.], [30., 0., 0.]])
        actual = capped_distance(left, right, 4.5, return_distances=False)
        distances = np.linalg.norm(left[:, None, :] - right[None, :, :], axis=2)
        expected = np.argwhere(distances <= 4.5)
        self.assertEqual({tuple(p) for p in actual}, {tuple(p) for p in expected})
        self.assertEqual(pair_counts(actual, [64, 64, 277], [549, 580, 581]),
                         {(64, 549), (277, 580)})

    def test_minimum_distance_reduction_and_cutoff_monotonicity(self):
        np = self.np
        keys, distance = minimum_pairs(np.array([[0, 0], [1, 0], [1, 1], [2, 2]]),
                                       np.array([4.8, 4.1, 3.8, 4.6]),
                                       np.array([64, 64, 277]), np.array([549, 549, 580]), 614)
        self.assertEqual(dict(zip(map(int, keys), distance)), {64 * 614 + 549: 3.8, 277 * 614 + 580: 4.6})
        self.assertEqual([int((distance <= c).sum()) for c in (4., 4.5, 5.)], [1, 1, 2])
        empty, _ = minimum_pairs(np.empty((0, 2), dtype=int), np.array([]), np.array([]), np.array([]), 614)
        self.assertEqual(len(empty), 0)


if __name__ == '__main__':
    unittest.main()
