"""Synthetic trajectory regression checks; no study inputs or production runs."""
import csv
import importlib.util
import json
from pathlib import Path
import unittest
from unittest.mock import patch

from analysis.structure import trajectory
from test_analysis import SharedFilesystemTemporaryDirectory


@unittest.skipUnless(importlib.util.find_spec('numpy') and importlib.util.find_spec('MDAnalysis'),
                     'NumPy and MDAnalysis are required for trajectory tests')
class TrajectoryTests(unittest.TestCase):
    def setUp(self):
        import numpy as np
        import MDAnalysis as mda
        from MDAnalysis.coordinates.memory import MemoryReader

        self.temp = SharedFilesystemTemporaryDirectory(dir=Path(__file__).parent)
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        _, cfg = trajectory.read_protocol(Path(__file__).parents[1] / 'config/analysis_protocol.json')
        self.cfg = dict(cfg, systems=['WT'], replicates=[1], prep_root='.',
                        reference_pdb='reference.pdb', start_frame=1, end_frame=5, stride=2)
        self.protocol = self.root / 'protocol.json'
        self.protocol.write_text(json.dumps(self.cfg))
        folder = self.root / 'WT/rep1'
        folder.mkdir(parents=True)
        # Placeholders supply provenance hashes; the patched reader supplies synthetic coordinates.
        (folder / 'dry.prmtop').write_text('synthetic topology')
        (folder / 'prod_dry.nc').touch()
        atom_names, residue_indices, positions, reference, residue_names = [], [], [], [], []
        for chain, first, count in [('A', 233, 212), ('B', 232, 213), ('C', 1, 171)]:
            for number in range(first, first + count):
                index = len(residue_names)
                name = 'ASP' if (chain, number) == ('A', 265) else 'LYS' if (chain, number) == ('C', 120) else 'ALA'
                residue_names.append(name)
                base = np.array([(index % 20) * 12., (index // 20) * 12., index % 7 * 3.])
                for atom, offset in [('N', (0., 0., 0.)), ('CA', (1., 0., 0.)), ('C', (1., 1., .5))]:
                    position = base + offset
                    atom_names.append(atom)
                    residue_indices.append(index)
                    positions.append(position)
                    reference.append(f'ATOM  {len(atom_names):5d} {atom:^4s} {name} {chain}{number:4d}    '
                                     f'{position[0]:8.3f}{position[1]:8.3f}{position[2]:8.3f}  1.00  0.00           C\n')
                if (chain, number) in [('A', 265), ('C', 120)]:
                    atom_names.extend(['OD1', 'OD2'] if chain == 'A' else ['NZ'])
                    n = 2 if chain == 'A' else 1
                    residue_indices.extend([index] * n)
                    positions.extend([[0., 0., 100.], [0., .5, 100.]] if chain == 'A' else [[3., 0., 100.]])
        for name in cfg['expected_glycan_names']:
            residue_indices.append(len(residue_names))
            residue_names.append(name)
            atom_names.append('C1')
            positions.append([(len(residue_names) - 597) * 12., 400., 0.])
        # A hydrogen near the receptor must never create a heavy-atom contact.
        atom_names.append('H1')
        residue_indices.append(0)
        positions.append([3., 0., 100.])
        (self.root / 'reference.pdb').write_text(''.join(reference))
        universe = mda.Universe.empty(len(atom_names), n_residues=len(residue_names),
                                     atom_resindex=residue_indices)
        universe.add_TopologyAttr('names', atom_names)
        universe.add_TopologyAttr('resnames', residue_names)
        universe.add_TopologyAttr('masses', [12.] * len(atom_names))
        coordinates = np.repeat(np.asarray(positions, dtype=np.float32)[None, :, :], 5, axis=0)
        nz_index = atom_names.index('NZ')
        for frame, separation in enumerate((3., 5., 4.2, 6., 8.)):
            coordinates[frame, nz_index, 0] = separation
            coordinates[frame] += frame * .25
        universe.load_new(coordinates, format=MemoryReader, dt=5.)
        self.universe = universe

    def run_analysis(self, implementation=trajectory, out=None):
        out = out or self.root / 'results'
        with patch('MDAnalysis.Universe', return_value=self.universe):
            implementation.run(self.protocol, self.cfg, 'WT', 1, out)
        return out

    def test_sampling_metrics_contacts_and_output_contract(self):
        out = self.run_analysis()
        with (out / 'frames.csv').open() as fh:
            frames = list(csv.DictReader(fh))
        self.assertEqual([int(row['frame']) for row in frames], [1, 3, 5])
        self.assertEqual([float(row['time_ns']) for row in frames], [0., .01, .02])
        for row in frames:
            for key, value in row.items():
                if 'rmsd' in key:
                    self.assertAlmostEqual(float(value), 0., places=10)
            self.assertGreater(float(row['complex_Rg_A']), 0.)
        with (out / 'pairs.csv').open() as fh:
            pairs = {(row['kind'], row['Fc_label'], row['receptor_label']): row
                     for row in csv.DictReader(fh)}
        key = ('A:265', 'C:120')
        self.assertEqual(int(pairs[('contact',) + key]['observed_frames']), 2)
        self.assertAlmostEqual(float(pairs[('contact',) + key]['occupancy']), 2 / 3)
        self.assertEqual(int(pairs[('salt_bridge',) + key]['observed_frames']), 1)
        self.assertAlmostEqual(float(pairs[('salt_bridge',) + key]['occupancy']), 1 / 3)
        self.assertNotIn(('contact', 'A:233', 'C:120'), pairs)
        metadata = json.loads((out / 'complete.json').read_text())
        self.assertEqual(metadata['n_frames'], 3)
        self.assertEqual(metadata['first_time_ns'], 0.)
        self.assertEqual(metadata['last_time_ns'], .02)
        self.assertTrue((out / 'residue_map.csv').is_file())
        with self.assertRaises(FileExistsError):
            self.run_analysis(out=out)

    def test_invalid_coordinates_do_not_produce_completed_results(self):
        self.universe.trajectory[0]
        self.universe.atoms[0].position = [float('nan'), 0., 0.]
        with self.assertRaisesRegex(ValueError, 'nonfinite'):
            self.run_analysis()
        self.assertFalse((self.root / 'results/complete.json').exists())

    def test_truncated_trajectory_is_rejected(self):
        self.cfg['end_frame'] = 6
        with self.assertRaisesRegex(ValueError, 'shorter than specified window'):
            self.run_analysis()
        self.assertFalse((self.root / 'results').exists())


if __name__ == '__main__':
    unittest.main()
