"""End-to-end report tests with tiny synthetic summaries; no trajectory reads."""
import base64
import csv
from html.parser import HTMLParser
import json
from pathlib import Path
import unittest

from analysis.build_report import build, figure_catalog
from test_analysis import SharedFilesystemTemporaryDirectory


class Links(HTMLParser):
    def __init__(self):
        super().__init__()
        self.links = []
        self.ids = set()
        self.images = 0

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if 'id' in attrs:
            if attrs['id'] in self.ids:
                raise AssertionError('duplicate ID')
            self.ids.add(attrs['id'])
        if tag == 'img' and attrs.get('src'):
            self.images += 1
        self.links.extend(attrs[k] for k in ('href', 'src') if attrs.get(k))


def write(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x', newline='') as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


class ReportTests(unittest.TestCase):
    def test_two_pages_retain_all_figures_and_derive_numbers_from_tables(self):
        with SharedFilesystemTemporaryDirectory(dir=Path(__file__).parent) as folder:
            root = Path(folder)
            project = root / 'renamed-source'
            for figure in figure_catalog(project):
                figure['source'].parent.mkdir(parents=True, exist_ok=True)
                figure['source'].write_bytes(base64.b64decode('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+jZqkAAAAASUVORK5CYII='))
                figure['source'].with_suffix('.pdf').write_bytes(b'%PDF-1.4\nsynthetic placeholder\n')
            systems = ['WT', 'D265A', 'A327S', 'R255A']
            metrics = ['receptor_common_reference_rmsd_A', 'receptor_self_fit_rmsd_A',
                       'partner_COM_distance_A', 'complex_Rg_A']
            write(project / 'defined_figures_review/replicate_means.csv',
                  [dict(system=s, replicate=n, **{m: 8+n for m in metrics}) for s in systems for n in (1, 2, 3)])
            write(project / 'defined_binding_figures/comparisons.csv',
                  [dict(system=s, ddG_vs_WT=11.23 if s == 'D265A' else 0,
                        ci95_lo=5, ci95_hi=15, p_holm=.0645) for s in systems])
            write(project / 'defined_binding_figures/replicates.csv',
                  [dict(system=s, replicate=n, EEL=100 if s == 'D265A' else 0,
                        EGB=-90 if s == 'D265A' else 0) for s in systems for n in (1, 2, 3)])
            pairs = []
            for kind, variant, fc, rc in [('salt_bridge', 'D265A', 'A:265', 'C:120'),
                                          ('contact', 'A327S', 'A:N297:glycan2', 'C:N162:glycan1'),
                                          ('contact', 'R255A', 'B:332', 'C:161')]:
                for system, occupancy in [('WT', .9), (variant, .1)]:
                    pairs.append(dict(kind=kind, system=system, Fc_label=fc, receptor_label=rc,
                                      mean_occupancy=occupancy, sd_between_replicates=.02,
                                      delta_occupancy_vs_WT=occupancy-.9))
            write(project / 'defined_figures_review/pair_comparisons.csv', pairs)
            summary = project / 'sensitivity_results/summary'
            write(summary / 'finding_ranges.csv',
                  [dict(kind=p['kind'], system=p['system'], Fc_label=p['Fc_label'],
                        receptor_label=p['receptor_label'], minimum_delta_pp=-95,
                        maximum_delta_pp=-70, direction_retained=True) for p in pairs if p['system'] != 'WT'])
            write(summary / 'binding_window_comparisons.csv',
                  [dict(system=s, start_frame=n, ddG_vs_WT=value) for s in ['D265A', 'S298A', 'A327S', 'E333A']
                   for n, value in [(301, -.5), (501, 1.9)]])
            write(summary / 'sampling_comparison.csv',
                  [dict(kind=k, max_abs_difference_pp=1) for k in ['contact', 'salt_bridge']])
            out = root / 'report'
            build(project, out)
            main = (out / 'index.html').read_text()
            self.assertIn('+11.23 kcal/mol', main)  # From supplied table, not a canned result.
            self.assertIn('10.00 ± 1.00', main)
            self.assertIn('No mutant comparison passes', main)
            docs = {}
            for name in ('index.html', 'supplemental.html'):
                docs[name] = Links()
                docs[name].feed((out / name).read_text())
            self.assertEqual(docs['index.html'].images, 1)
            self.assertEqual(docs['supplemental.html'].images, 23)
            for name, document in docs.items():
                for link in document.links:
                    self.assertFalse(link.startswith(('http:', 'https:', '//', '/')))
                    path, _, fragment = link.partition('#')
                    target = out / (path or name)
                    self.assertTrue(target.is_file(), link)
                    if fragment:
                        self.assertIn(fragment, docs[target.name].ids)
            manifest = json.loads((out / 'manifest.json').read_text())
            self.assertEqual(len(manifest['figures']), 23)
            self.assertTrue(all(not Path(f['source']).is_absolute() for f in manifest['figures']))
            with self.assertRaises(FileExistsError):
                build(project, out)


if __name__ == '__main__':
    unittest.main()
