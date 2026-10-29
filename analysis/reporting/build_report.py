#!/usr/bin/env python3
"""Build a portable, offline annotated figure report. No server or external assets."""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
from html import escape
import json
from pathlib import Path
import zipfile


def read_csv(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh))


def figure_catalog(project):
    structural = project / 'defined_figures_review'
    binding = project / 'defined_binding_figures'
    sensitivity = project / 'sensitivity_results/summary'
    figures = []
    def add(folder, name, category, title, question, reading, interpretation, limitation):
        figures.append(dict(source=folder / (name + '.png'), category=category, title=title,
                            question=question, reading=reading, interpretation=interpretation,
                            limitation=limitation))
    add(structural, 'structure_and_interface', 'structure', 'Complex structure and interface overview',
        'Do the variants differ in complex geometry, internal deformation, or interface contacts?',
        'Each dot is one replicate mean. Bars show the system mean ± between-replicate SD. '
        'RMSDs use a common WT backbone reference; the panel titles distinguish CH3 alignment from self-fitting.',
        'System means for complex size and partner separation are similar, while domain/receptor RMSDs and local interactions vary.',
        'SD bars are not confidence intervals. Similar global dimensions do not establish equivalent structures or interfaces.')
    add(binding, 'binding_comparison', 'energies', 'Whole-partner MM/GBSA comparison',
        'Which variants have a reproducible shift in the calculated binding-associated energy?',
        'Left: replicate means with between-replicate SD. Right: variant minus WT with unadjusted 95% Welch intervals. '
        'Positive differences are less favorable within this model.',
        'D265A has the clearest unfavorable shift (+10.02 kcal/mol); its Holm-adjusted p-value is 0.0645. '
        'All other unadjusted intervals cross zero.',
        'These are whole-partner MM/GBSA estimates without entropy, not measured affinities or interface-only energies.')
    add(binding, 'energy_components', 'energies', 'Electrostatic compensation and other components',
        'Which calculated components account for the mean differences from WT?',
        'The left panel separates direct electrostatics (EEL) and GB solvation (EGB). '
        'The right panel shows their sum beside van der Waals and surface contributions. Positive values are unfavorable differences.',
        'For D265A, approximately +148.4 kcal/mol in EEL is offset by −141.0 kcal/mol in EGB, leaving +7.39 kcal/mol net polar.',
        'Bars are descriptive means without uncertainty intervals. They do not assign energies to particular interface residue pairs.')
    add(binding, 'contacts_vs_binding', 'energies', 'Contact counts versus whole-partner energy',
        'Do replicate-level contact counts covary with the calculated energy?',
        'Each point is a replicate; colors identify systems. The two panels use total heavy-atom contact pairs and salt-bridge candidate pairs.',
        'Use this as an exploratory comparison of structural and energetic variation, rather than fitting a mechanism from the pooled points.',
        'Between-system and within-system variation are mixed; contact counts do not measure interaction strength.')
    for page in range(1, 5):
        add(structural, f'contact_changes_{page}', 'interface', f'Heavy-atom contact changes — page {page} of 4',
            'Which protein and glycan residue pairs gain or lose sampled contact occupancy?',
            'Rows identify original chain/residue labels; columns are variants. Red means more frequent contact than WT; '
            'blue means less frequent. Values are percentage-point changes, with a fixed ±100 scale. '
            'Only pairs reaching 20% mean occupancy in any system are displayed.',
            'Pages continue the same list, ordered by largest absolute primary change. '
            'Page 1 includes the large A327S change in the Fc A N297 glycan residue 2–receptor N162 NAG contact.',
            'A: and B: are Fc chains; C: is receptor. Glycan numbers label individual rings, not the whole glycan. '
            'These primary heatmaps use 4.5 Å and 50 ps sampling; use the sensitivity figures to judge dependence on those choices.')
    add(structural, 'salt_bridge_changes_1', 'interface', 'Salt-bridge candidate occupancy changes',
        'Which charged interface residue pairs change their geometric salt-bridge occupancy?',
        'Colors show variant minus WT occupancy. A:265–C:120 is Fc chain A Asp265–receptor Lys120; '
        'B:265–C:161 is Fc chain B Asp265–receptor Lys161.',
        'The chain A Asp265–Lys120 candidate occurs in 97.4% of primary WT samples and is absent in D265A.',
        'Loss after removing Asp charge is an expected chemical consequence, not an independent estimate of its energetic contribution. '
        'The criterion is distance-based and assumes fixed protonation.')
    timecourses = [
        ('CH3_common_reference_rmsd_A', 'Fc CH3 backbone RMSD',
         'Residual displacement of the CH3 domains after fitting them jointly to the WT reference.'),
        ('CH2_common_reference_rmsd_A', 'Fc CH2 backbone RMSD',
         'CH2 displacement in the CH3 reference frame; includes domain motion and internal deformation.'),
        ('receptor_common_reference_rmsd_A', 'Receptor backbone RMSD after CH3 alignment',
         'Receptor displacement relative to Fc CH3; includes relative motion and receptor deformation.'),
        ('Fc_self_fit_rmsd_A', 'Fc internal deformation',
         'Backbone RMSD after fitting Fc to WT independently; overall Fc translation/rotation is removed.'),
        ('receptor_self_fit_rmsd_A', 'Receptor internal deformation',
         'Backbone RMSD after fitting receptor independently, separating deformation from overall receptor motion.'),
        ('partner_COM_distance_A', 'Glycosylated partner separation',
         'Mass-weighted COM distance between complete modeled partners, including their glycans.'),
        ('complex_Rg_A', 'Complex radius of gyration',
         'Mass-weighted size of the whole dry complex; a global compactness measure.'),
        ('contact_residue_pairs', 'Interface heavy-atom contact count',
         'Number of distinct inter-partner residue pairs with at least one heavy-atom distance within 4.5 Å.'),
        ('salt_bridge_residue_pairs', 'Interface salt-bridge candidate count',
         'Number of distinct charged residue pairs meeting the 4.0 Å geometric criterion.')]
    for name, title, meaning in timecourses:
        add(structural, name + '_timecourses', 'structure' if 'residue_pairs' not in name else 'interface',
            title + ' — individual trajectories', 'Are there trends, shifts, or persistent differences among replicates?',
            meaning + ' Each panel is a system and each line is a replicate. Stored time includes equilibration; '
            'the primary samples span 3.16–10.11 ns at 50 ps intervals.',
            'Inspect late drift and replicate disagreement before interpreting the corresponding mean. '
            'Compare this panel with the window-sensitivity figure.',
            'Visual stability does not prove equilibration or convergence. Short-lived events between sampled frames can be missed.')
    sensitivity_cards = [
        ('contact_count_sensitivity', 'Contact-count sensitivity to cutoffs and windows',
         'Do total interface differences retain their direction when cutoffs or start frames change?',
         'Rows separate heavy-atom contacts and salt-bridge candidates. Left panels vary cutoff at the primary window; '
         'right panels vary start frame at the primary cutoff. All use 50 ps samples. Lines identify variants relative to WT.',
         'Changes across cutoffs indicate definition dependence; changes across start times indicate temporal dependence.'),
        ('sampling_sensitivity', '10 ps versus 50 ps occupancy comparisons',
         'Does the coarser primary sampling change residue-pair occupancy differences?',
         'Each point is a variant/residue-pair comparison. The diagonal indicates agreement between cadences. '
         'Panel titles give the largest absolute deviation, in percentage points.',
         'Points close to the diagonal support the primary sampling for occupancy, but do not establish contact lifetimes.'),
        ('finding_robustness', 'Largest primary contact changes across analysis choices',
         'Which highlighted findings persist across the full cutoff/window/cadence grid?',
         'Rows were selected by largest primary changes, not by sensitivity performance. Each row uses its interaction-specific cutoffs. '
         'Column labels are start frame / stride / cutoff index. Red and blue retain the same meaning as the primary heatmaps.',
         'A stable sign is evidence of directional consistency across these settings; the CSV gives minimum and maximum effects.'),
        ('binding_window_sensitivity', 'MM/GBSA sensitivity to the start of the window',
         'Does the mean energetic ordering depend on how much early production is excluded?',
         'Each colored line uses a different start frame with the same ending frame and existing 50 ps energy snapshots. '
         'Positive values indicate less favorable energy than WT.',
         'Check whether effects retain their sign and how much their magnitude changes. No new energy calculations are performed.'),
        ('structure_window_sensitivity', 'Structural sensitivity to the start of the window',
         'Do structural mean differences change in later subsets of the trajectories?',
         'Each panel is a structural observable; lines identify starting frames and show variant minus WT means at 50 ps sampling.',
         'Substantial movement across windows means a single full-window average is incomplete.')]
    for name, title, question, reading, interpretation in sensitivity_cards:
        add(sensitivity, name, 'sensitivity', title, question, reading, interpretation,
            'These are nested subsets of the same three trajectories per system, not additional independent replicates. '
            'Sensitivity checks are descriptive and do not prove convergence, equivalence, or a causal mechanism.')
    return figures


def sensitivity_summary(project):
    folder = project / 'sensitivity_results/summary'
    sampling = read_csv(folder / 'sampling_comparison.csv')
    binding = read_csv(folder / 'binding_window_comparisons.csv')
    findings = read_csv(folder / 'finding_ranges.csv')
    lines = []
    for kind, label in [('contact', 'Heavy-atom'), ('salt_bridge', 'Salt-bridge candidate')]:
        values = [float(r['max_abs_difference_pp']) for r in sampling if r['kind'] == kind]
        lines.append(f'{label} occupancy differences: the largest 10 ps versus 50 ps deviation is {max(values):.2f} percentage points.')
    for system in ('D265A', 'S298A', 'A327S'):
        values = [float(r['ddG_vs_WT']) for r in binding if r['system'] == system]
        lines.append(f'{system} mean MM/GBSA difference from WT ranges from {min(values):+.2f} to {max(values):+.2f} kcal/mol across start frames.')
    kept = sum(r['direction_retained'] == 'True' for r in findings)
    lines.append(f'{kept} of {len(findings)} highlighted primary variant/pair changes retain direction across all 18 settings per interaction type. '
                 'This is a descriptive sign check, not a significance test.')
    return lines


def prepare_report_inputs(project):
    """Locate and validate every figure and table before creating outputs."""
    figures = figure_catalog(project)
    tables = [('primary_structural_means.csv', project / 'defined_figures_review/replicate_means.csv'),
              ('primary_pair_comparisons.csv', project / 'defined_figures_review/pair_comparisons.csv'),
              ('primary_binding_comparisons.csv', project / 'defined_binding_figures/comparisons.csv')]
    tables += [(p.name, p) for p in sorted((project / 'sensitivity_results/summary').glob('*.csv'))]
    for i, fig in enumerate(figures, 1):
        fig.update(png=f'{i:02d}_{fig["source"].stem}.png', pdf=f'{i:02d}_{fig["source"].stem}.pdf')
        for source in (fig['source'], fig['source'].with_suffix('.pdf')):
            if not source.is_file():
                raise FileNotFoundError(source)
    for _, source in tables:
        if not source.is_file():
            raise FileNotFoundError(source)
    return figures, tables


def render_report_pages(project, figures, tables):
    """Render both offline pages from validated summary tables."""
    from .report_pages import home, supplemental

    # Read/validate summary inputs before creating outputs.
    main_html = home(project, figures)
    tablelinks = [f'<li><a href="assets/{escape(name)}">{escape(name)}</a></li>' for name, _ in tables]
    supplemental_html = supplemental(figures, tablelinks)
    return main_html, supplemental_html


def copy_report_assets(project, out, figures, tables):
    """Copy selected assets and compute their manifest entries."""
    import shutil

    out.mkdir(parents=True, exist_ok=False)
    assets = out / 'assets'
    assets.mkdir()
    manifest = []
    for i, fig in enumerate(figures, 1):
        for suffix, key in (('.png', 'png'), ('.pdf', 'pdf')):
            with fig['source'].with_suffix(suffix).open('rb') as src, (assets / fig[key]).open('xb') as dst:
                shutil.copyfileobj(src, dst)
        manifest.append(dict(number=i, **{k: str(v) for k, v in fig.items() if k != 'source'},
                             source=fig['source'].relative_to(project).as_posix(),
                             png_sha256=hashlib.sha256((assets / fig['png']).read_bytes()).hexdigest()))
    for name, source in tables:
        with source.open('rb') as src, (assets / name).open('xb') as dst:
            shutil.copyfileobj(src, dst)
    return manifest


def write_report_files(out, main_html, supplemental_html, manifest, summary):
    """Write pages and the timestamped report manifest."""
    for name, content in [('index.html', main_html), ('supplemental.html', supplemental_html)]:
        with (out / name).open('x') as fh:
            fh.write(content)
    with (out / 'manifest.json').open('x') as fh:
        json.dump(dict(created_utc=datetime.now(timezone.utc).isoformat(), figures=manifest,
                       sensitivity_summary=summary, local_only=True,
                       pages=['index.html', 'supplemental.html']), fh, indent=2)


def archive_report(out, archive):
    """Create a new ZIP archive without replacing an existing file."""
    with zipfile.ZipFile(archive, 'x', compression=zipfile.ZIP_DEFLATED) as z:
        for path in sorted(out.rglob('*')):
            if path.is_file():
                z.write(path, arcname=Path(out.name) / path.relative_to(out))


def build(project, out, archive=None):
    """Validate, render, copy assets, and write the portable report."""
    project = Path(project).resolve()
    out = Path(out)
    if out.exists() or (archive and Path(archive).exists()):
        raise FileExistsError('report directory/archive must be new')
    figures, tables = prepare_report_inputs(project)
    main_html, supplemental_html = render_report_pages(project, figures, tables)
    summary = sensitivity_summary(project)
    manifest = copy_report_assets(project, out, figures, tables)
    write_report_files(out, main_html, supplemental_html, manifest, summary)
    if archive:
        archive_report(out, archive)
    print(f'Concise home page + {len(figures)} supplemental figures -> {out / "index.html"}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--zip', type=Path)
    args = parser.parse_args()
    build(args.project, args.out, args.zip)


if __name__ == '__main__':
    main()
