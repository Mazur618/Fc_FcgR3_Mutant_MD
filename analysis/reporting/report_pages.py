"""Concise offline report pages derived from existing summary tables."""
import csv
from html import escape
from pathlib import Path
import statistics

STYLE = '''
:root{font-family:system-ui,sans-serif;color:#20333c;background:#f4f5f3;line-height:1.6}
*{box-sizing:border-box}body{margin:0}header,main,footer{max-width:1120px;margin:auto;padding:24px}
header{padding-top:40px}h1{font-size:clamp(1.8rem,4vw,2.7rem);line-height:1.16;max-width:900px}
h2{font-size:1.5rem;line-height:1.3}h3{font-size:1.1rem}p{max-width:88ch}a{color:#145c75}
nav{display:flex;gap:20px;flex-wrap:wrap;margin:22px 0}nav a{font-weight:650;text-underline-offset:5px}
.eyebrow,.tag{text-transform:uppercase;letter-spacing:.09em;font-size:.78rem;color:#52646b}
.lead{font-size:1.15rem}.study{border-top:1px solid #cdd6d7;border-bottom:1px solid #cdd6d7;padding:12px 0}
section,article{background:#fff;border:1px solid #dce2e0;border-radius:8px;padding:26px;margin:0 0 24px}
.takeaway{border-left:4px solid #26736e}.note{font-size:.92rem;color:#4a6069}.table-scroll{overflow:auto}
table{border-collapse:collapse;width:100%;font-size:.91rem}th,td{text-align:left;padding:10px 12px;border-bottom:1px solid #e0e6e4;vertical-align:top}
th{background:#f0f4f2;font-weight:650}tbody tr:first-child{font-weight:600}td.number{white-space:nowrap;font-variant-numeric:tabular-nums}
figure{margin:16px 0}figure img{display:block;width:100%;height:auto}figcaption{font-size:.9rem;color:#4a6069;margin-top:10px}
details{margin:16px 0}summary{cursor:pointer;font-weight:600}ul{padding-left:22px}.finding{padding:12px 0;border-top:1px solid #e0e6e4}
.toolbar{display:flex;gap:8px;flex-wrap:wrap;margin:16px 0}button,input{font:inherit;border:1px solid #aebfbd;padding:8px 12px;border-radius:4px;background:#fff}
button{cursor:pointer}button.active{background:#215a58;color:white}input{flex:1;min-width:180px}
dialog{border:0;border-radius:8px;max-width:96vw;max-height:95vh;padding:16px}dialog::backdrop{background:#172b3acc}
dialog img{max-width:none}dialog .scroll{overflow:auto;max-height:80vh;max-width:90vw}[hidden]{display:none!important}
@media(max-width:650px){header,main,footer{padding:16px}section,article{padding:18px}th,td{padding:8px}}
@media print{nav,.toolbar,dialog{display:none}body{background:white}section,article{break-inside:avoid}}
'''
SCRIPT = '''
<script>
const zoom=document.getElementById('zoom');
document.querySelectorAll('a.figure-link').forEach(a=>a.addEventListener('click',event=>{
 if(!zoom.showModal)return;event.preventDefault();const img=a.querySelector('img');
 document.getElementById('zoom-image').src=a.href;document.getElementById('zoom-image').alt=img.alt;zoom.showModal();}));
document.getElementById('close').addEventListener('click',()=>zoom.close());
let category='all';const search=document.getElementById('search');
function filter(){const q=search.value.toLowerCase();document.querySelectorAll('article[data-category]').forEach(a=>{
 a.hidden=!(category==='all'||a.dataset.category===category)||!a.textContent.toLowerCase().includes(q);});}
document.querySelectorAll('[data-filter]').forEach(b=>b.addEventListener('click',()=>{category=b.dataset.filter;
 document.querySelectorAll('[data-filter]').forEach(x=>x.classList.toggle('active',x===b));filter();}));
if(search)search.addEventListener('input',filter);
</script>'''


def page(title, body, supplemental=False):
    nav = '<a href="index.html">Main findings</a><a href="supplemental.html">Supplemental figures</a>'
    if not supplemental:
        nav += '<a href="#structure">Structure</a><a href="#interfaces">Interfaces</a><a href="#energies">MM/GBSA</a>'
    return ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f'<title>{escape(title)}</title><style>{STYLE}</style></head><body><header>'
            '<p class="eyebrow">Fc–FcγRIIIa VAL158 · preliminary exploration</p>'
            f'<h1>{escape(title)}</h1><nav aria-label="Report navigation">{nav}</nav></header>'
            f'<main>{body}</main><footer class="note">Offline report · Existing simulations and energies · '
            'No external assets or uploads. Keep both pages and the assets folder together.</footer>'
            '<dialog id="zoom"><button id="close">Close enlarged figure</button>'
            '<div class="scroll"><img id="zoom-image" alt="Enlarged figure"></div></dialog>'
            f'{SCRIPT}</body></html>')


def table(headers, rows):
    body = []
    for row in rows:
        cells = []
        for i, value in enumerate(row):
            tag = 'th scope="row"' if i == 0 else 'td'
            close = 'th' if i == 0 else 'td'
            cells.append(f'<{tag}>{escape(str(value))}</{close}>')
        body.append('<tr>' + ''.join(cells) + '</tr>')
    return ('<div class="table-scroll"><table><thead><tr>' +
            ''.join(f'<th scope="col">{escape(h)}</th>' for h in headers) +
            '</tr></thead><tbody>' + ''.join(body) + '</tbody></table></div>')


def read(path):
    with Path(path).open() as fh:
        return list(csv.DictReader(fh))


def figure_html(fig):
    return (f'<figure><a class="figure-link" href="assets/{escape(fig["png"])}">'
            f'<img loading="lazy" src="assets/{escape(fig["png"])}" alt="{escape(fig["title"])}"></a>'
            f'<figcaption>{escape(fig["reading"])} '
            f'<a href="assets/{escape(fig["pdf"])}">PDF</a></figcaption></figure>')


def prepare_home_content(project, figures):
    """Read summaries and prepare tables and highlighted values for the main page."""
    structural = read(project / 'defined_figures_review/replicate_means.csv')
    binding = read(project / 'defined_binding_figures/comparisons.csv')
    pairs = read(project / 'defined_figures_review/pair_comparisons.csv')
    components = read(project / 'defined_binding_figures/replicates.csv')
    ranges = read(project / 'sensitivity_results/summary/finding_ranges.csv')
    windows = read(project / 'sensitivity_results/summary/binding_window_comparisons.csv')
    systems = [r['system'] for r in binding]
    def mean_sd(system, metric):
        values = [float(r[metric]) for r in structural if r['system'] == system]
        return f'{statistics.mean(values):.2f} ± {statistics.stdev(values):.2f}'
    metrics = ['receptor_common_reference_rmsd_A', 'receptor_self_fit_rmsd_A',
               'partner_COM_distance_A', 'complex_Rg_A']
    structure_table = table(['System', 'Receptor after CH3 fit (Å)', 'Receptor self-fit (Å)',
                             'Partner COM distance (Å)', 'Complex Rg (Å)'],
                            [[s] + [mean_sd(s, m) for m in metrics] for s in systems])
    highlights = [('salt_bridge', 'D265A', 'A:265', 'C:120', 'Fc A Asp265–receptor Lys120'),
                  ('contact', 'A327S', 'A:N297:glycan2', 'C:N162:glycan1', 'Fc A N297 glycan ring 2–receptor N162 NAG'),
                  ('contact', 'R255A', 'B:332', 'C:161', 'Fc B residue 332–receptor residue 161')]
    interface_rows = []
    sensitivity_rows = []
    for kind, system, fc, receptor, label in highlights:
        selected = {r['system']: r for r in pairs if (r['kind'], r['Fc_label'], r['receptor_label']) == (kind, fc, receptor)}
        wt, mutant = selected['WT'], selected[system]
        def percent(row):
            return f'{100*float(row["mean_occupancy"]):.1f} ± {100*float(row["sd_between_replicates"]):.1f}%'
        interface_rows.append([system, label, 'Salt-bridge candidate' if kind == 'salt_bridge' else 'Heavy-atom contact',
                               percent(wt), percent(mutant), f'{100*float(mutant["delta_occupancy_vs_WT"]):+.1f} pp'])
        if system in ('D265A', 'A327S'):
            result = next(r for r in ranges if (r['kind'], r['system'], r['Fc_label'], r['receptor_label']) == (kind, system, fc, receptor))
            sensitivity_rows.append([system, f'{float(result["minimum_delta_pp"]):+.1f} to {float(result["maximum_delta_pp"]):+.1f} pp'])
    interface_table = table(['Variant', 'Interface pair', 'Definition', 'WT occupancy', 'Variant occupancy', 'Change'], interface_rows)
    delta = next(r for r in binding if r['system'] == 'D265A')
    def component(system, name):
        return statistics.mean(float(r[name]) for r in components if r['system'] == system)
    eel = component('D265A', 'EEL') - component('WT', 'EEL')
    egb = component('D265A', 'EGB') - component('WT', 'EGB')
    energy = next(f for f in figures if f['source'].stem == 'binding_comparison')
    e333 = [float(r['ddG_vs_WT']) for r in windows if r['system'] == 'E333A']
    return structure_table, interface_table, sensitivity_rows, delta, eel, egb, energy, e333


def home(project, figures):
    """Render the main findings page from prepared tables and highlighted values."""
    (structure_table, interface_table, sensitivity_rows, delta,
     eel, egb, energy, e333) = prepare_home_content(project, figures)
    body = '<section class="takeaway"><h2>Local interface changes are clearer than mutant energy rankings</h2>'
    body += ('<p class="lead">Most mutants cannot be reliably ordered by whole-partner MM/GBSA with three short replicates. '
             'D265A has the clearest unfavorable shift, while selected protein and glycan contacts change consistently across tested settings.</p>'
             '<p class="study">8 systems · 3 × 10 ns per system · 140 analyzed snapshots per run · stored times 3.16–10.11 ns</p>'
             '<p class="note">The independent sample size is three trajectories per system. No mutant comparison passes '
             'the seven-comparison Holm-adjusted 0.05 threshold. This does not establish equivalence.</p></section>')
    body += ('<section id="structure"><p class="tag">01 · Structure</p><h2>Similar overall dimensions; receptor motion varies</h2>'
             '<p>Mean ± between-replicate SD. The CH3-aligned receptor RMSD includes motion relative to Fc; '
             'self-fitting removes its overall rotation and translation.</p>' + structure_table +
             '<p class="note">Similar COM distance and radius of gyration do not rule out local changes or establish structural equivalence. '
             'A shared WT backbone reference is used. No formal structural group test or equivalence test is reported.</p>'
             '<a href="supplemental.html#figure-1">Structural overview and individual trajectory figures →</a></section>')
    body += ('<section id="interfaces"><p class="tag">02 · Specific interfaces</p><h2>Charged and glycan contacts respond differently</h2>'
             '<p>Two chemically interpretable changes and one large contact gain are highlighted; the full contact map is supplemental. '
             'Occupancies show mean ± between-replicate SD; changes are percentage points (pp).</p>' + interface_table +
             '<p class="note">A/B are Fc chains; C is receptor. Salt-bridge candidates use 4.0 Å, heavy-atom contacts 4.5 Å. '
             'D265A removes the charged Asp group, so its salt-bridge loss is expected. Contact gains can reflect rearrangement '
             'and do not establish stronger binding. These are selected descriptive observations, without contact-level significance tests.</p>'
             '<a href="supplemental.html#figure-5">All residue-pair heatmaps and salt-bridge changes →</a></section>')
    body += ('<section id="energies"><p class="tag">03 · Whole-partner MM/GBSA</p><h2>D265A shifts unfavorably; other rankings remain uncertain</h2>'
             f'<p>D265A − WT: <strong>{float(delta["ddG_vs_WT"]):+.2f} kcal/mol</strong>; unadjusted Welch 95% interval '
             f'{float(delta["ci95_lo"]):+.2f} to {float(delta["ci95_hi"]):+.2f}; Holm-adjusted p = {float(delta["p_holm"]):.4f}. '
             'Every other mutant’s unadjusted interval crosses zero.</p>' + figure_html(energy) +
             f'<p>For D265A, direct electrostatics ({eel:+.2f} kcal/mol) and GB solvation ({egb:+.2f}) largely offset, '
             f'leaving {eel+egb:+.2f} kcal/mol net polar change.</p>'
             '<p class="note">Entropy is omitted. These are whole-partner estimates, not measured affinities, Kd folds, '
             'or an energy decomposition of a particular residue pair.</p>'
             '<a href="supplemental.html#figure-3">Component panels and contact–energy comparisons →</a></section>')
    body += ('<section id="sensitivity"><p class="tag">Analysis sensitivity</p><h2>Key contact directions persist; exact effects depend on settings</h2>' +
             table(['Highlighted occupancy change', 'Range across 18 tested settings'], sensitivity_rows) +
             f'<p>E333A’s mean energetic difference ranges from {min(e333):+.2f} to {max(e333):+.2f} kcal/mol across later windows, '
             'including a sign change. Nested windows add no independent replicates and do not demonstrate convergence.</p>'
             '<a href="supplemental.html#figure-19">Sensitivity figures, settings, and data tables →</a></section>')
    body += ('<section><h2>Scope of this exploration</h2><p class="note">Short runs, fixed protonation, truncated receptor glycans, '
             'and the chosen force field/solvent model limit interpretation. Imaging validation remains incomplete; '
             'legacy BSA is excluded because of invalid values. The original WT minimization reached its cycle limit '
             'without meeting its gradient threshold. See the repository provenance and reproduction guides for details.</p></section>')
    return page('Structure, interfaces, and binding estimates', body)


def supplemental(figures, tablelinks):
    toolbar = ('<div class="toolbar" aria-label="Figure filters">' + ''.join(
        f'<button class="{"active" if key == "all" else ""}" data-filter="{key}">{label}</button>'
        for key, label in [('all', 'All'), ('structure', 'Structure'), ('interface', 'Interfaces'),
                           ('energies', 'MM/GBSA'), ('sensitivity', 'Sensitivity')]) +
        '<input id="search" aria-label="Search supplemental figures" placeholder="Search figures"></div>')
    body = ('<section><h2>Complete figure set</h2><p>All 23 figures, including the main energy comparison, '
            'are retained here for detailed review. Open a figure to enlarge it. Definitions and caveats are expandable.</p>' + toolbar +
            '<details><summary>CSV data tables</summary><ul>' + ''.join(tablelinks) + '</ul></details></section>')
    for i, fig in enumerate(figures, 1):
        body += (f'<article id="figure-{i}" data-category="{escape(fig["category"])}">'
                 f'<p class="tag">{escape(fig["category"])} · Figure {i}</p><h2>{escape(fig["title"])}</h2>'
                 f'<p>{escape(fig["interpretation"])}</p>' + figure_html(fig) +
                 '<details><summary>Question and interpretation limits</summary>'
                 f'<p>{escape(fig["question"])}</p><p>{escape(fig["limitation"])}</p></details></article>')
    return page('Supplemental figures and analysis details', body, supplemental=True)
