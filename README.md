# Fc–FcγRIIIa mutation exploration

Preliminary molecular dynamics exploration of how Fc mutations affect the sampled
structure, interface contacts, and whole-partner MM/GBSA energies of a glycosylated
Fc–FcγRIIIa (CD16a, VAL158) complex.

Three 10 ns trajectories per system do not reliably resolve the energetic ordering
of most mutants. D265A shows a consistent unfavorable MM/GBSA shift, but none of
the seven mutant comparisons passes a Holm-adjusted 0.05 significance threshold.
Some local contact changes persist across tested cutoffs and sampling windows.
Nonsignificance does not establish equivalence or absence of a mutation effect.

## Results report

Open [the main findings](report/index.html) for structural comparisons, selected
protein/glycan interface changes, and the MM/GBSA comparison. The
[supplemental page](report/supplemental.html) contains all 23 figures, detailed
time courses, contact heatmaps, sensitivity checks, and CSV tables. Both pages
work offline; keep the `report/assets` folder with them.

## Systems and simulation design

| System | Expected binding class recorded in this project |
|---|---|
| WT | Reference |
| D265A, A327S | Reduced |
| S337A, R255A | WT-like |
| S298A, K334A, E333A | Improved |

Mutations are applied to both Fc heavy chains (A/B); chain C is the receptor.
These class labels are prior expectations, not conclusions from these simulations.
Their experimental sources and assay conditions still need to be documented.

Preparation uses `WT.pdb` as the protein reference and `7URU.cif` as the glycan
source. Receptor C158 is changed from PHE to VAL (F158V). PyRosetta mutates and repacks nearby protein side chains (default 8 Å).
Amber preparation uses ff19SB, GLYCAM06j-1, OPC3 water, a truncated octahedral box
with an 8 Å buffer, and neutralizing ions. Each Fc N297 branch contains eight
sugars; receptor N45/N162 glycans are single-NAG stubs.

Each system has three replicates: one minimization stage (maximum 5,000 cycles,
switching from steepest descent to conjugate gradient after 2,500 cycles), 50 ps
heating, 50 ps restrained equilibration, 50 ps free equilibration, and 10 ns
production. Minimization may stop early on convergence; actual completed cycles
are recorded in `docs/run-history.csv`: all 24 reached 5,000 cycles.
The original WT final gradient remained above its requested threshold. Production uses a 2 fs timestep, SHAKE on
hydrogen bonds, 310.15 K Langevin thermostatting, 1 atm pressure coupling, an 8 Å
nonbonded cutoff, and 10 ps saved frames. Random seeds use Amber `ig=-1`; the actual seeds are recovered in
`docs/run-history.csv`. Exact regeneration of stochastic trajectories is not promised.

## Preliminary observations

Primary analysis uses original one-based frames 301–1000 with stride 5:
140 snapshots per replicate, ending at frame 996. Stored times are 3.16–10.11 ns,
including the equilibration clock. Replicates receive equal weight.

* WT mean whole-partner MM/GBSA is −49.07 kcal/mol; D265A is −39.05 kcal/mol.
  D265A minus WT is +10.02 kcal/mol (unadjusted Welch 95% interval +5.30 to
  +14.74; Holm-adjusted p = 0.0645 across seven comparisons).
* Every other mutant's unadjusted difference interval crosses zero; detailed
  energetic rankings are poorly resolved.
* WT Fc A Asp265–receptor Lys120 salt-bridge candidate occupancy is 97.4%, and
  disappears in D265A. This follows removal of the charged Asp group; it does
  not assign an energy to that residue pair.
* Fc A N297 glycan residue 2–receptor N162 NAG heavy-atom contact occupancy
  decreases from 99.5% in WT to 7.9% in A327S at the primary settings.
* These two contact changes retain their direction across tested windows,
  cadences, and cutoffs. Their magnitudes depend on the contact definition.
  E333A's energetic difference changes sign across analysis windows.

These numbers summarize existing local results. Raw inputs and generated results
are excluded from Git. A source-only clone cannot independently reproduce these
findings until the separately supplied input dataset is available. All required
inputs have been located locally; see the [provenance and export guide](docs/INPUT_PROVENANCE.md)
and [checksummed input manifest](docs/input-manifest.json).

## Setup and tests

Use the repository root as the working directory. Python 3.10, NumPy 2.2.6, and
MDAnalysis 2.9.0 match recorded original structural-analysis metadata. Other pinned
versions are reproducibility targets, not recovered historical versions.
On a machine where dependency installation is permitted:

```bash
python3.10 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-analysis.txt
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -B -m unittest discover -s tests -p 'test_*.py' -v
```

Alternatively, use `conda env create -f environment.yml` from the repository root.
Tests use synthetic inputs and need no trajectories, reference structure, Amber,
or PyRosetta. Numerical tests skip when dependencies are absent; CI installs the
analysis dependencies to exercise them. Scripts never install packages automatically.

See [the reproduction guide](docs/REPRODUCIBILITY.md) for required inputs, portable
path configuration, scheduler use, and figure regeneration. See
[analysis definitions](analysis/README.md) for exact observables. Run trajectory
processing, sensitivity calculations, and simulations in scheduler allocations,
rather than on shared login nodes.

## Repository layout

| Path | Purpose |
|---|---|
| `analysis/` | Current structural/contact analysis, comparisons, sensitivity, offline report |
| `tests/` | Synthetic workflow, numerical, and report tests |
| `config/analysis_protocol.json` | Scientific settings and expected topology layout |
| `slurm/` | Portable protocol writer, Amber environment helper, Slurm analysis template |
| `MD/` | Amber minimization, heating, equilibration, production inputs |
| `preparation/` | Protein extraction/mutation, glycosylated complex construction, dry topology |
| `simulation/` | Replicate submission and MM/GBSA postprocessing |
| `legacy/` | Archived exploration; definitions differ from `analysis`; reference copies only |
| `report/` | Selected two-page results report with figure/PDF/CSV assets |
| `docs/` | Reproduction and release notes |

The Git allowlist includes selected source files. Raw data, generated structures, run directories, and unselected results stay
local. The curated `report/` is explicitly included; other reports are excluded. Review contents and
intended availability before force-adding any ignored files.

## Interpretation limits

* The independent sample size is three trajectories per system, not 140 snapshots.
  Short trajectories and a shared starting structure limit conformational sampling.
  Nested-window sensitivity neither proves convergence nor adds independent runs.
* MM/GBSA uses a single-trajectory protocol, GB `igb=5`, mbondi2 radii, and
  0.150 M implicit salt. It omits entropy and measures whole-partner energy
  differences, not interface-only or per-residue energies. Do not interpret the
  legacy implied Kd folds as quantitative affinity predictions. EEL and GB
  solvation strongly compensate.
* Fixed protonation, truncated receptor glycans, neutralizing-only explicit ions,
  and the GLYCAM/OPC3 combination are modeling limitations not tested here.
* Contact/salt-bridge occupancies are geometric observations, not lifetimes,
  interaction energies, or multiplicity-adjusted contact significance tests.
* Similar global size/separation does not establish structural equivalence.
  Domain orientation, consistent window-matched hydrogen bonds, and residue energy
  decomposition remain extensions.
* Legacy `molsurf` outputs are labeled SASA and include impossible inferred BSA
  values. Exclude legacy BSA and its energy correlations from main conclusions
  until the observable is validated.
* Full imaging validation remains outstanding. Sampled image-distance diagnostics
  do not certify all analyzed frames.
* Energy CSVs have snapshot ordinals, not original trajectory frame IDs; their
  window alignment with structural observations retains a provenance assumption.

## Release status

This is a source scaffold for review, not an archived dataset release. See
[release notes](docs/RELEASE.md) and [repository layout](docs/REPOSITORY_LAYOUT.md). A license must be chosen with the appropriate
rights holders; no license is assigned here. Amber and PyRosetta have separate
access/licensing requirements and are not redistributed. The source and curated report are published on GitHub; the raw dataset is separate.

For repository creation and upload instructions, see [GitHub handoff](docs/GITHUB_SETUP.md).
