# Repository source selection

The publication tree separates supported work, historical exploration, and selected
results without deleting original local files or changing data paths.

| Included path | Role |
|---|---|
| `analysis/` | Supported structural/interface/MMGBSA comparisons, report generation, sensitivity |
| `slurm/` | Mutation tool, workflow automation, replicate submission, Slurm job templates |
| `MD/` | Simulation input templates |
| `preprocessing/` | Protein extraction, complex building, dry topology, input configuration/export |
| `analysis/structure/` | Trajectory geometry, interface contacts, structural plots |
| `analysis/mmgbsa/` | MM/GBSA calculation and binding-energy comparisons |
| `analysis/reporting/` | Offline report generation |
| `config/`, `variants.txt` | Analysis protocol and variant definitions |
| `tests/` | Synthetic workflow, numerical, and report tests |
| `docs/` | Provenance, run history, input manifest, reproduction/release instructions |
| `legacy/` | Historical analysis/plotting/checking scripts, archived for reference |
| `report/index.html` | Main findings: structural summaries, selected contacts, MM/GBSA, brief sensitivity |
| `report/supplemental.html` | Complete 23-figure gallery and detailed definitions |
| `report/assets/` | Selected PNG/PDF/CSV sources for both pages |

## Excluded local material

Original root-level copies of scripts now archived under `legacy/` are ignored.
So are raw structures/topologies/trajectories, all simulation run folders,
`partial_runs/`, scratch inspection scripts, installer snapshots, logs, temporary
source-check exports, old `analysis_review/`, and duplicate/older figure directories.
The new input manifest identifies the actual dataset; export that separately
rather than force-adding raw data to Git.

No files were deleted from the study. The root Git allowlist determines what a
fresh checkout contains. Archival scripts can assume the old working layout and
are not the supported reproduction route. The main README and reproduction guide
point to `analysis` instead.

## Updating the report

Run `python -m analysis.reporting.build_report` into a new output directory and ZIP. Review the
new pages before selecting them for a later release. This avoids overwriting the
current report or earlier analysis results. Captions remain specific to this
exploration; changed datasets require scientific review of the interpretation.

The source repository is published on GitHub. GitHub Pages hosting is not configured;
the HTML works offline. Licensing remains an owner decision.

## Modules and commands

Run all commands from the repository root. Python analysis commands use
`python -m analysis.SUBPACKAGE.COMMAND` so imports resolve consistently without modifying
`PYTHONPATH` or installing the repository. The package initializer has no side effects.

Shared definitions live in `preprocessing/system_definition.py` (glycans, chains,
and disulfides) and `config/analysis_protocol.json` (analysis settings).
`analysis/io.py` supplies common checksum and CSV helpers. These modules play
the shared-definition role that header files serve in C/C++. Numerical behavior
and input/output dataset filenames are unchanged.

| Task | Command |
|---|---|
| Extract protein | `python preprocessing/extract_protein.py` |
| Mutate protein | `python slurm/mutate_protein.py --help` |
| Build complex | `python preprocessing/build_complex.py --help` |
| Prepare variant | `bash slurm/prepare_variant.sh NAME Fc:NUM:AA3` |
| Prepare dry topology | `bash preprocessing/prepare_dry_topology.sh NAME` |
| Submit replicates | `bash slurm/submit_replicates.sh NAME FIRST LAST` |
| Build and submit all variants | `bash slurm/submit_all_variants.sh NREP` |
| MM/GBSA postprocessing | `bash analysis/mmgbsa/run_mmgbsa.sh NAME REP` |
| Analyze trajectory / audit inputs | `python -m analysis.structure.trajectory --help` |
| Compare structures | `python -m analysis.structure.plot_structure --help` |
| Compare binding energies | `python -m analysis.mmgbsa.plot_binding --help` |
| Sensitivity analysis | `python -m analysis.sensitivity --help` |
| Generate report | `python -m analysis.reporting.build_report --help` |
| Run tests | `python -B -m unittest discover -s tests -v` |

Preparation and simulation commands require their licensed tools and a fresh
working copy containing the input dataset. Submission commands submit jobs;
review them and obtain authorization before running on the cluster.

The former `analysis_v2/` commands are now package commands above.
`defined_analysis.py` became `trajectory.py`, `plot_defined.py` became
`plot_structure.py`, and `build_review.py` became `build_report.py`.
Root preparation and simulation scripts moved to the task folders with descriptive
names. Existing generated reports and historical provenance retain their original
source references. Archived `legacy/` scripts retain the old layout assumptions.

## Workflow stages

1. **Preprocessing:** extract the protein, apply mutations, assemble the glycosylated
   complex, prepare dry topologies, and configure paths. `slurm/prepare_variant.sh`
   automates complex mutation/building through `preprocessing/build_complex.py`.
   The standalone `slurm/mutate_protein.py` performs a single-site mutation,
   including the receptor F158V preparation described in the provenance guide.
   These tools are distinct; their scientific mutation/repacking behavior is unchanged.
2. **MD runs:** `slurm/submit_replicates.sh` generates and submits the dependent
   minimization, heating, equilibration, production, and postprocessing jobs using
   the Amber input files in `MD/`. `slurm/submit_all_variants.sh` wraps preparation
   and submission for the variant list.
3. **Postprocessing:** `analysis/structure/trajectory.py` computes structural
   observables and contacts; `analysis/mmgbsa/run_mmgbsa.sh` images/reorders
   trajectories and calculates energies. Their companion plotting modules compare
   replicates. `analysis/sensitivity.py` checks both structure and energy settings;
   `analysis/reporting/` builds the combined report.

Additional supporting pieces are `config/analysis_protocol.json` for scientific
settings, `preprocessing/system_definition.py` for shared molecular definitions,
`preprocessing/configure_protocol.py` for dataset paths, and
`preprocessing/export_inputs.py` plus `slurm/export_inputs.slurm` for dataset
packaging. Tests, provenance documentation, curated report assets, and archived
legacy code are kept separate from executable workflow stages.

The current MM/GBSA command combines imaging, atom reordering, topology preparation,
and energy calculation. Splitting those into independent reusable commands would
be a further refactor; this layout change preserves the existing calculation.
