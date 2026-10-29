# Repository source selection

The publication tree separates supported work, historical exploration, and selected
results without deleting original local files or changing data paths.

| Included path | Role |
|---|---|
| `analysis/` | Supported structural/interface/MMGBSA comparisons, report generation, sensitivity, tests |
| `slurm/` | Path configuration, Amber setup, dataset exporter, Slurm templates |
| `MD/` | Simulation input templates |
| `preparation/` | Protein/complex preparation and shared system definitions |
| `simulation/` | Simulation submission and MM/GBSA commands |
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

Run `python -m analysis.build_report` into a new output directory and ZIP. Review the
new pages before selecting them for a later release. This avoids overwriting the
current report or earlier analysis results. Captions remain specific to this
exploration; changed datasets require scientific review of the interpretation.

The source repository is published on GitHub. GitHub Pages hosting is not configured;
the HTML works offline. Licensing remains an owner decision.

## Modules and commands

Run all commands from the repository root. Python analysis commands use
`python -m analysis.COMMAND` so imports resolve consistently without modifying
`PYTHONPATH` or installing the repository. The package initializer has no side effects.

Shared definitions live in `preparation/system_definition.py` (glycans, chains,
and disulfides) and `config/analysis_protocol.json` (analysis settings).
`analysis/io.py` supplies common checksum and CSV helpers. These modules play
the shared-definition role that header files serve in C/C++. Numerical behavior
and input/output dataset filenames are unchanged.

| Task | Command |
|---|---|
| Extract protein | `python preparation/extract_protein.py` |
| Mutate protein | `python preparation/mutate_protein.py --help` |
| Build complex | `python preparation/build_complex.py --help` |
| Prepare variant | `bash preparation/prepare_variant.sh NAME Fc:NUM:AA3` |
| Prepare dry topology | `bash preparation/prepare_dry_topology.sh NAME` |
| Submit replicates | `bash simulation/submit_replicates.sh NAME FIRST LAST` |
| Build and submit all variants | `bash simulation/submit_all_variants.sh NREP` |
| MM/GBSA postprocessing | `bash simulation/run_mmgbsa.sh NAME REP` |
| Analyze trajectory / audit inputs | `python -m analysis.trajectory --help` |
| Compare structures | `python -m analysis.plot_structure --help` |
| Compare binding energies | `python -m analysis.plot_binding --help` |
| Sensitivity analysis | `python -m analysis.sensitivity --help` |
| Generate report | `python -m analysis.build_report --help` |
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
