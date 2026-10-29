# Repository source selection

The publication tree separates supported work, historical exploration, and selected
results without deleting original local files or changing data paths.

| Included path | Role |
|---|---|
| `analysis_v2/` | Supported structural/interface/MMGBSA comparisons, report generation, sensitivity, tests |
| `scripts/` | Path configuration, Amber setup, dataset exporter, Slurm templates |
| `templates/` | Simulation input templates |
| Root preparation/submission scripts, `variants.txt` | Preparation and simulation setup |
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
point to `analysis_v2` instead.

## Updating the report

Run `analysis_v2/build_review.py` into a new output directory and ZIP. Review the
new pages before selecting them for a later release. This avoids overwriting the
current report or earlier analysis results. Captions remain specific to this
exploration; changed datasets require scientific review of the interpretation.

No GitHub Pages deployment, remote repository, push, or public upload is configured.
The HTML is ready for local review; hosting and licensing remain owner decisions.
