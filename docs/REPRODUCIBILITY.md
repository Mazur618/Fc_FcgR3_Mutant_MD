# Reproduction guide

Run commands from this repository root, regardless of its directory name.
Use new output paths for every attempt. Existing analysis directories, including
partial runs, are rejected. Preserve them and choose a new path after a failure.

## Input contract

Supply the original WT protein reference `WT.pdb`, preserving chain A/B/C residue
numbering, and these files for all eight systems and replicates 1, 2, and 3:

```text
data-root/
  WT/rep1/dry.prmtop
  WT/rep1/prod_dry.nc
  WT/rep1/perframe.csv
  ... WT/rep2, WT/rep3, D265A/rep1, ... E333A/rep3
```

The topology must use the original 614-residue order (425 Fc protein residues,
171 receptor protein residues, 18 sugars), not the reordered MM/GBSA topology.
The trajectory must contain the complete, correctly assembled/imaged complex
and at least 1,000 saved frames. The energy CSV must contain 140 MM/GBSA Delta
Energy rows corresponding to original frames 301, 306, ..., 996. CSV ordinals
alone cannot certify this window identity.

Historical local WT replicate 1 uses links outside this repository. A portable
dataset must supply the actual files or valid links within its own data tree.
The repository setup does not copy these inputs. Production inputs alone cannot
reconstruct the original stochastic trajectories. The [input manifest](input-manifest.json) locates every required file, with hashes
for non-trajectory inputs. See [recovered provenance](INPUT_PROVENANCE.md) and the
export script to materialize the dataset and hash trajectories in an allocation.
An external archive/access location remains a release decision.

## Configure paths

Keep `config/analysis_protocol.json` as the scientific template. Write a separate
configured protocol, whose data/reference paths are relative to its own location:

```bash
mkdir -p work
python slurm/configure_protocol.py \
  --data-root /path/to/data-root --reference /path/to/WT.pdb \
  --out work/protocol.json
python -m analysis.trajectory --protocol work/protocol.json \
  inventory --out work/inventory.json
```

Inventory is a lightweight mapping/path audit, not trajectory validation. Inspect
its existence flags; writing an inventory does not require every run to exist.
Generated inventories contain resolved local paths and stay out of Git.
For a deliberately smaller study, use a separate template listing its population,
always including WT. Do not silently omit missing runs from a comparison.

## Compute allocation

Activate the analysis environment and set the numerical thread limits shown in
the main README. Review `slurm/analyze.slurm`, adding your cluster's account
and partition. It requests one CPU per task. Set exported `PROTOCOL` to the
absolute configured protocol path and `OUTPUT_ROOT` to a new absolute output path.
The default array 0–23 follows protocol system/replicate order; adjust its range
if deliberately using a different population. Submit from the repository root,
only after reviewing the script and obtaining required authorization.
No setup command submits jobs.

For one run inside an allocation:

```bash
python -m analysis.trajectory --protocol work/protocol.json run \
  --system WT --replicate 1 --out work/defined_results/WT/rep1
```

Repeat for all configured runs. Successful directories contain `complete.json`,
with source paths, reference/topology/protocol hashes, versions, and frame/time
counts. Whole trajectories are not hashed by the analysis workflow.

After all runs complete, use a suitable allocation for:

```bash
python -m analysis.plot_structure work/defined_results/*/rep* \
  --out work/defined_figures_review
python -m analysis.plot_binding work/defined_results/*/rep* \
  --out work/defined_binding_figures
python -m analysis.sensitivity work/defined_results/*/rep* \
  --out work/sensitivity_results
python -m analysis.build_report --project work \
  --out work/analysis_review --zip work/analysis_review.zip
```

Report `--project` names the directory containing the generated figure/result
folders, without an assumed `prep/` parent. The report has study-specific captions
and expects this study's figure set; reassess captions for changed data/settings.
Generated metadata can contain absolute paths: review before sharing.

Historical completion records contain absolute input paths. Moving only their
output folders does not make energy/sensitivity reanalysis portable. Regenerate
runs against the configured data tree; preserve original provenance unchanged.

The generated `index.html` is the concise home page; `supplemental.html` retains
the full gallery. Copy the whole output folder, including assets, for offline use.
The committed `report/` is a selected presentation artifact; regenerate into a NEW
folder, review differences, and deliberately update a future release rather than
overwriting existing results in place.

## Preparation and new simulations

Preparation additionally needs Gemmi, licensed PyRosetta, AmberTools (original
postprocessing used version 25), and compatible Amber `pmemd.cuda` (original
submission referenced version 24). These are separate from the analysis-only
environment; no recovered preparation environment lock is supplied.
Supply the original prepared VAL158 `WT.pdb` and glycan source `7URU.cif`.
Their preparation history, edits, protonation, and checksums must accompany a
released dataset. `preparation/build_complex.py` does not create the receptor allotype by itself.

Load tools through your site environment, or export `AMBER_ENV` as the path to an
Amber setup script. Leave it unset to use tools on PATH. GPU jobs need
`pmemd.cuda`; postprocessing needs `cpptraj`, `ante-MMPBSA.py`, `parmed`, and
`MMPBSA.py`. `PY` overrides the active Python. `PARTITION` and `GPU_RESOURCE`
override submission defaults (`gpuq`, `gpu:1`). Review resource directives for
your site and run preparation in an allocation.

Build/postprocessing/submission scripts can overwrite generated outputs.
Use a fresh working directory with copies of tracked source and supplied inputs;
do not rerun them over existing results. `simulation/submit_all_variants.sh`, `simulation/submit_replicates.sh`, and
`legacy/queue_analysis.sh` submit jobs and require explicit review/authorization.
They have not been run during repository setup.

Exact regeneration of stochastic trajectories is not supported: Amber seeds are now recovered in `docs/run-history.csv`, but
PyRosetta packing seeds and complete historical environments remain unknown. The reproducible target is the defined analysis of supplied
trajectories, with fresh provenance and new outputs.
