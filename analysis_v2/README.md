# Defined structural and interface analyses

This is a separate analysis workflow for existing Fc–FcγRIIIa VAL158 trajectories.
It leaves production inputs, Slurm scripts, and previous results unchanged.
The earlier scripts remain available for comparison; they have not been corrected
or replaced by this implementation.

## Scientific questions

1. Do mutations change the sampled complex geometry or internal deformation?
2. Which Fc–receptor contacts change, including contacts involving glycans?
3. Do the contact changes help interpret whole-partner MM/GBSA differences?

Similar global observables do not establish structural equivalence or unchanged
local contacts. This workflow does not measure an electrostatic potential or
perform energy decomposition. Distance-based contacts are structural observations,
not residue interaction energies.

## Fixed protocol

`protocol.json` defines all systems, replicate identities, the reference, window,
cutoffs, and expected topology layout. Paths are relative to the protocol file.
Every comparison requires the complete configured input set; missing runs cannot
silently change the population being compared. For a deliberately smaller study,
make a new protocol listing its systems/replicates, always including WT.

The primary analysis uses one-based trajectory frames 301–1000 with stride 5:
140 snapshots, ending at frame 996. This matches existing MM/GBSA sampling.
Times come from NetCDF via MDAnalysis, converted from ps to ns; frame 1 is not
artificially placed at zero. For denser contact sampling or window sensitivity,
copy the protocol, change its window/stride, and use new output directories.
The primary contact analysis is sampled occupancy, not a contact lifetime estimate.

### Residue identity and validation

The protein mapping is read from `WT.pdb`, preserving original chain/residue
numbers and insertion codes. Glycan labels distinguish the A/B N297 branches
and the receptor N162/N45 stubs. Each run must have the expected 614-residue
original dry topology. Protein residue identities are checked against WT with
the expected mutation applied to both Fc chains; glycan identities are checked
individually. Reordered MM/GBSA topologies must not be supplied here.

`WT/rep1` links to the original root-level WT files. The inventory resolves these
links and rejects repeated trajectory paths. It does not hash/read trajectories.
Run metadata records protocol, topology, and reference SHA256 hashes, source
paths, package versions, analyzed frame count, and first/last sampled times.

### Structural observables

All backbone metrics use N, CA, and C atoms and the same WT protein reference.
They differ from the older CA-only or per-run-reference observables and should
not be appended to older tables as though they were identical measurements.

* Fit both Fc CH3 domains together; measure CH3, CH2, and receptor backbone RMSD
  in that frame. CH3 masks retain the existing domain boundaries.
* Fit Fc and receptor separately to measure internal deformation relative to WT.
* Measure mass-weighted COM distance between complete glycosylated partners and
  mass-weighted radius of gyration of the dry complex.

The current implementation does not yet add domain orientation, glycan RMSD,
RMSF, or a validated surface-area observable. Those remain separate extensions.

### Interface observables

* **Heavy-atom contact:** at least one inter-partner atom pair within 4.5 Å;
  atoms named H* are excluded. This naming assumption matches these Amber
  topologies and must be revisited for other hydrogen naming conventions.
* **Salt bridge candidate:** at least one distance within 4.0 Å between
  ARG NE/NH1/NH2, LYS NZ, or HIP ND1/NE2 and ASP OD1/OD2 or GLU OE1/OE2 on the
  opposite partner. This uses fixed topology protonation and excludes terminal
  charges. It is a geometric criterion, with no hydrogen-bond angle requirement.
* Each residue pair contributes at most once per sampled frame. Occupancy is
  observed frames divided by all analyzed frames, separately for every replicate.
* Distances are measured within the already imaged complex without applying
  minimum-image wrapping. Inputs must be the whole, correctly assembled
  `prod_dry.nc` trajectories. Existing periodic-image diagnostics remain necessary;
  this implementation does not certify imaging or absence of image interactions.

Comparative tables use the union of observed pairs across all completed runs.
Missing pairs in a completed run count as zero. Heatmaps display pairs reaching
20% mean occupancy in any system, including WT. This threshold controls display
only; all observed pairs remain in CSV tables. The heatmaps paginate all qualifying
pairs and use a fixed ±100 percentage-point color scale. Replicates are weighted
equally. Uncertainty bars show between-replicate SD, not confidence intervals.

## Usage

From the repository root (the former `prep` folder), audit paths and mapping (lightweight, no trajectory reads):

```bash
python analysis_v2/defined_analysis.py inventory --out analysis_v2/input_inventory.json
python -B -m unittest discover -s analysis_v2 -p 'test_*.py' -v
```

The trajectory command requires NumPy and MDAnalysis in an existing analysis
environment. Run trajectory processing in a scheduler allocation, not on a shared login node.
Limit numerical libraries to one thread; see `scripts/analyze.slurm`.
The Python command itself does not submit jobs:

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python analysis_v2/defined_analysis.py run --system WT --replicate 1 --out defined_results/WT/rep1
```

Repeat for each configured system/replicate. The output directory must not
already exist. Each successful run writes `frames.csv`, `residue_map.csv`,
`pairs.csv`, and finally `complete.json`. Partial runs are rejected by plotting;
after a failure choose a new output directory rather than reusing partial output.

Once all runs are complete, produce tables and PNG/PDF figures (also requires
Matplotlib):

```bash
python analysis_v2/plot_defined.py defined_results/*/rep* --out defined_figures
```

No dependencies are installed automatically. Outputs refuse existing directories
and files. Plot generation is separate from trajectory processing.

### Existing MM/GBSA comparisons

`plot_binding.py` reads the existing per-frame energy CSVs and requires SciPy,
NumPy, and Matplotlib. It performs no new energy calculations. It writes
replicate/variant tables, a replicate-energy and Welch-interval comparison,
separate direct-electrostatic/GB and net-component panels, and exploratory
contact-count versus energy plots. Component bars show descriptive means.

```bash
python analysis_v2/plot_binding.py defined_results/*/rep* --out defined_binding_figures
```

Snapshot ordinals and counts are checked against structural results. Energy CSVs
contain snapshot ordinals rather than original trajectory frame IDs: matching
the original MM/GBSA window remains a provenance assumption. All current files
have 140 snapshots. Confidence intervals are unadjusted Welch 95% intervals;
the table separately provides Holm-adjusted p-values across mutant comparisons.
Stored times in these trajectories include equilibration: selected frames run
from 3.16 to 10.11 ns, not 3.00 to 9.95 ns on the stored clock.

## Sensitivity checks

`sensitivity.py` reads frames 301–1000 once per trajectory at the saved 10 ps
cadence. It reduces atom distances to minimum residue-pair distances and evaluates
all combinations of start frame (301/501/701), stride (1/5), and distance cutoff
(contacts 4.0/4.5/5.0 Å; salt-bridge candidates 3.5/4.0/4.5 Å). These give 18
conditions per interaction type and 36 per run, with 700/500/300 dense frames or
140/100/60 coarse frames. Later windows are nested, not independent replicates.

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python analysis_v2/sensitivity.py defined_results/*/rep* --out sensitivity_results
```

Each run must reproduce the primary pair counts exactly at the original settings.
New per-run CSVs and completion markers retain results for every tested condition.
Summary tables include all observed pairs with explicit zeros in completed runs,
and plots show cutoff, window, sampling, and highlighted-finding sensitivity.
The sampling comparison figure quantifies cadence effects at the primary window
and cutoffs. The finding-range table covers the full 18-condition grid.

Structural and MM/GBSA window comparisons subset existing 50 ps observations;
no new energies or dense structural RMSDs are calculated. Energy snapshot ordinals
inherit the primary window provenance assumption. The plots are descriptive:
they neither test equivalence nor demonstrate convergence. Window changes do not
resolve the limitations of short runs or fixed starting structures.

## Annotated offline image viewer

```bash
python analysis_v2/build_review.py --out analysis_review --zip analysis_review.zip
```

The report has two pages: `index.html` condenses structural summaries, selected
protein/glycan interfaces, and MM/GBSA, with a short sensitivity note.
`supplemental.html` retains all 18 primary and five sensitivity figures, detailed
captions, filtering, image enlargement, PDFs, and CSV tables. Image files are
linked through relative paths rather than embedded repeatedly. Keep both pages
and the `assets/` directory together, or extract the generated ZIP.

The reviewed current build is in `report/`; the original local `analysis_review/`
remains unchanged. No server, external scripts, fonts, or uploads are used.
These pages include unpublished results; building them does not publish them.

Outputs require new paths. To preserve a revised report, choose a new directory
and ZIP filename rather than overwriting an existing report.

## Remaining work and interpretation limits

* Numerical tests check known rotations/translations, rejection of reflections,
  and capped contacts against brute-force distances. They require NumPy and
  MDAnalysis. Mapping, counting, inventory, windows, and failure handling have
  separate lightweight tests. Compare real-run observables with cpptraj as well.
* Add angle-qualified hydrogen bonds over the same window; salt-bridge candidates
  are not a substitute for hydrogen-bond analysis.
* Investigate invalid legacy surface values before including BSA; distinguish
  molecular surface from solvent-accessible area and document the area convention.
* Add decomposition on existing snapshots as separate compute work. Existing
  MM/GBSA values represent whole-partner energy differences, not interface-only
  energies; no Kd fold prediction should be inferred from them.
* Record fixed protonation, neutralizing-only explicit ions versus 0.150 M GB salt,
  truncated receptor glycans, omitted entropy, and limited simulation sampling.
* Contact-cutoff/window sensitivity checks are complete in the local exploration;
  inspect imaging independently. Similarity
  within these three short replicates does not demonstrate absence of an effect.
* Environment specifications and repository documentation are provided at the root.
  Select a license and arrange data access before publication. Generated inventory
  JSON contains local absolute paths and should not be published unchanged.
