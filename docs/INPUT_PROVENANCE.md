# Recovered input provenance

This audit used the original `Fc_FcR3` directory and its `prep` child. Original
files and results were left unchanged. Sources below are relative to `Fc_FcR3`,
not to a specific university filesystem. This account distinguishes evidence in
current files from historical steps that cannot be established exactly.

Historical filenames below describe the original study layout. Current commands use
`preprocessing/extract_protein.py`, `slurm/mutate_protein.py`, and
`preprocessing/build_complex.py`; see [the layout guide](REPOSITORY_LAYOUT.md).

## Authoritative inputs and provenance chain

1. `prep/7URU.cif` and `7URU.cif` are byte-identical. The source structure has
   receptor author-chain C residue 158 PHE. The mutation is **F158V**, not P158V.
2. `prep/make_p.py` extracts protein, removes hydrogens, alternative conformations,
   ligands and waters, and writes `protein_only.pdb`. The existing product has
   A/B/C counts 212/213/171 and retains C158 PHE.
3. `prep/WT.pdb` is the authoritative prepared protein reference: original chain
   labels/numbering, the same 596 protein residues, and C158 VAL. `prep/mutate.py`
   supports this operation by chain/number, with configurable repacking radius.
   The parent `scripts/pyrosetta_mutate.py` is another candidate historical tool:
   it takes a pose index and repacks the entire protein. The exact original
   invocation, selected tool/radius, and PyRosetta seed are not recorded.
   Pose index 579 maps to C158 in the retained 596-residue ordering.
   A documented reconstruction route (not the recovered historical invocation),
   run only in a fresh directory with licensed preparation dependencies, is:

   ```bash
   python preprocessing/extract_protein.py
   python slurm/mutate_protein.py protein_only.pdb --chain C --num 158 --aa VAL --radius 8 -o WT.pdb
   ```

   Radius 8 is an explicit reconstruction choice. Repacking is stochastic; use the
   supplied hashed `prep/WT.pdb` for exact analysis/reference reproduction.
4. `prep/build.py` assembles protein with two eight-residue Fc N297 glycans and
   NAG stubs at receptor N162/N45, declares six disulfides and glycosidic bonds,
   then writes a GLYCAM PDB and LEaP input. The prepared reference, original CIF,
   current generator, and actual generated simulation inputs are retained locally.
5. Parent and `prep` versions of `WT.prmtop` and `WT.rst7` are byte-identical.
   Their current topology contains VAL at original dry residue index 579 and
   includes CMAP sections. Their SHA256 hashes appear in the input manifest.
6. WT replicate 1 reuses the original parent run. Its dry trajectory is
   `prod1_dry.nc`, original-order topology is `dry.WT.prmtop`, and matching energy
   CSV is `WT_perframe_s301.csv`. The three `prep/WT/rep1` links resolve to these.
   Other replicates use their own files in `prep/<system>/rep<N>/`.
7. Parent/prep `dry.WT.prmtop` whole-file hashes differ, but every named `%FLAG`
   section has identical content. Preserve each exact source hash rather than
   silently replacing the historical file with the other copy.
8. `mmgbsa_s301.in`, its submission script, and the header in `FINAL_WT_s301.dat`
   establish the primary WT window: start frame 301, interval 5, `igb=5`,
   salt 0.150 M, Fc mask `:1-441`, receptor mask `:442-614`, using `prod1_mm.nc`.
   The original-order dry topology must not be replaced by MM-reordered inputs.

## Historical files that must not be substituted

- Parent `WT.pdb` has blank chain labels, renumbered residues, and **PHE**, not
  VAL, at the receptor-158 equivalent position. Its protein sequence differs from
  `prep/WT.pdb` at that site (after normalizing histidine/Cys naming). Their CA
  coordinates match, which does not make the files interchangeable.
- Parent `tleap.in` and `leap.log` describe early protein-only attempts, including
  failures and a nonglycosylated successful attempt. They do not reconstruct the
  current glycosylated system.
- `prep/build_WT.out` is stale: it reports ff14SB, while current `build_WT.in`
  specifies ff19SB and current WT topology/logs include CMAP. Mutant build logs
  explicitly report ff19SB. The precise historical WT build invocation cannot be
  proven from this old output; use the hashed actual topology for analysis.
- `WT_perframe.csv`, `FINAL_WT.dat`, and `mmgbsa.in` belong to an earlier WT
  calculation starting at frame 157. Use the `_s301` versions for this study.
- `prod_1_imaged.nc`, `snap.nc`, protein-only test topologies, and assorted old
  `*_old`/nobox/MM files are not replacements for the manifest's analysis inputs.

## Inventory and run history

`input-manifest.json` lists **113 existing files**, grouped as:

- `analysis`: 73 files — prepared WT reference plus three inputs for each of
  24 runs (original dry topology, imaged dry trajectory, matching energy CSV).
- `preparation`: one file — original `7URU.cif` glycan/structure source.
- `simulation`: 39 files — actual solvated topology, coordinates, solute atom
  count, assembled GLYCAM PDB for each system, plus seven mutant protein PDBs.

Logical destination paths are independent of original source paths. All sources
are local and present. Non-trajectory inputs have SHA256 checksums. The 24 dry
trajectories total **2,845,123,872 bytes**; only their metadata was inspected.
Trajectory hashes are deliberately deferred to the export copy in an allocation.
The manifest contains no absolute cluster paths and no coordinates.

`run-history.csv` records source filenames, seeds recovered from bounded log
headers, completion status, and minimization counts for all 120 stages. All have
normal timing completion records and all 120 headers record a seed (the
minimization seed is not evidence of stochastic minimization). All 24
minimizations completed 5,000 cycles. For original WT, the final RMS gradient
was 0.091507 versus requested `drms=0.00010`: the run reached its iteration cap,
rather than meeting the gradient criterion. This observation alone does not
establish that the subsequent equilibrated production ensemble is invalid.
Fc mutant preparation logs identify PyRosetta 2026.39 (Python 3.10, release
49e05e797ca4383c6038f090e5db6c2593cb4963). This does not establish the environment
used for the original receptor F158V preparation.
Logged seeds help describe historical runs but do not ensure identical GPU
trajectories, particularly without original environments and packing seeds.

## Compile a self-contained dataset

`preprocessing/export_inputs.py` materializes the manifest into a NEW directory, follows
no unresolved parent-folder links in the output, checks known hashes, and computes
all output hashes while copying (including trajectories). A successful export ends
with `dataset-manifest.json`. Failure leaves a partial directory for inspection;
choose a new destination for the next attempt.

For the analysis-only dataset, in a compute allocation:

```bash
python preprocessing/export_inputs.py --source-root /path/to/Fc_FcR3 \
  --out /path/to/new-dataset --groups analysis preparation
```

To include actual prepared simulation inputs, add `simulation` to `--groups`.
`slurm/export_inputs.slurm` provides the full-copy job template. Review it and
set exported `SOURCE_ROOT`/`DATASET_OUT` before any authorized submission.
No export, trajectory hashing, or job submission was run on the login node.

Configure the analysis against the resulting dataset root and its `WT.pdb` using
`preprocessing/configure_protocol.py`. Remaining historical uncertainty is the precise
protein mutation/repacking invocation and full software environment, not missing
analysis files. A public dataset location, licensing, and experimental-label
references still need owner decisions.
