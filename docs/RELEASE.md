# Release review

The local repository selects source/documentation and the curated two-page `report/`. No remote, push, public
upload, or license grant is configured. Review `git status --short` and
`git diff --cached` before the first commit. Git author identity belongs to the
repository owner and is not configured by setup.

Before public release:

- Approve the README and its distinction between legacy and current analysis.
- Choose a license with rights holders; document separate input/data rights.
- Cite the experimental basis for mutation classes and the structure, force
  fields, software, and analysis methods.
- Decide how unpublished inputs and selected outputs will be shared. Supply a
  dataset identifier/access instructions. Local input checksums and recovered
  historical provenance are compiled in `input-manifest.json` and `INPUT_PROVENANCE.md`;
  trajectory hashes will be produced by the dataset exporter.
- Validate imaging and investigate BSA before using them in conclusions. These
  extensions need not block a clearly labeled preliminary code release.
- Run the pinned numerical tests and full workflow against the released dataset
  in a compute allocation. Synthetic tests do not establish end-to-end reproduction.
- Review site configuration and output-writing behavior before submitting jobs.

The root Git ignore file is an allowlist. New source names need an explicit entry.
Raw data, unselected results, local inventories, archives, and scratch scripts
remain untouched on disk. Only the explicitly included report assets are selected for presentation; other
generated tables/figures are not automatically selected for publication. Record a source commit and input checksums with archived results.
