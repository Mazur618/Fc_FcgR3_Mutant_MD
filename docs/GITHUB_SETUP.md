# GitHub handoff

The publication selection contains the supported workflow, archived legacy code,
provenance/reproduction guides, and the two-page results report. Raw trajectories,
structures, topologies, simulation outputs, and local scratch files are excluded.
The repository owner must choose a license before granting reuse rights; absence
of a license does not prevent creating a repository but does not grant a license.

## Existing local repository

The original `prep` folder already has Git initialized on `main`. Check the
selection with `git status --short` and `git diff --cached --stat`.
Set your own author identity locally if it has not already been configured:

```bash
git config --local user.name "YOUR NAME"
git config --local user.email "YOUR GITHUB EMAIL OR NOREPLY EMAIL"
git commit -m "Prepare preliminary Fc-FcGR3 exploration and report"
```

These identity settings are not authentication credentials. Use GitHub's normal
SSH or HTTPS authentication; never place tokens/passwords in source files or chat.

## Clean folder or ZIP

The handoff export is generated from the Git index, not by copying the study
folder. It contains no `.git` metadata or simulation data. The ZIP keeps hidden
files such as `.github/workflows/tests.yml`. Extract the complete folder to keep
the report's relative links intact. To start Git from that exported folder:

```bash
git init -b main
git config --local user.name "YOUR NAME"
git config --local user.email "YOUR GITHUB EMAIL OR NOREPLY EMAIL"
git add .
git commit -m "Prepare preliminary Fc-FcGR3 exploration and report"
```

## Connect to GitHub when ready

Create an empty repository in your own GitHub account or organization. If pushing
this prepared history, leave GitHub's README/license/gitignore initialization
options unchecked. Substitute its repository URL below, then run these commands
from the local Git repository you intend to publish:

```bash
git remote add origin REPOSITORY_URL
git push -u origin main
```

No remote or push is configured by preparation. The report is at
`report/index.html`, with all detailed figures at `report/supplemental.html`.
The pages work offline. Creating the repository does not configure website hosting.

## Validation status and scope

The last local clean-export check passed 16 tests with four numerical tests
skipped because the active Python lacks the numerical dependencies. CI installs
the specified dependencies before running the suite. The full trajectory workflow
has not been rerun during repository preparation. Input provenance, export tools,
and scientific limits are documented separately in this folder.
