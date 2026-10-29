#!/usr/bin/env bash
#SBATCH --job-name=fc-interface
#SBATCH --array=0-23
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=04:00:00
# Review resource requests and add your cluster's partition/account before submission.
# Activate the analysis environment first, then submit from the repository root.
# Required exported variables: PROTOCOL (absolute path), OUTPUT_ROOT (new directory).
set -euo pipefail
: "${PROTOCOL:?Set PROTOCOL to the configured protocol JSON}"
: "${OUTPUT_ROOT:?Set OUTPUT_ROOT to a new analysis directory}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-${TMPDIR:-/tmp}}/fc-mpl-${SLURM_JOB_ID}"
PY=${PY:-python}
# Array order follows the protocol, rather than hardcoding system names.
read -r SYSTEM REP < <("$PY" - "$PROTOCOL" "$SLURM_ARRAY_TASK_ID" <<'PYCODE'
import json, sys
cfg = json.load(open(sys.argv[1]))
runs = [(s, r) for s in cfg['systems'] for r in cfg['replicates']]
print(*runs[int(sys.argv[2])])
PYCODE
)
"$PY" analysis_v2/defined_analysis.py --protocol "$PROTOCOL" run \
  --system "$SYSTEM" --replicate "$REP" --out "$OUTPUT_ROOT/$SYSTEM/rep$REP"
