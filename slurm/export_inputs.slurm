#!/usr/bin/env bash
#SBATCH --job-name=fc-input-export
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=01:00:00
# Review partition/account/storage requirements. Submit from the repository root.
# Activate a Python 3 environment first. No numerical packages are required.
set -euo pipefail
: "${SOURCE_ROOT:?Set SOURCE_ROOT to the original Fc_FcR3 directory}"
: "${DATASET_OUT:?Set DATASET_OUT to a new directory whose parent exists}"
"${PY:-python}" slurm/export_inputs.py --source-root "$SOURCE_ROOT" \
  --out "$DATASET_OUT" --groups analysis preparation simulation
