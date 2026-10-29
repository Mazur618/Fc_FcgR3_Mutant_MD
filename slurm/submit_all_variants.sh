#!/bin/bash
# usage: bash slurm/submit_all_variants.sh [NREP]      builds every variant in variants.txt, submits NREP replicates each
NREP=${1:-3}
grep -vE '^\s*(#|$)' variants.txt | while read NAME SPEC CLASS; do
  echo "=== $NAME ($SPEC, $CLASS)"
  bash slurm/prepare_variant.sh $NAME $SPEC || { echo "BUILD FAILED for $NAME, skipping"; continue; }
  bash slurm/submit_replicates.sh $NAME 1 $NREP
done
