#!/bin/bash
# usage: bash simulation/submit_all_variants.sh [NREP]      builds every variant in variants.txt, submits NREP replicates each
NREP=${1:-3}
grep -vE '^\s*(#|$)' variants.txt | while read NAME SPEC CLASS; do
  echo "=== $NAME ($SPEC, $CLASS)"
  bash preparation/prepare_variant.sh $NAME $SPEC || { echo "BUILD FAILED for $NAME, skipping"; continue; }
  bash simulation/submit_replicates.sh $NAME 1 $NREP
done
