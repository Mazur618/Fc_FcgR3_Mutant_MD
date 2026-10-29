#!/bin/bash
# usage: bash run_all.sh [NREP]      builds every variant in variants.txt, submits NREP replicates each
NREP=${1:-3}
grep -vE '^\s*(#|$)' variants.txt | while read NAME SPEC CLASS; do
  echo "=== $NAME ($SPEC, $CLASS)"
  bash build_variant.sh $NAME $SPEC || { echo "BUILD FAILED for $NAME, skipping"; continue; }
  bash submit_variant.sh $NAME 1 $NREP
done
