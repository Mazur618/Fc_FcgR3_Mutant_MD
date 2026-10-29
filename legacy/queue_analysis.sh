#!/bin/bash
# usage: bash queue_analysis.sh            queue analysis for every finished replicate
#        FORCE=1 bash queue_analysis.sh    redo ones that already have summary.txt
# Run from the prep folder.
ROOT=$(pwd)
for d in */rep*/; do
  NAME=${d%%/*}; R=${d##*rep}; R=${R%/}
  [ -f $NAME/rep$R/prod_dry.nc ] || { echo "skip $NAME rep$R (no prod_dry.nc yet)"; continue; }
  if [ -f $NAME/rep$R/summary.txt ] && [ -z "$FORCE" ]; then echo "skip $NAME rep$R (done)"; continue; fi
  JID=$(sbatch --parsable --job-name=${NAME}_r${R}_an --nodes=1 --ntasks-per-node=1 --mem=16G --time=02:00:00 \
    -o $ROOT/$NAME/rep$R/analyze.txt -e $ROOT/$NAME/rep$R/analyze.err \
    --wrap="cd $ROOT; bash analyze_variant.sh $NAME $R")
  echo "$NAME rep$R -> job $JID"
done
