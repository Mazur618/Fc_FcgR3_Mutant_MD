#!/bin/bash
# usage: bash analyze_variant.sh NAME REP
# Run from the prep folder, after the replicate's post step has made prod_dry.nc.
# Writes into NAME/repN/: rmsd*.dat rg.dat comdist.dat sasa.dat rmsf.dat pucker.dat
#   hb_interface_avg.dat  stability.png  mm.png  summary.txt
NAME=$1; R=$2
ROOT=$(pwd); D=$ROOT/$NAME/rep$R
PY=${PY:-python}
START=${START:-301}; INTERVAL=${INTERVAL:-5}     # must match post_variant.sh
[ -d "$D" ] || { echo "no folder $D"; exit 1; }
cd $D
for f in prod_dry.nc dry.prmtop; do
  [ -f $f ] || { echo "$f missing in $D (post step not finished?)"; exit 1; }
done
cp $ROOT/analysis/analysis.cpptraj $ROOT/analysis/pucker.cpptraj .

source "$ROOT/scripts/amber_env.sh"
cpptraj -i analysis.cpptraj > analysis.log 2>&1 || echo "analysis.cpptraj failed: see $D/analysis.log"
cpptraj -i pucker.cpptraj   > pucker.log   2>&1 || echo "pucker.cpptraj failed: see $D/pucker.log"

# interface H-bonds only: one side Fc(+glycans, residues 1-425/597-612), other CD16a(+stubs)
if [ -f hb_all_avg.dat ]; then
awk 'NR==1{print; next}
{
  split($1,a,/[_@]/); split($2,d,/[_@]/);
  ra=a[2]+0; rd=d[2]+0;
  fa=(ra<=425 || (ra>=597 && ra<=612)); fd=(rd<=425 || (rd>=597 && rd<=612));
  if (fa != fd) print
}' hb_all_avg.dat > hb_interface_avg.dat
fi

{
  echo "== $NAME rep$R"
  $PY $ROOT/analysis/stability_plot.py . --mm perframe.csv --start $START --step $INTERVAL -o stability.png
  echo; echo "== MM-GBSA DELTA TOTAL"
  [ -f perframe.csv ] && $PY $ROOT/plot_mm.py . --start $START --step $INTERVAL
  echo; echo "== pucker"
  [ -f pucker.dat ] && $PY $ROOT/analysis/pucker_summary.py
  echo; echo "== distance to periodic image (every 10th frame)"
  $PY $ROOT/analysis/image_dist.py dry.prmtop prod_dry.nc 10 | tail -1
  echo; echo "== top interface H-bonds"
  [ -f hb_interface_avg.dat ] && head -12 hb_interface_avg.dat
} > summary.txt 2>&1
echo "done: $D/summary.txt  $D/stability.png  $D/mm.png"
