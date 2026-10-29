#!/bin/bash
# usage: bash check_replicate.sh NAME REP        (run from the prep folder)
# One-shot sanity check of a finished replicate. Prints PASS / FAIL / INFO lines.
NAME=$1; R=${2:-1}
ROOT=$(pwd); D=$ROOT/$NAME/rep$R
START=${START:-301}; INTERVAL=${INTERVAL:-5}
PY=${PY:-python}
NF=0
pass() { echo "[PASS] $*"; }
fail() { echo "[FAIL] $*"; NF=$((NF+1)); }
info() { echo "[INFO] $*"; }
[ -d "$D" ] || { echo "no folder $D"; exit 1; }
cd $D

echo "=== $NAME rep$R ==="
# 1. every stage finished normally (timing block only prints on a clean exit)
for st in min heat eq1 eq2 prod; do
  if [ ! -f $st.out ]; then fail "$st.out missing"
  elif grep -q "Total wall time" $st.out; then pass "$st finished"
  else fail "$st.out has no timing block (crashed or still running?)"; fi
  grep -qE "NaN|\*{6,}|vlimit exceeded" $st.out && fail "$st.out contains NaN / ****** / vlimit"
done

# 2. minimization energy
if [ -f min.out ]; then
  awk '/NSTEP +ENERGY/{getline; if(!f){f=$2; s=$1} l=$2; ls=$1}
       END{printf "[INFO] min: energy %s (step %s) -> %s (step %s)\n", f, s, l, ls}' min.out
fi

# 3. heating reached the target (first NSTEP line at the final step = real value, not the averages)
if [ -f heat.out ]; then
  T=$(awk '/NSTEP =/{ if (match($0,/NSTEP = +[0-9]+/)) {s=substr($0,RSTART+8,RLENGTH-8)+0}
           if (match($0,/TEMP\(K\) = +[0-9.]+/)) {t=substr($0,RSTART+10,RLENGTH-10)+0}
           if (s>m) {m=s; v=t} } END{print v}' heat.out)
  TARGET=$(grep -oE "temp0 *= *[0-9.]+" heat.in | head -1 | grep -oE "[0-9.]+$")
  info "heat final temperature ${T} K (target ${TARGET:-?})"
  awk -v t="$T" -v g="${TARGET:-310.15}" 'BEGIN{exit !(t>g-8 && t<g+8)}' && pass "heat reached target" || fail "heat did not reach target"
fi

# 4. density in the NPT stages: real series excludes the last two values (average, RMS)
for st in eq1 eq2 prod; do
  [ -f $st.out ] || continue
  awk -v st=$st '/Density/{n++; v[n]=$3}
    END{ if(n<4){print "[INFO] " st ": too few density values"; exit}
         lo=hi=v[1]; for(i=1;i<=n-2;i++){if(v[i]<lo)lo=v[i]; if(v[i]>hi)hi=v[i]}
         printf "[INFO] %s density: min %.4f max %.4f, run average %.4f +- %.4f\n", st, lo, hi, v[n-1], v[n]
         if (st=="prod" && (v[n-1]<0.99 || v[n-1]>1.04)) print "[FAIL] prod average density outside 0.99-1.04"
         else if (st=="prod") print "[PASS] prod average density in range" }' $st.out
done | tee /tmp/dens_$$.txt | sed 's/^\[FAIL\]/[FAIL]/'
grep -q "^\[FAIL\]" /tmp/dens_$$.txt && NF=$((NF+1)); rm -f /tmp/dens_$$.txt

# 5. trajectory contents versus ntwprt, nstlim, ntwx
EXPECT_ATOMS=$(cat $ROOT/$NAME.ntwprt 2>/dev/null)
IN_NTWPRT=$(grep -oE "ntwprt *= *[0-9]+" prod.in | grep -oE "[0-9]+$")
[ "$IN_NTWPRT" = "$EXPECT_ATOMS" ] && pass "prod.in ntwprt = $IN_NTWPRT" || fail "prod.in ntwprt ($IN_NTWPRT) != $NAME.ntwprt ($EXPECT_ATOMS)"
if command -v ncdump >/dev/null; then
  A=$(ncdump -h prod.nc | grep -E "^\s+atom =" | grep -oE "[0-9]+")
  FR=$(ncdump -h prod.nc | grep -oE "\(([0-9]+) currently\)" | grep -oE "[0-9]+")
  NS=$(grep -oE "nstlim *= *[0-9]+" prod.in | grep -oE "[0-9]+$"); WX=$(grep -oE "ntwx *= *[0-9]+" prod.in | grep -oE "[0-9]+$")
  [ "$A" = "$EXPECT_ATOMS" ] && pass "prod.nc has $A atoms" || fail "prod.nc atoms $A != expected $EXPECT_ATOMS"
  [ -n "$NS" ] && [ -n "$WX" ] && [ "$FR" = "$((NS/WX))" ] && pass "prod.nc has $FR frames (nstlim/ntwx)" || fail "prod.nc frames $FR vs nstlim/ntwx ${NS:-?}/${WX:-?}"
else info "ncdump not found: trajectory atom/frame check skipped"; FR=""; fi

# 6. MM-GBSA
if [ -f FINAL.dat ]; then
  grep -q "residues: complex 614 receptor 441 ligand 173" post.txt && pass "topology sizes 614/441/173" || fail "topology size line missing in post.txt"
  awk '/Differences/{d=1} d && /^(BOND|ANGLE|DIHED) /{v=$2; if (v<0) v=-v; if (v>0.01) bad=1; n++}
       END{ if(n==3 && !bad) print "[PASS] bonded terms in Differences are ~0 (topologies consistent)"; else print "[FAIL] bonded terms in Differences are not ~0" }' FINAL.dat | tee /tmp/bond_$$.txt
  grep -q "^\[FAIL\]" /tmp/bond_$$.txt && NF=$((NF+1)); rm -f /tmp/bond_$$.txt
  PROC=$(grep -oE "^[0-9]+ frames were processed" post.txt | grep -oE "^[0-9]+")
  if [ -n "$FR" ]; then EXP=$(( (FR-START)/INTERVAL + 1 )); [ "$PROC" = "$EXP" ] && pass "MM-GBSA used $PROC frames (expected $EXP)" || fail "MM-GBSA used ${PROC:-?} frames, expected $EXP"
  else info "MM-GBSA used ${PROC:-?} frames"; fi
  awk '/DELTA TOTAL/{l=$0} END{split(l,a," "); printf "[INFO] DELTA TOTAL mean %s  SD %s  SEM %s kcal/mol\n", a[3], a[4], a[5]}' FINAL.dat
else fail "FINAL.dat missing"; fi

# 7. distance to periodic image (every 50th frame; needs MDAnalysis, skipped if unavailable)
if [ -z "$NOIMG" ] && [ -f prod_dry.nc ] && [ -f $ROOT/analysis/image_dist.py ]; then
  M=$(timeout 600 $PY $ROOT/analysis/image_dist.py dry.prmtop prod_dry.nc 50 2>/dev/null | tail -1 | grep -oE "[0-9.]+ A" | grep -oE "[0-9.]+")
  if [ -z "$M" ]; then info "image-distance check skipped (MDAnalysis or script unavailable)"
  else awk -v m=$M 'BEGIN{exit !(m>12)}' && pass "min distance to periodic image ${M} A (> 12)" || fail "min distance to periodic image ${M} A (< 12)"; fi
fi

echo "=== $NF failure(s) ==="
