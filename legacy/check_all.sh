#!/bin/bash
# usage: bash check_all.sh           (run from the prep folder)
#        NOIMG=1 bash check_all.sh   skip the slow image-distance test
# Runs check_replicate.sh on every NAME/repN folder and prints one summary line each.
ROOT=$(pwd); HERE=$(cd "$(dirname "$0")" && pwd)
[ -f "$HERE/check_replicate.sh" ] || { echo "check_replicate.sh must be next to check_all.sh"; exit 1; }
if command -v sacct >/dev/null; then
  echo "== failed / cancelled / timed-out Slurm jobs in the last 7 days"
  sacct -u $USER --starttime=$(date -d '7 days ago' +%F) --format=JobID,JobName%24,State,ExitCode,Elapsed -X 2>/dev/null \
    | grep -E "FAILED|CANCELLED|TIMEOUT|OUT_OF_ME|NODE_FAIL" || echo "none"
  echo
fi
printf "%-16s %-11s %6s %10s %10s\n" replicate status fails dG density
nok=0; nfail=0; ninc=0; declare -A reps
for d in */rep*/; do
  NAME=${d%%/*}; R=${d##*rep}; R=${R%/}
  reps[$NAME]=$(( ${reps[$NAME]:-0} + 1 ))
  if [ ! -f $NAME/rep$R/prod.out ] || ! grep -q "Total wall time" $NAME/rep$R/prod.out; then
    printf "%-16s %-11s\n" "$NAME rep$R" INCOMPLETE; ninc=$((ninc+1)); continue
  fi
  out=$(bash "$HERE/check_replicate.sh" $NAME $R 2>&1)
  if ! echo "$out" | grep -q "failure(s)"; then
    printf "%-16s %-11s %s\n" "$NAME rep$R" SCRIPT_ERR "$(echo "$out" | head -1)"; nfail=$((nfail+1)); continue
  fi
  nf=$(echo "$out" | grep -c '^\[FAIL\]')
  dg=$(echo "$out" | grep -oE "DELTA TOTAL mean -?[0-9.]+" | grep -oE "\-?[0-9.]+$")
  de=$(echo "$out" | grep "prod density" | grep -oE "run average [0-9.]+" | grep -oE "[0-9.]+$")
  if [ -z "$dg" ]; then st=NO_MMGBSA; elif [ "$nf" = 0 ]; then st=OK; else st=FAIL; fi
  printf "%-16s %-11s %6s %10s %10s\n" "$NAME rep$R" $st $nf "${dg:--}" "${de:--}"
  echo "$out" | grep '^\[FAIL\]' | sed 's/^/      /'
  case $st in OK) nok=$((nok+1));; FAIL|NO_MMGBSA) nfail=$((nfail+1));; esac
done
echo
echo "OK: $nok   problems: $nfail   incomplete: $ninc"
for n in "${!reps[@]}"; do [ "${reps[$n]}" -lt 3 ] && echo "note: $n has ${reps[$n]} replicate folder(s) (expected 3)"; done | sort
