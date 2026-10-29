#!/bin/bash
# usage: bash simulation/submit_replicates.sh NAME REP_FIRST REP_LAST      e.g. bash simulation/submit_replicates.sh D265A 1 3
# Needs MD/<stage>.in for each stage (copy your WT inputs there).
# Each replicate: min -> heat -> eq1 -> eq2 -> prod (dependent GPU jobs) -> post (CPU job).
NAME=$1; R1=${2:-1}; R2=${3:-1}
STAGES=${STAGES:-"min heat eq1 eq2 prod"}
ROOT=$(pwd)
NTW=$(cat $NAME.ntwprt)
PARTITION=${PARTITION:-gpuq}
GPU_RESOURCE=${GPU_RESOURCE:-gpu:1}
for R in $(seq $R1 $R2); do
  D=$ROOT/$NAME/rep$R
  mkdir -p $D
  PREV=$ROOT/$NAME.rst7
  DEP=""
  for st in $STAGES; do
    [ -f MD/$st.in ] || { echo "missing MD/$st.in"; exit 1; }
    grep -qE "ntwprt *= *[0-9]+" MD/$st.in || [ "$st" != "prod" ] || echo "WARNING: prod.in has no ntwprt line"
    sed -E "s/ntwprt *= *[0-9]+/ntwprt=$NTW/" MD/$st.in > $D/$st.in
    REF=""; grep -qE "ntr *= *1" $D/$st.in && REF="-ref $PREV"
    TRAJ=""; grep -qE "ntwx *= *[1-9]" $D/$st.in && TRAJ="-x $D/$st.nc"
    if [ "$st" = "prod" ]; then T=${TIME_PROD:-08:00:00}; else T=${TIME_OTHER:-02:00:00}; fi
    cat > $D/$st.sbatch << SB
#!/bin/bash
#SBATCH --job-name=${NAME}_r${R}_${st}
#SBATCH --partition=$PARTITION
#SBATCH --gres=$GPU_RESOURCE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=$T
#SBATCH -o $D/$st.txt
#SBATCH -e $D/$st.err
#SBATCH --mail-type=FAIL

source "$ROOT/slurm/amber_env.sh"
cd $D

pmemd.cuda -O \\
  -i $D/$st.in \\
  -p $ROOT/$NAME.prmtop \\
  -c $PREV $REF \\
  -o $D/$st.out \\
  -r $D/$st.rst7 \\
  -inf $D/$st.mdinfo $TRAJ
SB
    JID=$(sbatch --parsable ${DEP:+--dependency=afterok:$DEP} $D/$st.sbatch)
    echo "$NAME rep$R $st -> job $JID"
    DEP=$JID; PREV=$D/$st.rst7
  done
  cat > $D/post.sbatch << SB
#!/bin/bash
#SBATCH --job-name=${NAME}_r${R}_mm
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH -o $D/post.txt
#SBATCH -e $D/post.err
#SBATCH --mail-type=FAIL

cd $ROOT
bash simulation/run_mmgbsa.sh $NAME $R
SB
  JID=$(sbatch --parsable --dependency=afterok:$DEP $D/post.sbatch)
  echo "$NAME rep$R post -> job $JID"
done
