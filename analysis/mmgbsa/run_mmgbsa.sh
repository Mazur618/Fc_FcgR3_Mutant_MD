#!/bin/bash
# usage: bash analysis/mmgbsa/run_mmgbsa.sh NAME REP   (submitted automatically after prod)
# Image -> reorder atoms -> MM-GBSA topologies -> MMPBSA.py. Output: NAME/repN/FINAL.dat
set -e
NAME=$1; R=$2
ROOT=$(pwd); D=$ROOT/$NAME/rep$R
START=${START:-301}; INTERVAL=${INTERVAL:-5}   # stored frame 301 = 3.16 ns; same rule for every variant
source "$ROOT/slurm/amber_env.sh"
cp dry.$NAME.prmtop $D/dry.prmtop
cd $D
rm -f mm.dry.prmtop receptor.prmtop ligand.prmtop complex.prmtop cplx_nobox.prmtop

cpptraj << 'EOT'
parm dry.prmtop
trajin prod.nc 1 last 1
autoimage anchor :1-212
trajout prod_dry.nc netcdf
run
EOT
cpptraj << 'EOT'
parm dry.prmtop
trajin prod_dry.nc
fixatomorder outprefix mm
trajout prod_mm.nc netcdf
run
EOT
[ -f mm.dry.prmtop ] || { echo "mm.dry.prmtop not written"; exit 1; }

ante-MMPBSA.py -p mm.dry.prmtop -r receptor.prmtop -l ligand.prmtop --radii mbondi2 -m ":1-441"
cpptraj << 'EOT'
parm mm.dry.prmtop
parmbox nobox
parmwrite out cplx_nobox.prmtop
EOT
parmed -p cplx_nobox.prmtop << 'EOT'
changeRadii mbondi2
outparm complex.prmtop
quit
EOT

nres() { cpptraj << EOT | grep -oE "[0-9]+ residues" | head -1 | grep -oE "[0-9]+"
parm $1
parminfo
EOT
}
C=$(nres complex.prmtop); RC=$(nres receptor.prmtop); L=$(nres ligand.prmtop)
echo "residues: complex $C receptor $RC ligand $L (expect 614 441 173)"
[ "$C" = 614 ] && [ "$RC" = 441 ] && [ "$L" = 173 ] || { echo "size check failed"; exit 1; }
cpptraj << 'EOT' | grep -q UYB || { echo "layout check failed: residue 213 is not a glycan"; exit 1; }
parm mm.dry.prmtop
resinfo :213
EOT

cat > mmgbsa.in << EOT
MM-GBSA single trajectory
&general
  startframe=$START, endframe=999999, interval=$INTERVAL,
  verbose=2, keep_files=0,
/
&gb
  igb=5, saltcon=0.150,
/
EOT
MMPBSA.py -O \
  -i mmgbsa.in \
  -o FINAL.dat \
  -eo perframe.csv \
  -cp complex.prmtop \
  -rp receptor.prmtop \
  -lp ligand.prmtop \
  -y prod_mm.nc
echo "done: $D/FINAL.dat"
