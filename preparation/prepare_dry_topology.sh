#!/bin/bash
# usage: bash preparation/prepare_dry_topology.sh NAME
# Makes dry.NAME.prmtop and NAME.ntwprt (solute atom count for ntwprt). Use for WT too.
set -e
NAME=$1
source "$(dirname "${BASH_SOURCE[0]}")/../slurm/amber_env.sh"
cpptraj > dry_$NAME.log << EOT
parm $NAME.prmtop
parmstrip :WAT,Na+,Cl-
parmwrite out dry.$NAME.prmtop
EOT
N=$(cpptraj << EOT | grep -oE "contains [0-9]+ atoms" | grep -oE "[0-9]+"
parm dry.$NAME.prmtop
parminfo
EOT
)
echo $N > $NAME.ntwprt
echo "$NAME: dry solute atoms = $N  (written to $NAME.ntwprt)"
