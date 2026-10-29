#!/bin/bash
# usage: bash slurm/prepare_variant.sh D265A Fc:265:ALA
# Run in the PyRosetta env, in the folder with preprocessing/build_complex.py, WT.pdb, 7URU.cif.
set -e
NAME=$1; SPEC=$2
IFS=: read CH NUM AA3 <<< "$SPEC"
declare -A O2T=( [A]=ALA [R]=ARG [N]=ASN [D]=ASP [C]=CYS [E]=GLU [Q]=GLN [G]=GLY [H]=HIS [I]=ILE [L]=LEU [K]=LYS [M]=MET [F]=PHE [P]=PRO [S]=SER [T]=THR [W]=TRP [Y]=TYR [V]=VAL )
EXPECT=${O2T[${NAME:0:1}]}
for c in A B; do
  FOUND=$(awk -v n=$NUM -v c=$c '$1=="ATOM" && $3=="CA" && $5==c && $6==n {print $4}' WT.pdb)
  if [ "$FOUND" != "$EXPECT" ]; then
    echo "ERROR: chain $c residue $NUM is '$FOUND' in WT.pdb, expected $EXPECT from name $NAME"; exit 1
  fi
done
python preprocessing/build_complex.py --name $NAME --mut $SPEC | tee prep_$NAME.txt
source "$(dirname "${BASH_SOURCE[0]}")/amber_env.sh"
tleap -f build_$NAME.in > build_$NAME.out
tail -3 build_$NAME.out
grep -q "Errors = 0" build_$NAME.out || { echo "tleap reported errors: see build_$NAME.out"; exit 1; }
grep -E "ions required|Placed" build_$NAME.out | head
ls -l $NAME.prmtop $NAME.rst7
bash preprocessing/prepare_dry_topology.sh $NAME
