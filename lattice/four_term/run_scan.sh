#!/bin/bash
# usage: ./run_scan.sh TAG AMAX CMAX DMAX FMAX GMAX N NJOBS [classfile]
set -e
cd "$(dirname "$0")"
TAG=$1; AMAX=$2; CMAX=$3; DMAX=$4; FMAX=$5; GMAX=$6; N=$7; NJOBS=${8:-12}; CF=${9:-classes.txt}
mkdir -p out
while read -r RN RD M J1 J2; do
  [ -z "$RN" ] && continue
  NAME="${TAG}_r${RN}o${RD}_M${M}_${J1}_${J2}"
  for ((s=0;s<NJOBS;s++)); do
    ./03_fscan $RN $RD $M $J1 $J2 1 $AMAX $CMAX $DMAX $FMAX $GMAX $N $NJOBS $s \
        > out/${NAME}.s${s}.txt 2> out/${NAME}.s${s}.log &
  done
  wait
  cat out/${NAME}.s*.txt > out/${NAME}.txt; rm -f out/${NAME}.s*.txt
  grep -h '^# deep' out/${NAME}.s*.log | awk -v n="$NAME" '{d+=$3;h+=$5} END{printf "%-28s deep %-14d hits %d\n", n, d, h}'
  rm -f out/${NAME}.s*.log
done < "$CF"
