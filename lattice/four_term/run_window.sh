#!/bin/bash
# positive-score window scan: |d| <= 2a+5, |C| <= (a+2)/M^2+1  (complete for |lambda_2|<1)
# usage: ./run_window.sh TAG AMAX CMAX FMAX N NJOBS
cd "$(dirname "$0")"; mkdir -p out
TAG=$1; AMAX=$2; CMAX=$3; FMAX=$4; N=$5; NJOBS=${6:-8}
while read -r RN RD M J1 J2; do
  [ -z "$RN" ] && continue
  NAME="${TAG}_r${RN}o${RD}_M${M}_${J1}_${J2}"
  for ((s=0;s<NJOBS;s++)); do
    ./03_fscan $RN $RD $M $J1 $J2 1 $AMAX $CMAX -1 $FMAX -1 $N $NJOBS $s \
      > out/${NAME}.s${s}.txt 2>/dev/null &
  done
  wait
  cat out/${NAME}.s*.txt > out/${NAME}.txt; rm -f out/${NAME}.s*.txt
  echo "$NAME $(wc -l < out/${NAME}.txt)"
done < classes.txt
