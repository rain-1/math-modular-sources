#!/bin/bash
cd "$(dirname "$0")"
mkdir -p out
AMAX=${AMAX:-4000}; BMAX=${BMAX:-200}; DMAX=${DMAX:-20000}; N=${N:-30}
NSH=6
while read M J1 J2; do
  [ -z "$M" ] && continue
  CMAX=$(( DMAX / (M*M) )); [ $CMAX -gt 4000 ] && CMAX=4000
  TAG="c_${M}_${J1}_${J2}"
  W=$(( (AMAX + NSH - 1) / NSH ))
  for i in $(seq 0 $((NSH-1))); do
    LO=$(( 1 + i*W )); HI=$(( LO + W - 1 )); [ $HI -gt $AMAX ] && HI=$AMAX
    ./02_hscan $M $J1 $J2 $LO $HI $BMAX $CMAX $N > out/${TAG}_$i.txt 2> out/${TAG}_$i.log &
  done
  wait
  cat out/${TAG}_*.txt > out/${TAG}.txt
  rm -f out/${TAG}_[0-9].txt
  echo "$TAG  CMAX=$CMAX  hits=$(wc -l < out/${TAG}.txt)"
done < classes.txt
echo SCAN_DONE
