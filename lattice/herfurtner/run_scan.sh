#!/bin/bash
cd "$(dirname "$0")"
mkdir -p out
AMAX=${AMAX:-3000}; BMAX=${BMAX:-150}; DMAX=${DMAX:-12000}; N=${N:-30}
NSH=12
while read M J1 J2; do
  [ -z "$M" ] && continue
  CMAX=$(( DMAX / (M*M) )); [ $CMAX -gt 4000 ] && CMAX=4000
  TAG="c_${M}_${J1}_${J2}"
  if [ -s out/${TAG}.txt ]; then echo "$TAG already done"; continue; fi
  W=$(( (AMAX + NSH - 1) / NSH ))
  for i in $(seq 0 $((NSH-1))); do
    LO=$(( 1 + i*W )); HI=$(( LO + W - 1 )); [ $HI -gt $AMAX ] && HI=$AMAX
    ./02_hscan $M $J1 $J2 $LO $HI $BMAX $CMAX $N > out/${TAG}_s$i.txt 2> out/${TAG}_s$i.log &
  done
  wait
  cat out/${TAG}_s*.txt > out/${TAG}.txt
  rm -f out/${TAG}_s*.txt
  echo "$TAG  CMAX=$CMAX  hits=$(wc -l < out/${TAG}.txt)  $(date +%H:%M:%S)"
done < classes.txt
echo SCAN_DONE
