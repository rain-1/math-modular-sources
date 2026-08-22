#!/bin/bash
cd "$(dirname "$0")"
RHOS="1/6 1/4 1/3 1/2 2/3 3/4 5/6 1/1 7/6 5/4 4/3 3/2 5/3 7/4 11/6 2/1 7/3 5/2 8/3 3/1"
for R in $RHOS; do
  P=${R%/*}; Q=${R#*/}
  CH=11; LO=1; HI=3000; W=$(( (HI-LO+1+CH-1)/CH ))
  for i in $(seq 0 $((CH-1))); do
    A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
    ./01_class_scan -2 $A $B 60 300 400 26 $P $Q > out/s_${P}_${Q}_$i.txt 2> out/s_${P}_${Q}_$i.log &
  done
  wait
  cat out/s_${P}_${Q}_*.txt | awk '$3!=0' > out/rho2_${P}_${Q}.txt
  rm -f out/s_${P}_${Q}_*.txt out/s_${P}_${Q}_*.log
  echo "rho=$P/$Q done $(wc -l < out/rho2_${P}_${Q}.txt)"
done
echo RHO2_DONE
