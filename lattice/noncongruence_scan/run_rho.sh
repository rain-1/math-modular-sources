#!/bin/bash
cd "$(dirname "$0")"
RHOS="1/1 2/1 3/1 1/2 3/2 5/2 1/3 2/3 4/3 5/3 7/3 8/3 1/4 3/4 5/4 7/4 9/4 11/4 \
1/5 2/5 3/5 4/5 6/5 7/5 8/5 9/5 11/5 12/5 13/5 14/5 1/6 5/6 7/6 11/6 13/6 17/6"
for R in $RHOS; do
  P=${R%/*}; Q=${R#*/}
  CH=11; LO=1; HI=3000; W=$(( (HI-LO+1+CH-1)/CH ))
  for i in $(seq 0 $((CH-1))); do
    A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
    ./01_class_scan -2 $A $B 60 300 400 26 $P $Q > out/r_${P}_${Q}_$i.txt 2> out/r_${P}_${Q}_$i.log &
  done
  wait
  cat out/r_${P}_${Q}_*.txt | awk '$3!=0' > out/rho_${P}_${Q}.txt
  rm -f out/r_${P}_${Q}_*.txt
  echo "rho=$P/$Q done $(wc -l < out/rho_${P}_${Q}.txt)"
done
echo RHO_DONE
