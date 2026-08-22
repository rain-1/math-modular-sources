#!/bin/bash
cd "$(dirname "$0")"
for E in 0 3 4 6; do
  CH=6; LO=1; HI=3000; W=$(( (HI-LO+1+CH-1)/CH ))
  for i in $(seq 0 $((CH-1))); do
    A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
    ./01_class_scan $E $A $B 150 150 300 26 > out/e${E}_$i.txt 2> out/e${E}_$i.log &
  done
  wait
  cat out/e${E}_*.txt > out/eclass_$E.txt
  echo "e=$E done $(wc -l < out/eclass_$E.txt)"
done
