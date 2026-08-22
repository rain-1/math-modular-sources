#!/bin/bash
# R3 class: al in [2,3000] even, |de|<=150, |ga|<=150, |ze|<=300, N=26
cd "$(dirname "$0")"
CH=12; LO=2; HI=3000
W=$(( (HI-LO+1+CH-1)/CH ))
for i in $(seq 0 $((CH-1))); do
  A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
  ./01_class_scan 0 $A $B 150 150 300 26 > out/r3_$i.txt 2> out/r3_$i.log &
done
wait
cat out/r3_*.txt > out/r3_all.txt
wc -l out/r3_all.txt
