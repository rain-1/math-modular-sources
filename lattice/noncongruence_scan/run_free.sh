#!/bin/bash
cd "$(dirname "$0")"
CH=4; LO=1; HI=300; W=$(( (HI-LO+1+CH-1)/CH ))
for i in $(seq 0 $((CH-1))); do
  A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
  ./01_class_scan -1 $A $B 60 30 60 24 > out/free_$i.txt 2> out/free_$i.log &
done
wait
cat out/free_*.txt > out/freeclass.txt
echo "free done $(wc -l < out/freeclass.txt)"
