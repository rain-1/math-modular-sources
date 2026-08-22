#!/bin/bash
cd "$(dirname "$0")"
CH=11; LO=1; HI=3000; W=$(( (HI-LO+1+CH-1)/CH ))
for i in $(seq 0 $((CH-1))); do
  A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
  ./01_class_scan -1 $A $B 60 300 400 26 > out/fb_$i.txt 2> out/fb_$i.log &
done
wait
cat out/fb_*.txt | awk "$3!=0" > out/freebig.txt
echo "freebig done $(wc -l < out/freebig.txt)"
