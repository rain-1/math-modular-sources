#!/bin/bash
cd ~/ncs
run(){ E=$1; AL=$2; AH=$3; DE=$4; GA=$5; ZE=$6; TAG=$7
  CH=8; W=$(( (AH-AL+1+CH-1)/CH ))
  for i in $(seq 0 $((CH-1))); do
    A=$((AL+i*W)); B=$((A+W-1)); if [ $B -gt $AH ]; then B=$AH; fi
    ./01_class_scan $E $A $B $DE $GA $ZE 26 > w_${TAG}_$i.txt 2> w_${TAG}_$i.log &
  done
  wait
  cat w_${TAG}_*.txt | awk '$3!=0' > wide_${TAG}.txt
  echo "$TAG done $(wc -l < wide_${TAG}.txt)"
}
run 2 1 1000 100 1000 200  e2_gamma
run 2 1 1000 100 100  2000 e2_zeta
run 2 1 30000 60 60   120  e2_alpha
run 0 1 1000 100 1000 200  e0_gamma
run 0 1 1000 100 100  2000 e0_zeta
run 0 1 30000 60 60   120  e0_alpha
echo ALLDONE
