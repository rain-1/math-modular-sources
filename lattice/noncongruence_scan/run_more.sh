#!/bin/bash
cd "$(dirname "$0")"
go(){ CLS=$1; LO=$2; HI=$3; DE=$4; GA=$5; ZE=$6; TAG=$7; CH=$8
  W=$(( (HI-LO+1+CH-1)/CH ))
  for i in $(seq 0 $((CH-1))); do
    A=$((LO+i*W)); B=$((A+W-1)); if [ $B -gt $HI ]; then B=$HI; fi
    ./01_class_scan $CLS $A $B $DE $GA $ZE 26 > out/m_${TAG}_$i.txt 2> out/m_${TAG}_$i.log &
  done
  wait
  cat out/m_${TAG}_*.txt > out/more_${TAG}.txt
  echo "$TAG done $(wc -l < out/more_${TAG}.txt)"
}
go 2 1 3000 400 150 300 e2deep 11
go -1 1 500 70 60 120 freewide 11
echo MORE_DONE
