#!/bin/bash
# usage: mk.sh p smax  -> emits a gp driver on stdout
p=$1; s=$2
for r in G C D; do
  echo "report(towerrun(\"$r\",$p,[1,2,3],$s,60));"
done
echo quit
