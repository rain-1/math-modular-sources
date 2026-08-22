#!/bin/bash
cd "$(dirname "$0")"
read i N lo hi <<< "$1"
printf 'n=scanlevel(NLEV,DIVS,TS,FS,%s,%s,"out/mod/j%s.txt");\n' "$lo" "$hi" "$i" > out/mod/d$i.gp
gp -q 03_modscan.gp 03b_exact.gp lev/L$N.gp out/mod/d$i.gp > out/mod/j$i.log 2>&1
rm -f out/mod/d$i.gp
