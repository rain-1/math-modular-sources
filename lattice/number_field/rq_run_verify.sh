#!/bin/bash
SD=/home/ubuntu/code/math-modular-sources/lattice/number_field
W=$(mktemp /tmp/rqw_XXXX.gp)
{ printf 'INFILE="%s";\n' "$1"; cat $SD/rq_verify.gp; } > $W
/usr/bin/gp -q $W > "$2" 2>&1
rm -f $W
