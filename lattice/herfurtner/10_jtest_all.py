#!/usr/bin/env python3
"""Emit a gp driver running jtest on every non-degenerate row of out/rows_full.json."""
import json
R=json.load(open('out/rows_full.json'))
rows=[r for r in R if not r['casdeg'] and not r['dbl']]
with open('11_jall.gp','w') as f:
    f.write('read("05_jtest.gp");\n')
    f.write('tst(M,j1,j2,A,B,C)={my(r=jtest(M,j1,j2,A,B,C));\n')
    f.write(' if(r==0, printf("%d %d %d %d %d %d NO 0 0\\n",M,j1,j2,A,B,C),\n')
    f.write('          printf("%d %d %d %d %d %d YES %d %d\\n",M,j1,j2,A,B,C,r[1],r[2]));}\n')
    for r in rows:
        f.write("tst(%d,%d,%d,%d,%d,%d);\n"%(r['M'],r['j1'],r['j2'],r['A'],r['B'],r['C']))
    f.write("quit;\n")
print("rows to test:",len(rows))
