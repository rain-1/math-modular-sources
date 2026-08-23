#!/usr/bin/env python3
"""emit a gp driver running jtest5 on every CANDIDATE row of an analysis json."""
import sys, json
A = json.load(open(sys.argv[1]))
rows = [r for r in A if r.get('verdict','').startswith('CANDIDATE')] if isinstance(A[0].get('verdict',None),str) else A
with open(sys.argv[2],'w') as fh:
    fh.write('\\r 06_jtest.gp\ndefault(parisizemax, 8000000000);\n')
    for r in rows:
        rn,rd,m,j1,j2 = r['cls']; a,c,d,f,C = r['row']
        lab = f"c{rn}_{rd}_{m}_{j1}_{j2}__{a}_{c}_{d}_{f}_{C}"
        fh.write(f'print("{lab} ", jtest5({rn},{rd},{m},{j1},{j2}, {a},{c},{d},{f},{C}));\n')
    fh.write('quit;\n')
print(len(rows), "rows")
