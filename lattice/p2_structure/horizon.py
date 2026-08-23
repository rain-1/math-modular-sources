#!/usr/bin/env python3
"""lattice/p2_structure/horizon.py -- the rational-surrogate horizon n_f(E).

Requires data/vectors_n120.txt (true G) and data/vec_E<E>.txt for each E, produced by

    echo 'EEXP=320;' > pre.gp
    cat lattice/positivity/rows_pos.gp lattice/p2_structure/p2core.gp \
        lattice/p2_structure/p2run.gp pre.gp lattice/p2_structure/run_surrogate.gp > run.gp
    gp -q run.gp

(the dumps are ~1.8 MB each and are not kept in the repository; the summary they
produce is data/horizon.csv).  Prints the first n at which v_1 differs from the
true-G one, per k, against the prediction n_f = log(1/delta)/(xi-E_1),
xi-E_1 = 14.43.
"""
import math, os, sys
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "data")
KS = ('22.4000', '23.0000', '23.9000')

def load(f):
    out = {}
    for line in open(f):
        p = line.split()
        out[(p[0][:7], int(p[1]))] = tuple(p[2:8])
    return out

base = load(os.path.join(D, "vectors_n120.txt"))
print("E,logdelta,pred_nf,nf_22.4,nf_23.0,nf_23.9")
for E in (40, 80, 160, 320):
    f = os.path.join(D, "vec_E%d.txt" % E)
    if not os.path.exists(f):
        print("%d,,, (dump missing -- rerun run_surrogate.gp)" % E); continue
    S = load(f); out = []
    for k in KS:
        nf = next((n for n in range(4, 121) if base.get((k, n)) != S.get((k, n))), None)
        out.append(nf)
    ld = 2 * E * math.log(10)
    print("%d,%.1f,%.1f,%s,%s,%s" % (E, -ld, ld / 14.43, *out))
