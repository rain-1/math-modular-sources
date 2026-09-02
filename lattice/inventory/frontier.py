"""Margin/entry frontier as a function of m, measured supply, g=1."""
import math
from final import best, ev, HOSTS
from optimise import greedy
BIG = 10**6
def bestm(k, ceil, BC, D, I, m, g=1, nu=1, q2=3, qk=None, I2=3):
    if qk is None: qk = 2 if k == 2 else 4
    sup1, sup2 = [nu], [q2, 2, 1]
    supk = [qk+D*g, I*g+2, I2*g+1, I2*g+1, 1]
    be = bm = None
    for a in range(0, m):
        n2r = range(0, m-1-a) if k == 3 else [0]
        for n2 in n2r:
            r = ev(k, m, a, n2, sup1, sup2, supk, ceil, BC)
            if not r: continue
            if be is None or r['ec'] > be['ec']: be = r
            if bm is None or r['marR'] > bm['marR']: bm = r
    return be, bm
print("| m | H1 entry@ceil | H1 margin | H2 entry@ceil | H2 entry@contour | H2 margin | "
      "H3 entry@ceil | H3 entry@contour | H3 margin |")
print("|---|---|---|---|---|---|---|---|---|")
for m in list(range(2, 21))+[24, 28, 32, 40]:
    cells = []
    for h, k, ceil, BC, D, I in HOSTS:
        be, bm = bestm(k, ceil, BC, D, I, m)
        if be is None: cells += ['-', '-', '-']; continue
        if h.startswith('H1'): cells += [f"{be['ec']:+.3f}", f"{bm['marR']:+.2f}"]
        else: cells += [f"{be['ec']:+.3f}", f"{be['er']:+.3f}", f"{bm['marR']:+.2f}"]
    print("| " + str(m) + " | " + " | ".join(cells) + " |")
