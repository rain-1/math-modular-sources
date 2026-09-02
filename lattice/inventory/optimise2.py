"""Task 4/5: the definitive scans.  Scenarios, in increasing order of optimism:

 A  CDT-transported            : CDT's own (m=14) inventory moved to the host
 B  measured supply            : nu=1 at L=1; q2=2 (+F_2) at (2,0), 2 at (2,1), 1 per e after;
                                 conditional D derivatives (e=0) + I integrals (e=1) per generator
 C  measured L=1, UNLIMITED top: nu=1 single-layer supply (the PROVED/VERIFIED constraint)
                                 but an unbounded supply of e=0 atoms at the top level
 D  nu=2 / nu=3                : what if the single-layer supply were 2 or 3 per exponent
 E  nu=inf                     : the absurd upper bound (u_2 = m-1)
"""
import math
from fasttau import tau_sharp_f
from optimise import greedy, evaluate, LOSS, HOSTS

BIG = 10**6
def scan(k, ceil, BC, D, I, I2=3, g=1, nu=1, top_unlimited=False, nu_inf=False,
         mmax=110, amax=40, q2=3, qk=None):
    if qk is None: qk = 2 if k == 2 else 4
    sup1 = [BIG] if nu_inf else [nu]
    sup2 = [BIG] if top_unlimited else [q2, 2, 1]
    supk = [BIG] if top_unlimited else [qk+D*g, I*g+2, I2*g+1, I2*g+1, 1]
    best_e = best_m = None
    for m in range(2, mmax+1):
        for a in range(0, min(amax, m-1)+1):
            n2r = (sorted(set(list(range(0, min(m-1-a, 17))) + [max(0, m-2-a)]))
                   if k == 3 else [0])
            for n2 in n2r:
                r = evaluate(k, m, a, n2, sup1, sup2, supk, ceil, BC)
                if r is None: continue
                if best_e is None or r[0] > best_e[0]: best_e = r
                if best_m is None or r[2] > best_m[2]: best_m = r
    return best_e, best_m

def row(host, lab, r):
    ec, er, mar, m, a, n2, ntop, tf, ts, T, se, me = r
    return (f"| {host} | {lab} | {m} | {a} | {n2} | {ntop} | {se} | {me} | {tf:.4f} | {ts:.4f} | "
            f"{T:.4f} | {ec:+.4f} | {er:+.4f} | {mar:+.3f} |")

SCEN = [('B measured supply, g=1', dict(g=1)),
        ('B measured supply, g=3', dict(g=3)),
        ('C nu=1, unlimited top',  dict(g=1, top_unlimited=True)),
        ('D nu=2, unlimited top',  dict(g=1, nu=2, top_unlimited=True)),
        ('D nu=3, unlimited top',  dict(g=1, nu=3, top_unlimited=True)),
        ('E nu=inf (absurd)',      dict(g=1, nu_inf=True, top_unlimited=True))]

print("| host | scenario | m | a(L1) | n2(L2) | top | sum e | max e | tau^flat | tau^# | tau | "
      "entry@ceil | entry@CDTcontour | margin |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for nm, k, ceil, BC, D, I in HOSTS:
    short = nm.split()[0]
    # scenario A: CDT transported
    ts, _ = tau_sharp_f(14, 6, 1); tf = 4 - (2.0/196)*(1+9)
    if k == 3:
        tf = 6 - (2.0/196)*(1+9+25); ts, _ = tau_sharp_f(14, 6, 1)
    T = tf+ts
    print(row(short, 'A CDT-transported (m=14)',
              (ceil-T, ceil+LOSS-T, 14*(ceil+LOSS-T)-BC, 14, 2, 0, 11, tf, ts, T, 6, 1)))
    for lab, kw in SCEN:
        be, bm = scan(k, ceil, BC, D, I, **kw)
        print(row(short, lab+' [maxE]', be))
        print(row(short, lab+' [maxM]', bm))
