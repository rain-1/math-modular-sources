"""Task 4/5 final tables."""
import math
from fasttau import tau_sharp_f
from optimise import greedy, LOSS

BIG = 10**6
def ev(k, m, a, n2, sup1, sup2, supk, ceil, BC):
    ntop = m-1-a-n2
    if ntop < 1: return None
    r1 = greedy(a, sup1) if a else (0, 0)
    r2 = greedy(n2, sup2) if n2 else (0, 0)
    rk = greedy(ntop, supk)
    if None in (r1, r2, rk): return None
    se = r1[0]+r2[0]+rk[0]; me = max(r1[1], r2[1], rk[1])
    u = [1, 1+a] if k == 2 else [1, 1+a, 1+a+n2]
    tf = 2*k - (2.0/(m*m))*sum(x*x for x in u)
    ts, _ = tau_sharp_f(m, se, me); T = tf+ts
    ec, er = ceil-T, ceil+LOSS-T
    return dict(ec=ec, er=er, marR=m*er-BC, marC=m*ec-BC, m=m, a=a, n2=n2, ntop=ntop,
                tf=tf, ts=ts, T=T, se=se, me=me)

def best(k, ceil, BC, D, I, I2=3, g=1, nu=1, q2=3, qk=None, top_unlimited=False,
         nu_inf=False, mmax=110, amax=40, key='marR'):
    if qk is None: qk = 2 if k == 2 else 4
    sup1 = [BIG] if nu_inf else [nu]
    sup2 = [BIG] if top_unlimited else [q2, 2, 1]
    supk = [BIG] if top_unlimited else [qk+D*g, I*g+2, I2*g+1, I2*g+1, 1]
    b = None
    for m in range(2, mmax+1):
        for a in range(0, min(amax, m-1)+1):
            n2r = (sorted(set(list(range(0, min(m-1-a, 17)))+[max(0, m-2-a)])) if k == 3 else [0])
            for n2 in n2r:
                r = ev(k, m, a, n2, sup1, sup2, supk, ceil, BC)
                if r and (b is None or r[key] > b[key]): b = r
    return b

HOSTS = [('H1 CDT lvl6',   2, math.log(256), 11.845,                4, 3),
         ('H2 Catalan lvl8', 2, math.log(64),  11.845+math.log(0.25), 4, 3),
         ('H3 X1(5) Sym^2',  3, math.log(256), 11.845,                6, 3)]

def line(h, lab, r):
    if r is None: return f"| {h} | {lab} | - |"
    return (f"| {h} | {lab} | {r['m']} | {r['a']} | {r['n2']} | {r['ntop']} | {r['se']} | {r['me']} | "
            f"{r['tf']:.4f} | {r['ts']:.4f} | {r['T']:.4f} | {r['ec']:+.4f} | {r['er']:+.4f} | "
            f"{r['marR']:+.3f} | {r['marC']:+.3f} |")

if __name__ == '__main__':
    print("| host | scenario / objective | m | a | n2 | top | Se | maxe | tau^flat | tau^# | tau | "
          "entry@ceil | entry@contour | margin@contour | margin@ceil |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    SC = [('B measured, g=1', dict(g=1)), ('B measured, g=2', dict(g=2)), ('B measured, g=3', dict(g=3)),
          ('C nu=1, unlimited top', dict(g=1, top_unlimited=True)),
          ('D nu=2, unlimited top', dict(g=1, nu=2, top_unlimited=True)),
          ('D nu=3, unlimited top', dict(g=1, nu=3, top_unlimited=True)),
          ('E nu=inf (absurd)',     dict(g=1, nu_inf=True, top_unlimited=True))]
    for h, k, ceil, BC, D, I in HOSTS:
        for lab, kw in SC:
            for key, tag in (('ec', 'max entry'), ('marR', 'max margin'), ('marC', 'max margin@ceil')):
                print(line(h, f"{lab} : {tag}", best(k, ceil, BC, D, I, key=key, **kw)))
