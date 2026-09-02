"""Task 4: optimise CDT's entry and margin over ALL inventories, per host.

tau depends on the inventory only through (m, u_1..u_k, sum e, max e):
    tau^flat = 2k - (2/m^2) sum_j u_j^2 ,   u_j = #{atoms with < j layers}
    tau^#    = (2/m^2) min_xi { xi*sum(e) + max(e)*I_xi^m(xi) }
so the search is over (m, a, n2) plus the exponent assignment, which is greedy
against the measured supply.

SUPPLY per level L (number of Q(y)-independent admissible atoms with exponent e):
  L=0 : one atom, e=0                                            [PROVED]
  L=1 : nu per exponent e=0,1,2,...   (measured nu=1: the arcsine tower F_{e+1})
  L=2 : s2 = [3,2,1,1,1,...] (measured N(2,e) increments; 4 at e=0 if B_5 counts)
  L=k : pure qk at e=0 plus the weight tower, and per conditional generator
        D atoms at e=0, I at e=1, I2 at e=2, e=3, ...            [VERIFIED D=4,I=3]
'UNLIMITED' replaces every cap by infinity: a rigorous upper bound on what ANY
inventory of admissible functions could achieve.
"""
import math
from fasttau import tau_sharp_f

LOSS = math.log(0.6292232680)

def greedy(count, supply):
    """cheapest total (sum e, max e) for `count` atoms drawn from supply[e] (a list;
    the last entry repeats forever).  Returns (sum_e, max_e) or None if impossible."""
    se, me, left, e = 0, 0, count, 0
    while left > 0:
        s = supply[min(e, len(supply)-1)]
        if s <= 0:
            e += 1
            if e > 400: return None
            continue
        t = min(left, s)
        se += t*e; me = max(me, e) if t else me
        left -= t; e += 1
        if e > 400: return None
    return se, me

def evaluate(k, m, a, n2, sup1, sup2, supk, ceil, BC):
    ntop = m - 1 - a - n2
    if ntop < 1: return None                    # at least the conditional generator
    r1 = greedy(a, sup1) if a else (0, 0)
    r2 = greedy(n2, sup2) if n2 else (0, 0)
    rk = greedy(ntop, supk)
    if r1 is None or r2 is None or rk is None: return None
    se = r1[0]+r2[0]+rk[0]
    me = max(r1[1], r2[1], rk[1])
    u = [1, 1+a] if k == 2 else [1, 1+a, 1+a+n2]
    tf = 2*k - (2.0/(m*m))*sum(x*x for x in u)
    ts, _ = tau_sharp_f(m, se, me)
    T = tf+ts
    ec, er = ceil-T, ceil+LOSS-T
    return (ec, er, m*er-BC, m, a, n2, ntop, tf, ts, T, se, me)

def scan(k, ceil, BC, D, I, I2=3, g=1, nu=1, unlimited=False, mmax=90, amax=20, q2=3, qk=None):
    if qk is None: qk = 2 if k == 2 else 4
    BIG = 10**6
    sup1 = [BIG] if unlimited else [nu]
    sup2 = [BIG] if unlimited else [q2, 2, 1]
    supk = [BIG] if unlimited else [qk+D*g, I*g+2, I2*g+1, I2*g+1, 1]
    best_e = best_m = None
    for m in range(2, mmax+1):
        for a in range(0, min(amax, m-1)+1):
            n2range = range(0, min(m-1-a, 41)) if k == 3 else [0]
            for n2 in n2range:
                r = evaluate(k, m, a, n2, sup1, sup2, supk, ceil, BC)
                if r is None: continue
                if best_e is None or r[0] > best_e[0]: best_e = r
                if best_m is None or r[2] > best_m[2]: best_m = r
    return best_e, best_m

def show(tag, r):
    ec, er, mar, m, a, n2, ntop, tf, ts, T, se, me = r
    print(f"   {tag:26s} m={m:3d} (1 + {a} L1 + {n2} L2 + {ntop} Lk)  sum_e={se:3d} max_e={me:2d}  "
          f"tf={tf:.4f} ts={ts:.4f} tau={T:.4f} | entryC={ec:+.4f} entryR={er:+.4f} margin={mar:+.3f}")

HOSTS = [
 ('H1 CDT level 6   (lam2=1, k=2, s=1)',    2, math.log(256), 11.845,                4, 3),
 ('H2 Catalan lvl 8 (lam2=4, k=2, s=1/4)',  2, math.log(64),  11.845+math.log(0.25), 4, 3),
 ('H3 X_1(5) Sym^2  (normalised, k=3)',     3, math.log(256), 11.845,                6, 3),
]

if __name__ == '__main__':
    print("CONTROL -- CDT's own inventory (m=14, a=2, top block 6 at e=0 + 5 at e=1)")
    ts, _ = tau_sharp_f(14, 6, 1)
    tf = 4 - (2.0/196)*(1+9)
    lp = math.log(256*0.6292232680)
    print(f"   tau^flat={tf:.6f} (191/49={191/49:.6f})  tau^#={ts:.6f} (27/80={27/80})  "
          f"tau={tf+ts:.6f} (16603/3920={16603/3920:.6f})")
    print(f"   entry={lp-tf-ts:.6f}  m <= {11.845/(lp-tf-ts):.4f} (CDT 13.9938)  "
          f"margin={14*(lp-tf-ts)-11.845:+.4f}\n")
    for nm, k, ceil, BC, D, I in HOSTS:
        print(f"=== {nm}  ceiling={ceil:.4f}  BC={BC:.4f}  D={D} I={I} ===")
        for lab, kw in (('measured supply, g=1', dict(g=1)),
                        ('measured supply, g=3', dict(g=3)),
                        ('UNLIMITED (upper bd)', dict(unlimited=True)),
                        ('nu=2 (2 per exponent)', dict(g=3, nu=2)),
                        ('nu=3 (3 per exponent)', dict(g=3, nu=3))):
            be, bm = scan(k, ceil, BC, D, I, **kw)
            show(lab+' [max entry]', be)
            show(lab+' [max margin]', bm)
        print()
