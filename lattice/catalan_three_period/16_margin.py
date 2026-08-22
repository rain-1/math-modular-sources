"""Adelic margins on the level-16 Catalan host with the corrected inventory.

New inputs (this task, lattice/catalan_three_period/):
  * fold cusp of the level-16 host = the cusp 1/2 of Gamma_0(16) (x=-1/4),
    width 4, gamma_1 = [-7,4;-16,9].
  * inner chi_{-4} classes  Phi=(c1+c2 V2+c4 V4)E:  fold-regular <=> 8c2+c4=0
    (2-dim), xi = -(1/2)P(2) G ; target zero: Phi_0^in=(1+4V2-32V4)E.
  * outer classes Phi=(c1+c2V2+c4V4)T:  the period polynomial always satisfies
    2 alp + bet = 0, but the companion carries the constant term a_0=-P(0)/4;
    fold-regular <=> P(0)=0 (2-dim), and then xi = (3/2) zeta(2) (P(1)-P(2)).
    Target zero: P(0)=0 and P(1)=P(2)  =>  Phi_0^out = (1+3V2-4V4)T.  NEW.
  * measured 2-adic slopes in y: pure 2, inner doubly-small 1, outer
    doubly-small 0, conditional 0.
  * Q(y)-rank: {g1,g2,g3} and their CDT 7-member orbits are independent (21/21
    up to degree 2); 6 of 7 candidate pure functions independent.
"""
import math, sys
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/adelic_holonomy')
from adelic_bound import adelic, fmt

ceil_ = math.log(256*0.25)                 # 4.158883  uniformisation ceiling
real_ = ceil_ + math.log(0.6292232680)     # 3.695614  CDT contour loss, TRANSPORTED (flagged)
rig_  = math.log(4.376879)                 # 1.476336  thrice-punctured, rigorous but useless
BC    = 11.845 + math.log(0.25)            # 10.458706

def ev(tag, m, cols, e, slopes, lp):
    r = adelic(m, cols, e, ({2: slopes} if slopes else {}), lp, BC)
    fmt(tag, r, m)
    return r

def eprof(npure, nsmallorbits, extra=0):
    """CDT-proportional integration profile: 3 of every 7 members carry n^1."""
    e = [0]*npure
    for j in range(min(3,npure)): e[-(j+1)] = 1
    for _ in range(nsmallorbits): e += [0,0,0,0,1,1,1]
    e += [0]*extra
    return e

cols = [(1,2),(3,2)]          # u1=1, u2=3, b1=b2=2  (CDT's, transported)
print(f"BC={BC:.6f}  ceiling={ceil_:.6f}  realised(transported)={real_:.6f}  rigorous={rig_:.6f}\n")

print("--- reference rows (previous documents) ---")
e14=[0,0,1,0,0,0,0,0,0,1,1,1,1,1]
ev("level 8/16, m=14, pure slope 2, ceiling ", 14, cols, e14, [2]*7+[0]*7, ceil_)
ev("level 16, m=15 (+1 doubly-small), ceiling",15, cols, e14+[0], [2]*7+[1]+[0]*7, ceil_)
print()

print("--- corrected, conservative: one member per small generator (m=16) ---")
ev("m=16: 7 pure + 7 cond + 2 dbl, ceiling ", 16, cols, e14+[0,0], [2]*7+[1,0]+[0]*7, ceil_)
ev("m=16, realised(transported)            ", 16, cols, e14+[0,0], [2]*7+[1,0]+[0]*7, real_)
print()

print("--- corrected, CDT orbit multiplier (7 members per small generator) ---")
print("    hypotheses (a) G in Q  and  (b) a+b zeta(2)+c G=0 give the SAME m.")
sl28 = [2]*7 + [1]*7 + [0]*14
e28  = eprof(7,3)
ev("m=28: 7 pure + 3x7 small, ceiling      ", 28, cols, e28, sl28, ceil_)
ev("m=28, realised(transported)            ", 28, cols, e28, sl28, real_)
ev("m=28, rigorous thrice-punctured        ", 28, cols, e28, sl28, rig_)
print()

print("--- the UNCONDITIONAL sub-inventory (must NOT give a positive margin) ---")
ev("m=21: 7 pure + 2x7 doubly-small, ceil  ", 21, cols, eprof(7,2), [2]*7+[1]*7+[0]*7, ceil_)
ev("m=20: 6 pure + 2x7 doubly-small, ceil  ", 20, cols, eprof(6,2), [2]*6+[1]*7+[0]*7, ceil_)
ev("m=14: 2x7 doubly-small only, ceiling   ", 14, [(0,2),(0,2)], eprof(0,2), [1]*7+[0]*7, ceil_)
ev("m=7 : inner doubly-small orbit, ceiling", 7,  [(0,2),(0,2)], eprof(0,1), [1]*7, ceil_)
print()

print("--- one-period vs three-period, conditional orbit only added ---")
for m_, tag, sl, e_ in [
    (21, "one-period (a): 7 pure+7 dbl_in+7 cond ", [2]*7+[1]*7+[0]*7, eprof(7,2)),
    (28, "one-period (a): + outer dbl orbit      ", sl28, e28),
    (28, "three-period (b): identical inventory  ", sl28, e28)]:
    ev(tag, m_, cols, e_, sl, ceil_)
print()

print("--- sensitivity: what m would be needed for margin 0 at the ceiling ---")
for m_ in (28, 34, 42, 49, 56, 70):
    k = (m_-7)//7
    ev(f"m={m_:3d} ({k} small orbits + 7 pure)      ", m_, cols, eprof(7,k, m_-7-7*k),
       [2]*7+[1]*7+[0]*(m_-14), ceil_)

print()
print("--- optimisation over admissible sub-inventories (ceiling) ---")
print("    available: 7 pure (slope 2, u=(1,3)), 7 inner doubly-small (slope 1),")
print("               7 outer doubly-small (slope 0), 7 conditional (slope 0).")
best = None; bestu = None
rows = []
for n2 in range(0, 8):
    for n1 in range(0, 8):
        for nd0 in range(0, 8):
            for nc in range(0, 8):
                m_ = n2+n1+nd0+nc
                if m_ < 2: continue
                u1 = 1 if n2 >= 1 else 0
                u2 = min(3, n2)
                c_ = [(u1,2),(u2,2)]
                ne = min(3, n2) + (3 if n1 == 7 else max(0,n1-4)) \
                     + (3 if nd0 == 7 else max(0,nd0-4)) + (3 if nc == 7 else max(0,nc-4))
                e_ = [1]*ne + [0]*(m_-ne)
                sl = [2]*n2 + [1]*n1 + [0]*(nd0+nc)
                r = adelic(m_, c_, e_, {2: sl}, ceil_, BC)
                rows.append((r['margin'], m_, n2, n1, nd0, nc, r))
                if nc >= 1 and (best is None or r['margin'] > best[0]):
                    best = (r['margin'], m_, n2, n1, nd0, nc, r)
                if best is not None: pass
                if bestu is None or (nc == 0 and r['margin'] > bestu[0]):
                    bestu = (r['margin'], m_, n2, n1, nd0, nc, r)
rows.sort(reverse=True)
print("  top 8 sub-inventories overall (n2,n1,nd0,nc = pure, inner-dbl, outer-dbl, conditional):")
for mg, m_, n2, n1, nd0, nc, r in rows[:8]:
    print(f"    m={m_:3d}  (n2,n1,nd0,nc)=({n2},{n1},{nd0},{nc})  tau={r['tau']:.4f} "
          f"gamma={r['gamma']:+.4f} entry={r['entry']:+.4f} margin={mg:+.3f}")
mg,m_,n2,n1,nd0,nc,r = best
print(f"  best WITH a conditional function: m={m_} ({n2},{n1},{nd0},{nc}) margin={mg:+.3f}")
mg,m_,n2,n1,nd0,nc,r = bestu
print(f"  best unconditional               : m={m_} ({n2},{n1},{nd0},{nc}) margin={mg:+.3f}"
      "   (must be <=0, else the inputs are inconsistent)")

print()
print("--- headline table: tau^flat, tau^sharp, gamma_2, entry, margin ---")
def row(tag, m_, n2, n1, nd0, nc, lp):
    u1 = 1 if n2 >= 1 else 0; u2 = min(3, n2)
    c_ = [(u1,2),(u2,2)]
    ne = min(3,n2) + (3 if n1==7 else max(0,n1-4)) + (3 if nd0==7 else max(0,nd0-4)) \
         + (3 if nc==7 else max(0,nc-4))
    e_ = [1]*ne + [0]*(m_-ne)
    sl = [2]*n2 + [1]*n1 + [0]*(nd0+nc)
    r = adelic(m_, c_, e_, {2: sl}, lp, BC)
    print(f"  {tag:<42s} m={m_:3d}  tf={r['tau_flat']:.4f} ts={r['tau_sharp']:.4f} "
          f"tau={r['tau']:.4f} g2={r['gamma']:+.4f} entry={r['entry']:+.4f} margin={r['margin']:+8.3f}")
    return r
for lp, lab in ((ceil_, "ceiling"), (real_, "realised(transported)"), (rig_, "rigorous 3-punctured")):
    print(f" [{lab}]")
    row("previous: 7 pure + 7 conditional",        14, 7,0,0,7, lp)
    row("previous: + 1 doubly-small",              15, 7,1,0,7, lp)
    row("(a)/(b) full: 7 pure + 3 small orbits",   28, 7,7,7,7, lp)
    row("(a)/(b) best with a conditional",         15, 7,7,0,1, lp)
    row("unconditional only, best",                14, 7,7,0,0, lp)
