"""Candidate table: for each projector F in the span, its singular set in z,
the width-law ceiling log(16|t2|), the denominator growth of F and of its fold
companion B_F = int_0^z F dt/(2-t), the p-adic slopes, and the entry number.

Fold: every F in the span has F(2) in Q*G (H(2)=2G etc.), and z=2 is a puncture
of the local system (the pure companion 2 arcsin sqrt(z/2) branches there), so
B_F - F(2)*A is fold-regular under G in Q, with A = -log(1-z/2).
"""
from fractions import Fraction as Fr
import math, json, sys, os
from emn_core import psi_basis, ratseries, solve_from_psi
from emn_cand import psi_series, build, vp

N = int(sys.argv[1]) if len(sys.argv) > 1 else 400
NS = [N//4, N//2, N]

def poles_of(comb):
    """returns sorted list of z-values of the poles (theta in (0,1))."""
    from emn_span import basis_residues
    A = max(a for (a,k) in comb)
    res = basis_residues(A)
    thetas = sorted({t for (a,k) in comb for t in res[(a,'0' if k=='0' else 'p3')]})
    out = []
    for t in thetas:
        r = sum(comb[(a,k)]*res[(a,'0' if k=='0' else 'p3')].get(t, Fr(0)) for (a,k) in comb)
        if r != 0:
            out.append((float(t), 1-math.cos(math.pi*float(t)), r))
    return out

def growth(c):
    return {n: (math.log(c[n].denominator)/n if c[n] != 0 else 0.0) for n in NS}

def slope(c, p):
    v = [vp(c[n], p) for n in NS]
    return v[-1]/NS[-1] if v[-1] is not None else 0.0

def analyse(name, comb):
    F = build(comb, N)
    pl = poles_of(comb)
    if not pl:
        return dict(name=name, const=True)
    t2 = min(p[1] for p in pl)
    # fold companion in the z coordinate: B = int_0^z F dt/(2-t)
    part = [Fr(0)]*(N+1); run = Fr(0)
    for n in range(N+1):
        part[n] = run/Fr(2)**1 if False else run
        run = run/Fr(2) + F[n]      # running sum of F_k / 2^{n-k}
    # redo cleanly: [z^{n}] (F/(2-z)) = sum_{k<=n} F_k / 2^{n-k+1}
    conv = [Fr(0)]*(N+1); acc = Fr(0)
    for n in range(N+1):
        acc = acc/Fr(2) + F[n]/Fr(2)
        conv[n] = acc
    B = [Fr(0)]+[conv[n-1]/Fr(n) for n in range(1, N+1)]
    gF, gB = growth(F), growth(B)
    ceil_ = math.log(16*t2)
    tau_crude = gB[N]
    return dict(name=name, t2=t2, npoles=len(pl), poles=[round(p[1],5) for p in pl][:8],
                sigF=gF[N], sigB=gB[N], v2F=slope(F,2), v3F=slope(F,3),
                v2B=slope(B,2), ceiling=ceil_, entry=ceil_-tau_crude)

CANDS = {
  'H = h(S)                     ': {(1,'0'): Fr(1)},
  'h(2S)                        ': {(2,'0'): Fr(1)},
  'h(3S)                        ': {(3,'0'): Fr(1)},
  'hpair(S, pi/3)               ': {(1,'p3'): Fr(1)},
  'K = h(3S)+3h(S)   [note s.13]': {(3,'0'): Fr(1), (1,'0'): Fr(3)},
  'h(2S)+2h(S) [quadratic cond] ': {(2,'0'): Fr(1), (1,'0'): Fr(2)},
  'h(2S)-2h(S)                  ': {(2,'0'): Fr(1), (1,'0'): Fr(-2)},
  'h(4S)+..  (quartic cond)     ': {(4,'0'): Fr(1), (2,'0'): Fr(2), (1,'0'): Fr(4)},
  'h(3S)-3 hpair(S,pi/3)  [=-3H]': {(3,'0'): Fr(1), (1,'p3'): Fr(-3)},
  'full cubic trace             ': {(3,'0'): Fr(1), (1,'0'): Fr(3), (1,'p3'): Fr(-3)},
  'h(5S)                        ': {(5,'0'): Fr(1)},
  'h(6S)                        ': {(6,'0'): Fr(1)},
  'hpair(2S,pi/3)               ': {(2,'p3'): Fr(1)},
}
if __name__ == '__main__':
    rows = []
    for k, c in CANDS.items():
        r = analyse(k, c)
        rows.append(r)
        if r.get('const'):
            print(f"{k}  CONSTANT (Psi=0)"); continue
        print(f"{k} t2={r['t2']:.5f} #poles={r['npoles']:2d} "
              f"sig(F)={r['sigF']:.3f} sig(B)={r['sigB']:.3f} "
              f"v2(F)={r['v2F']:+.3f} v2(B)={r['v2B']:+.3f} "
              f"ceil={r['ceiling']:+.4f} ENTRY={r['entry']:+.4f}")
    json.dump([r for r in rows], open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'table.json'),'w'), indent=1)
