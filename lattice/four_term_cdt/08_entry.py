#!/usr/bin/env python3
"""08_entry.py -- the CDT entry condition and margin for the mixed-exponent
four-term Catalan / L(2,chi_-3) hosts, in the conventions of CDT_FINDER.md 1-4
and ADELIC_HOLONOMY.md 2.6 (Theorem A).

   tau^flat(m,u,b) = sum_j b_j - (1/m^2) sum_j u_j^2 b_j
   tau^sharp       = 0.3375        (CDT's value; depends on e and m only)
   entry           = log|varphi'(0)| + (1 - 1/m) sum_p log R_p - tau
   margin          = m * entry - BC(varphi)

Two accountings are reported for every host:
  ceiling  : log|varphi'(0)| = log(size), the exact orbifold cusp size computed in
             06_sizes.py -- a genuine supremum, not attained.
  contour  : log(0.62922 * size), CDT's own realised contour loss transported;
             BC = 6.763 + log(0.62922*size), also transported.  Both flagged
             [estimate] -- the shape of our 4-point orbifold is not CDT's.

The four-term hosts are UNSYMMETRISED (no involution of P^1 preserves their
singular set -- checked in 08_entry.py:involutions), so b = (1,1), sigma_m = 2.
Level 8 (Zagier E) is symmetrised, b = (2,2), sigma_m = 4: that is where its
tau = 4.2355 and its ceiling 256/4 = 64 come from.
"""
import os, sys, json, math
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))

TAU_SHARP = 27.0/80          # CDT, m = 14, e = (0,0,1;0^6;1,1,1,1,1)
CDT_LOSS = 0.62922           # |psi'(0)| of CDT's 4-slit+lune Riemann map
CDT_SHAPE = 11.845 - 5.081908  # BC - log|varphi'(0)| on CDT's own contour


def tau_flat(m, u, b):
    return sum(b) - sum(uj*uj*bj for uj, bj in zip(u, b))/float(m*m)


def involutions(pts):
    """all Moebius involutions of P^1 fixing 0 and permuting {t_1,t_2,t_3,infinity}.
    w(x) = s x/(x - s) fixes 0, swaps s <-> infinity; x -> -x fixes 0 and infinity."""
    ts = [complex(p['tre'], p['tim']) for p in pts]
    found = []
    for s in ts:
        img = []
        for t in ts:
            if abs(t - s) < 1e-12:
                img.append('inf')
            else:
                img.append(s*t/(t - s))
        rest = [t for t in ts if abs(t - s) > 1e-12]
        got = [v for v in img if v != 'inf']
        ok = all(any(abs(v - t) < 1e-9 for t in rest) for v in got)
        if ok:
            found.append(("w(x)=%s x/(x-%s)" % (s, s)))
    if all(any(abs(-t - u) < 1e-9 for u in ts) for t in ts):
        found.append("x -> -x")
    return found


def main():
    geo = {g['label']: g for g in json.load(open(os.path.join(HERE, 'out', 'geometry.json')))}
    sizes = {s['label']: s for s in json.load(open(os.path.join(HERE, 'out', 'sizes.json')))}
    padic = {p['label']: p for p in json.load(open(os.path.join(HERE, 'out', 'padic.json')))}
    XI = {"R1": "G/4", "R2": "G/2 - 3z(2)/16", "R3": "3z(2)/8 - G/2",
          "R4": "(2z(2)+15L)/32", "R5": "(15L-6z(2))/16", "R6": "(2z(2)-3L)/8",
          "R7": "(no arch. limit)"}
    m = 14
    cfgA = dict(name="CDT-transported inventory, unsymmetrised",
                u=(1, 3), b=(1, 1))
    cfgB = dict(name="best conceivable inventory (u_j = m/2), unsymmetrised",
                u=(m/2., m/2.), b=(1, 1))
    tauA = tau_flat(m, cfgA['u'], cfgA['b']) + TAU_SHARP
    tauB = tau_flat(m, cfgB['u'], cfgB['b']) + TAU_SHARP
    tauSym = tau_flat(m, (1, 3), (2, 2)) + TAU_SHARP
    print("tau (A: CDT inventory, b=(1,1)) = %.6f" % tauA)
    print("tau (B: best inventory, b=(1,1)) = %.6f" % tauB)
    print("tau (level-8 symmetrised, b=(2,2)) = %.6f   [CDT's 16603/3920]" % tauSym)
    print()
    rows = []
    for lab in ["R1", "R2", "R3", "R4", "R5", "R6", "R7"]:
        if lab not in sizes:
            continue
        g, s, pa = geo[lab], sizes[lab], padic[lab]
        # adelic: R_p = min over the inventory {A, B, H}
        sl2 = min(pa['p2']['slope_a'], pa['p2']['slope_b'], pa['p2']['slope_h'])
        adel = max(0.0, sl2) * math.log(2)
        adel = adel if adel > 0.02 else 0.0
        gain = (1 - 1.0/m)*adel
        lsz = math.log(s['size'])
        for tag, tau in (("A", tauA), ("B", tauB)):
            eC = lsz + gain - tau
            eR = lsz + math.log(CDT_LOSS) + gain - tau
            BC = CDT_SHAPE + lsz + math.log(CDT_LOSS)
            rows.append(dict(label=lab, cfg=tag, xi=XI[lab], size=s['size'],
                             logsize=lsz, ceil16=s['ceil16'], loss=s['loss'],
                             adelic=gain, tau=tau, entry_ceiling=eC,
                             entry_contour=eR, margin=m*eR - BC))
        inv = involutions(g['pts'])
        print("%-3s xi = %-18s  size %.6f (log %+0.5f)  ceil16 %.4f  4th-point loss %.4f  "
              "adelic +%.4f  involutions: %s"
              % (lab, XI[lab], s['size'], lsz, s['ceil16'], s['loss'], gain, inv or "NONE"))
    print()
    print("%-3s %-4s %-9s %-9s %-9s %-9s" % ("row", "cfg", "log|f'(0)|", "tau", "entry(ceil)", "margin"))
    for r in rows:
        print("%-3s %-4s %+9.5f %9.5f %+11.5f %+10.3f"
              % (r['label'], r['cfg'], r['logsize'] + r['adelic'], r['tau'],
                 r['entry_ceiling'], r['margin']))
    # benchmarks
    print()
    print("benchmarks (same conventions):")
    for nm, sz, tau in (("level 8 (Zagier E), symmetrised, k=2", 64.0, tauSym),
                        ("level 8 (Zagier E), unsymmetrised", 4.0, tauA),
                        ("level 16 Catalan, symmetrised", 16.0*0.25*4, tauSym),
                        ("CDT X_0(6), symmetrised", 256.0, tauSym)):
        print("   %-40s size %8.3f  entry(ceiling) = %+0.5f" % (nm, sz, math.log(sz) - tau))
    json.dump(rows, open(os.path.join(HERE, 'out', 'entry.json'), 'w'), indent=1)
    print("\nwrote out/entry.json")


if __name__ == '__main__':
    main()
