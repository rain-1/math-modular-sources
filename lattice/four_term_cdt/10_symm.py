#!/usr/bin/env python3
"""10_symm.py -- is a CDT-style symmetrisation available on these hosts?

CDT's +log 16 comes from descending P^1 - {0,s,infinity} along the involution
w(x) = s x/(x-s) (fixes 0, swaps s <-> infinity), quotient coordinate
y = x + w(x) = x^2/(x-s): the cusp at 0 acquires half its width, so the cusp size
goes 16 s -> 256 s = 64 * (4s), 4s being the position of the order-2 branch point
in y.  (Verified in 05_uniformise: the orbifold {0 cusp, infinity cusp, c order-2}
has cusp size exactly 64 c, and 64 * 4s = 256 s reproduces CDT's ceiling.)
The price is that y ~ x^2, so every denominator layer [1..n] becomes [1..2n]:
sigma_m doubles and tau goes 2.2865 -> 4.2355.

Here we enumerate the involutions of P^1 that fix 0 and permute the singular set
type-preservingly, BEFORE and AFTER the hypothesis, and price the survivors.
"""
import os, sys, json, math
HERE = os.path.dirname(os.path.abspath(__file__))


def invols(pts):
    """pts = [(z, kind)]; kind in {'cusp','orb2'}; 0 and infinity are cusps.
    Involutions fixing 0: w_s(x) = s x/(x-s) (swaps s <-> infinity) for s a CUSP,
    and x -> -x."""
    out = []
    for s, ks in pts:
        if ks != 'cusp':
            continue                      # w_s sends s to infinity, a cusp
        img = []
        ok = True
        for t, kt in pts:
            if abs(t - s) < 1e-12:
                continue
            v = s*t/(t - s)
            hit = [(u, ku) for u, ku in pts if abs(u - v) < 1e-9]
            if not hit or hit[0][1] != kt:
                ok = False; break
            img.append((t, v))
        if ok:
            out.append(("w(x) = %g x/(x - %g)" % (s.real, s.real), s, img))
    return out


def main():
    geo = json.load(open(os.path.join(HERE, 'out', 'geometry.json')))
    for g in geo:
        pts = [(complex(p['tre'], p['tim']), p['tag']) for p in g['pts']]
        post = pts[1:]
        pre_i = invols(pts)
        post_i = invols(post)
        print("%s" % g['label'])
        print("   pre-hypothesis  (5 points): %s"
              % ("; ".join(x[0] for x in pre_i) or "none"))
        print("   post-hypothesis (4 points): %s"
              % ("; ".join(x[0] for x in post_i) or "none"))
        for nm, s, img in pre_i:
            y = lambda x: (x*x/(x - s)) if abs(x - s) > 1e-12 else complex('inf')
            fold = pts[0][0]
            if abs(fold - s) < 1e-12:
                print("      %s : the FOLD IS the involution centre s = %g, which the"
                      " quotient identifies with the cusp at infinity -- the hypothesis"
                      " cannot delete it.  Symmetrisation unavailable." % (nm, s.real))
                continue
            print("      %s : quotient y = x^2/(x-%g); branch point y=%g (order 2);"
                  " fold %g and its partner -> y = %g"
                  % (nm, s.real, (4*s).real, fold.real, y(fold).real))
            print("      -> Sym^+ H is still branched at y = %g (H is regular at the fold"
                  " but NOT at its involution partner), so the symmetrised"
                  " post-hypothesis geometry keeps FOUR points:" % y(fold).real)
            print("         {0 cusp, oo cusp, %g order-2 (branch), %g order-2 (fold pair)}"
                  % ((4*s).real, y(fold).real))
    print()
    print("Only R1 and R4 admit an involution at all, and in both cases the fold's")
    print("partner survives the quotient.  For R1 the symmetrised post-hypothesis")
    print("geometry is {0, oo cusps; +-1/2 order-2}, whose cusp size is exactly 4")
    print("(the z -> z^2 cover of {0, oo cusps; 1 order-2}, size 64/8 = 8, at scale 1/2)")
    print("-- verified numerically in 05_uniformise.  Entry = log 4 - 4.235459 = %+0.5f,"
          % (math.log(4) - 4.235459))
    print("worse than the unsymmetrised -1.914.  For R4 the fold IS the involution")
    print("centre s = 1/8, which the quotient identifies with the cusp at infinity, so")
    print("the hypothesis cannot delete it: symmetrisation is unavailable there.")


if __name__ == '__main__':
    main()
