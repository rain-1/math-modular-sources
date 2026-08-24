#!/usr/bin/env python3
"""01_geometry.py -- exact geometry of the six mixed-exponent four-term Catalan /
L(2,chi_-3) hosts of FOUR_TERM_DEEP.md 6.4, in the CDT variables.

For each row (class (rho_p,rho_r;M,j1,j2) = (-1/2,0;1,0,0), R(n)=C n^2):
  * characteristic cubic  lam^3 - a lam^2 + d lam - g,  g = C M^2
  * the three finite singular points t_i = 1/lam_i  (plus 0 and infinity)
  * local exponent differences: 0 at the rational root (log / cusp, I_n),
    1/2 at the two roots of the quadratic factor (order-2 orbifold, Kodaira III)
  * the fold = the singular point nearest 0 (the one governing lim b_n/a_n)
  * the post-hypothesis singular set = the rest
  * the rigorous monotonicity ceiling  16 * min |t| over the post set
Everything here is exact (sympy / Fraction), no floating-point decisions.
"""
import os, sys, json
from fractions import Fraction
import sympy as sp

HERE = os.path.dirname(os.path.abspath(__file__))

# (label, r, a, c, d, f, C, xi_string, family)
ROWS = [
    ("R1", 8,  16,  8,  48,  0, -128, "G/4",                       "sqrt2"),
    ("R2", 2,  14,  8,  28,  4,    8, "G/2 - 3 zeta(2)/16",        "sqrt2"),
    ("R3", -2,  6,  4, -32, -8,   32, "3 zeta(2)/8 - G/2",         "sqrt2"),
    ("R4", 8,  16,  8,  68,  8,   32, "(2 zeta(2) + 15 L)/32",     "sqrt3"),
    ("R5", 1,  17, 10,  32,  8,   16, "(15 L - 6 zeta(2))/16",     "sqrt3"),
    ("R6", -1, 13,  8, -13, -1,   -1, "(2 zeta(2) - 3 L)/8",       "sqrt3"),
    # the seventh primitive row of the class: complex dominant pair, no arch. limit
    ("R7", 4,   8,  4,  32,  8,   64, "no archimedean limit",      "sqrt3"),
]

RHO_P = Fraction(-1, 2)
RHO_R = Fraction(0, 1)


def row_coeffs(r, a, c, d, f, C, M=1, j1=0, j2=0):
    """Theorem D3: rebuild P,Q,R of the four-term recurrence and the operator data."""
    g = C * M * M
    assert g % r == 0
    pi = g // r          # product of the two quadratic roots
    sig = a - r          # sum of the two quadratic roots
    b = (1 - RHO_R) * r + (1 - RHO_P) * sig
    e = -RHO_P * (2 * pi + r * sig) - RHO_R * r * sig
    h = -(1 + 2 * RHO_P + RHO_R) * Fraction(g)
    jj = Fraction(C * j1 * j2)
    assert d == sig * r + pi, (d, sig * r + pi)
    for x in (b, e, h, jj):
        assert x.denominator == 1
    return dict(a=a, b=int(b), c=c, d=d, e=int(e), f=f, g=g, h=int(h), j=int(jj),
                r=r, sigma=sig, pi=pi, C=C)


def analyse(label, r, a, c, d, f, C, xi, fam):
    co = row_coeffs(r, a, c, d, f, C)
    lam = sp.symbols('lam')
    chi = lam**3 - a*lam**2 + d*lam - co['g']
    quad = sp.Poly(lam**2 - co['sigma']*lam + co['pi'], lam)
    disc = sp.discriminant(chi, lam)
    roots = sp.Poly(chi, lam).all_roots()
    # classify each root: rational r -> exponent difference 0 (I_n, cusp);
    # quadratic pair -> exponent difference 1/2 (Kodaira III, order-2 orbifold)
    pts = []
    for rt in roots:
        rt = sp.nsimplify(sp.radsimp(rt))
        t = sp.simplify(1/rt)
        if sp.simplify(rt - r) == 0 and not any(p['tag'] == 'cusp' for p in pts):
            tag, ediff = 'cusp', sp.Rational(0)
        else:
            tag, ediff = 'orb2', sp.Rational(1, 2)
        pts.append(dict(lam=rt, t=t, tag=tag, ediff=ediff,
                        absl=complex(sp.N(sp.Abs(t), 40)).real,
                        tnum=complex(sp.N(t, 40))))
    pts.sort(key=lambda p: p['absl'])
    fold = pts[0]
    post = pts[1:]
    ceil16 = 16 * min(p['absl'] for p in post)
    return dict(label=label, xi=xi, family=fam, coeffs=co,
                chi=sp.srepr(chi), disc=int(disc),
                quad=str(quad.as_expr()),
                pts=pts, fold=fold, post=post, ceil16=ceil16)


def fmt(x, n=18):
    return sp.nstr if False else ("%.*g" % (n, x))


if __name__ == '__main__':
    out = []
    for row in ROWS:
        A = analyse(*row)
        out.append(A)
        co = A['coeffs']
        print("=" * 78)
        print("%s   xi = %s   [%s]" % (A['label'], A['xi'], A['family']))
        print("  (a,c,d,f,C) = (%d,%d,%d,%d,%d),  r = %d,  g = %d" %
              (co['a'], co['c'], co['d'], co['f'], co['C'], co['r'], co['g']))
        print("  P(n) = %dn^2 + %dn + %d ; Q(n) = %dn^2 + %dn + %d ; R(n) = %dn^2" %
              (co['a'], co['b'], co['c'], co['d'], co['e'], co['f'], co['g']))
        print("  chi(lam) = lam^3 - %d lam^2 + %d lam - %d ;  disc = %d" %
              (co['a'], co['d'], co['g'], A['disc']))
        print("  quadratic factor: %s" % A['quad'])
        print("  singular points (|t| increasing):")
        for i, p in enumerate(A['pts']):
            role = 'FOLD' if p is A['fold'] else 'outer'
            print("    t = %-28s = %-14.10f  |t| = %-12.9f  exp.diff %s  %-5s  %s"
                  % (sp.sstr(p['t']), p['tnum'].real, p['absl'], p['ediff'],
                     p['tag'], role))
        print("  plus t = 0 (log, I_n) and t = oo (log, delta_inf = 0)")
        print("  post-hypothesis finite singular set: %s"
              % ", ".join(sp.sstr(p['t']) + ("[%s]" % p['tag']) for p in A['post']))
        print("  monotonicity ceiling 16*min|t_post| = %.9f   log = %+.6f"
              % (A['ceil16'], __import__('math').log(A['ceil16'])))
        print("  score log|t_post,min| - 2 = %+.6f"
              % (__import__('math').log(min(p['absl'] for p in A['post'])) - 2))
    # machine-readable
    ser = []
    for A in out:
        ser.append(dict(label=A['label'], xi=A['xi'], family=A['family'],
                        coeffs={k: int(v) for k, v in A['coeffs'].items()},
                        disc=A['disc'],
                        pts=[dict(t=sp.sstr(p['t']), tre=p['tnum'].real,
                                  tim=p['tnum'].imag, absl=p['absl'],
                                  tag=p['tag'], ediff=str(p['ediff'])) for p in A['pts']],
                        ceil16=A['ceil16']))
    json.dump(ser, open(os.path.join(HERE, 'out', 'geometry.json'), 'w'), indent=1)
    print("\nwrote out/geometry.json")
