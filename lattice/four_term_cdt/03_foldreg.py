#!/usr/bin/env python3
"""03_foldreg.py -- the fold-regularity test on the six mixed-exponent hosts.

Under the hypothesis xi = lim b_n/a_n in Q, the conditional function
H = B - xi A must be single-valued and holomorphic at the fold t_fold.  The
operational test (CATALAN_THREE_PERIOD.md 1) is RADIUS ENLARGEMENT: the
coefficients h_n = b_n - xi a_n must satisfy |h_n|^{1/n} -> 1/|t_2^post|, the
reciprocal of the SECOND singular point, instead of 1/|t_fold|.

We also report the residual growth rate of b_n - xi' a_n for xi' a nearby wrong
constant (control: it must revert to 1/|t_fold|), and the enlargement factor
log|t_2^post / t_fold| -- the number of nats the hypothesis buys.
"""
import os, sys, json
from fractions import Fraction
from mpmath import mp, mpf, mpc, catalan, zeta, nstr, log as mlog

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
mp.dps = 260

import importlib.util
spec = importlib.util.spec_from_file_location("d4", os.path.join(HERE, "04_denom.py"))
d4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d4)

NMAX = int(os.environ.get('NMAX', '400'))


def Lchi3():
    """L(2, chi_{-3}) = sum chi_{-3}(n)/n^2 via the Hurwitz zeta."""
    return (mp.zeta(2, mpf(1)/3) - mp.zeta(2, mpf(2)/3)) / 9


def xis():
    G = mp.catalan
    z2 = mp.zeta(2)
    L = Lchi3()
    return {
        "R1": G/4,
        "R2": G/2 - 3*z2/16,
        "R3": 3*z2/8 - G/2,
        "R4": (2*z2 + 15*L)/32,
        "R5": (15*L - 6*z2)/16,
        "R6": (2*z2 - 3*L)/8,
    }


def main():
    geo = json.load(open(os.path.join(HERE, 'out', 'geometry.json')))
    G = {g['label']: g for g in geo}
    XI = xis()
    out = []
    for lab, r, a, c, d, f, C in d4.ROWS:
        if lab not in XI:
            continue
        co = d4.coeffs(r, a, c, d, f, C)
        A = d4.seq(co, NMAX, 1)
        B = d4.seq(co, NMAX, 0, u1=1)
        # the Apery limit in THIS normalisation (b_0=0, b_1=1)
        rat = mpf(B[NMAX].numerator)/mpf(B[NMAX].denominator) / \
              (mpf(A[NMAX].numerator)/mpf(A[NMAX].denominator))
        xi_lit = XI[lab]
        # the scan's xi may differ by the normalisation of b_1: find the scalar
        scal = rat / xi_lit
        xi_use = xi_lit * scal
        pts = G[lab]['pts']
        tf = abs(pts[0]['absl']); t2 = abs(pts[1]['absl'])
        rows = []
        for n in (NMAX//4, NMAX//2, 3*NMAX//4, NMAX):
            an = mpf(A[n].numerator)/mpf(A[n].denominator)
            bn = mpf(B[n].numerator)/mpf(B[n].denominator)
            hn = bn - xi_use*an
            rows.append((n, float(mp.log(abs(hn))/n), float(mp.log(abs(an))/n)))
        # control: a wrong constant
        wr = xi_use*(1 + mpf(10)**-8)
        ctrl = []
        for n in (NMAX//2, NMAX):
            an = mpf(A[n].numerator)/mpf(A[n].denominator)
            bn = mpf(B[n].numerator)/mpf(B[n].denominator)
            ctrl.append((n, float(mp.log(abs(bn - wr*an))/n)))
        rec = dict(label=lab, scal=float(scal), t_fold=tf, t_2post=t2,
                   growth=rows, control=ctrl,
                   expect_H=float(-mp.log(t2)), expect_A=float(-mp.log(tf)),
                   gain=float(mp.log(t2/tf)))
        out.append(rec)
        print("%s  b_1 normalisation factor xi_scan/xi_here = %s" % (lab, nstr(scal, 12)))
        print("   |t_fold| = %.10f   |t_2^post| = %.10f   enlargement %+.6f nats"
              % (tf, t2, rec['gain']))
        print("   n      log|h_n|/n     log|a_n|/n     (targets: %+.6f and %+.6f)"
              % (rec['expect_H'], rec['expect_A']))
        for n, lh, la in rows:
            print("   %-6d %+.8f    %+.8f" % (n, lh, la))
        print("   control (xi perturbed by 1e-8): %s   -> must revert to %+.6f"
              % (["%d: %+.6f" % x for x in ctrl], rec['expect_A']))
    json.dump(out, open(os.path.join(HERE, 'out', 'foldreg.json'), 'w'), indent=1)
    print("\nwrote out/foldreg.json")


if __name__ == '__main__':
    main()
