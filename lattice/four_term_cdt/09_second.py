#!/usr/bin/env python3
"""09_second.py -- is there a SECOND conditional direction on these hosts?

The Picard-Fuchs operator L has order 2; the companion B solves L B = t.  The
solutions of L y = t^k (k = 2,3,...) are further power series in Q[[x]],
holomorphic at 0, each with its own fold constant xi^(k) = lim b^(k)_n / a_n.
If the xi^(k) span more than a line in Q + Q zeta(2) + Q G then a Q-linear
RELATION among 1, zeta(2), G supplies a fold-regular combination even when no
single xi is rational -- this is exactly CDT's H_A, H_B, H_C mechanism, and it is
what decides whether the THREE-PERIOD hypothesis buys anything here.

In recurrence terms y^(k) is the solution with b_0 = ... = b_{k-1} = 0, b_k = 1
(the recurrence is violated exactly at n = k-1).
"""
import os, sys, json
from fractions import Fraction
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location("d4", os.path.join(HERE, "04_denom.py"))
d4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(d4)
from mpmath import mp, mpf, pslq, nstr
mp.dps = 220

NMAX = int(os.environ.get('NMAX', '3000'))


def seq_k(co, N, k):
    """b_0=..=b_{k-1}=0, b_k=1, then the four-term recurrence for n >= k."""
    P = lambda n: co['a']*n*n + co['b']*n + co['c']
    Q = lambda n: co['d']*n*n + co['e']*n + co['f']
    R = lambda n: co['g']*n*n + co['h']*n + co['j']
    u = {-2: mpf(0), -1: mpf(0)}
    for i in range(0, k):
        u[i] = mpf(0)
    u[k] = mpf(1)
    for n in range(k, N):
        u[n+1] = (P(n)*u[n] - Q(n)*u[n-1] + R(n)*u[n-2]) / mpf((n+1)**2)
    return u


def seq_a(co, N):
    P = lambda n: co['a']*n*n + co['b']*n + co['c']
    Q = lambda n: co['d']*n*n + co['e']*n + co['f']
    R = lambda n: co['g']*n*n + co['h']*n + co['j']
    u = {-2: mpf(0), -1: mpf(0), 0: mpf(1)}
    for n in range(0, N):
        u[n+1] = (P(n)*u[n] - Q(n)*u[n-1] + R(n)*u[n-2]) / mpf((n+1)**2)
    return u


def main():
    G = mp.catalan; z2 = mp.zeta(2)
    L3 = (mp.zeta(2, mpf(1)/3) - mp.zeta(2, mpf(2)/3))/9
    basis = {"sqrt2": [mpf(1), z2, G], "sqrt3": [mpf(1), z2, L3]}
    names = {"sqrt2": ["1", "zeta(2)", "G"], "sqrt3": ["1", "zeta(2)", "L(2,chi_-3)"]}
    fam = {"R1": "sqrt2", "R2": "sqrt2", "R3": "sqrt2",
           "R4": "sqrt3", "R5": "sqrt3", "R6": "sqrt3"}
    out = []
    for lab, r, a, c, d, f, C in d4.ROWS:
        if lab not in fam:
            continue
        co = d4.coeffs(r, a, c, d, f, C)
        A = seq_a(co, NMAX)
        print("%s" % lab)
        rec = dict(label=lab, xis=[])
        for k in (1, 2, 3):
            Bk = seq_k(co, NMAX, k)
            xi = Bk[NMAX]/A[NMAX]
            xi2 = Bk[NMAX-200]/A[NMAX-200]
            agree = int(-mp.log10(abs(xi - xi2)/abs(xi))) if xi != xi2 else 999
            v = pslq([xi] + basis[fam[lab]], tol=mpf(10)**(-min(agree, 60)+5), maxcoeff=10**12,
                     maxsteps=10**5)
            print("   xi^(%d) = %s   (agree ~%d digits)" % (k, nstr(xi, 30), agree))
            if v:
                num = "-(" + " + ".join("%s*%s" % (v[i+1], names[fam[lab]][i])
                                        for i in range(3) if v[i+1]) + ")/%s" % v[0]
                print("        pslq: %s  ->  xi^(%d) = %s" % (v, k, num))
            else:
                print("        pslq: no relation with 1, zeta(2), %s at height 1e12"
                      % names[fam[lab]][2])
            rec['xis'].append(dict(k=k, xi=str(nstr(xi, 40)), agree=agree,
                                   relation=[int(x) for x in v] if v else None))
        out.append(rec)
    json.dump(out, open(os.path.join(HERE, 'out', os.environ.get('OUTJSON', 'second.json')), 'w'), indent=1)
    print("\nwrote out/%s" % os.environ.get("OUTJSON", "second.json"))


if __name__ == '__main__':
    main()
