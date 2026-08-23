"""Q(x)-linear independence of the mu_4 host-A inventory (numerical rank check
mod a large prime, in the style of lattice/cdt_finder/indep_check2.py).

We test that no relation  sum_i P_i(x) f_i(x) = 0  with deg P_i <= D exists,
by computing the rank of the (m(D+1)) x Ntrunc matrix of the series x^j f_i.
This is a check of, not a substitute for, a proof.
"""
from fractions import Fraction as F
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mu4_series import build_x, mul, integ_dx, conditional_u

P = (1 << 61) - 1
N = 200


def tomod(ser):
    out = []
    for a in ser:
        num, den = a.numerator % P, a.denominator % P
        out.append(num*pow(den, P-2, P) % P)
    return out


def rank_mod(rows):
    rows = [r[:] for r in rows]
    n = len(rows); m = len(rows[0]); r = 0
    for c in range(m):
        piv = None
        for i in range(r, n):
            if rows[i][c]:
                piv = i; break
        if piv is None:
            continue
        rows[r], rows[piv] = rows[piv], rows[r]
        inv = pow(rows[r][c], P-2, P)
        rows[r] = [v*inv % P for v in rows[r]]
        for i in range(n):
            if i != r and rows[i][c]:
                f = rows[i][c]
                rows[i] = [(a - f*b) % P for a, b in zip(rows[i], rows[r])]
        r += 1
        if r == n:
            break
    return r


if __name__ == '__main__':
    Sx = build_x(N)
    names = list(Sx.keys())
    fx = [Sx[k] for k in names]
    # add J(-x^2) and the conditional function H(x^2) with a generic rational G
    from mu4_series import build_u
    Su = build_u(N)
    J = Su['J(-u) (CDT fifth fn)']
    Jx = [F(0)]*N
    for n in range(N//2):
        Jx[2*n] = J[n]
    fx.append(Jx); names.append('J(-x^2)')
    HH, A, g = conditional_u(N)
    a, b = 3, 7                       # generic "hypothesis" G = a/b
    Hc = [b*HH[n] - a*A[n] for n in range(N)]
    Hcx = [F(0)]*N
    for n in range(N//2):
        Hcx[2*n] = Hc[n]
    fx.append(Hcx); names.append('H_cond(x^2)')
    # its derivatives
    cur = Hcx
    for k in (1, 2, 3):
        cur = [cur[n+1]*(n+1) for n in range(len(cur)-1)] + [F(0)]
        fx.append(cur[:]); names.append(f'H_cond^({k})')
    # added integrations of H_cond
    prim = integ_dx([Hcx[n] for n in range(N)], N)
    fx.append(prim); names.append('int H_cond dx')

    m = len(fx)
    print(f'  m = {m} functions:', ', '.join(names))
    mods = [tomod(f) for f in fx]
    for D in range(0, 6):
        rows = []
        for f in mods:
            for j in range(D+1):
                rows.append([0]*j + f[:N-j])
        need = m*(D+1)
        if need > N:
            print(f'  deg P_i <= {D}: skipped (needs {need} > {N} coefficients)')
            continue
        r = rank_mod(rows)
        print(f'  deg P_i <= {D}:  rank {r} / {need}  ' + ('OK' if r == need else '*** RELATION ***'))
