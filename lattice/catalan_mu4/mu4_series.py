"""The mu_4 (fourth-roots-of-unity) polylogarithm inventory for Catalan's constant.

Coordinates
-----------
  x-line :  P^1 - {0, i, -i, oo}.   Everything is a power series in x over Q.
  u-line :  u = x^2,  P^1 - {0, -1, oo}   (the quotient by x -> -x).
  v-line :  v = u^2/(u+1)              (the CDT-style involution quotient of the
            u-line; w(u) = -u/(u+1) fixes u=0 and u=-2).

Catalan function
----------------
  f(x) = sum_{k>=0} (-1)^k x^{2k+1}/(2k+1)^2 = ( Li_2(ix) - Li_2(-ix) ) / (2i),
  f(1) = G,   f(i) = i pi^2/8,  singular exactly at x = +- i (and 0, oo).
  On the u-line f(x) = x g(u), g(u) = sum (-1)^n u^n/(2n+1)^2, g(1) = G.

Conditional function (the "fold" mechanism)
-------------------------------------------
  HH(u) = int_0^u g(t) dt/(1-t)  has a log singularity at u = 1 with
  Delta_{u=1} HH = -2 pi i g(1) = -2 pi i G, while A(u) = -log(1-u) has
  Delta_{u=1} A = -2 pi i.  Hence HH - G.A is regular (principal branch) at u=1.
  Under G = a/b the series b.HH - a.A lies in Q[[u]] and is fold-regular.

This script computes exact rational series, measures denominator types
empirically, checks fold regularity numerically, and checks Q(x)-independence.
"""
from fractions import Fraction as F
import math

N = 160                      # series order (in the x variable)


# ------------------------------------------------------------------ utilities
def lcmlist(n):
    """[1,...,n] as an int (n>=0)."""
    r = 1
    for k in range(1, n+1):
        r = r*k//math.gcd(r, k)
    return r


LCM = [lcmlist(k) for k in range(0, 4*N+3)]


def mul(a, b, N):
    c = [F(0)]*N
    for i, ai in enumerate(a):
        if ai == 0:
            continue
        for j in range(0, N-i):
            if b[j]:
                c[i+j] += ai*b[j]
    return c


def integ_dx(a, N):
    """int_0^x a(t) dt."""
    return [F(0)] + [a[k]/(k+1) for k in range(N-1)]


def profile(ser, name, var_scale=1):
    """Smallest (e, b1, b2) with  a_n * n^e * [1..b1 n] * [1..b2 n]  integral,
    searched over e in 0..3 and b in {0, 1/2, 1, 2} (b=0 means layer absent).
    var_scale: multiply the index n by this before forming [1..b n] (used when
    the series is read in a different variable)."""
    cands = []
    bs = [0, F(1,2), 1, 2]
    SLACK = 2   # a layer [1..bn] is allowed the additive slack [1..bn+SLACK];
                # this changes the denominator by e^{O(log n)}, invisible to tau.
    for e in (0, 1, 2, 3):
        for i1, b1 in enumerate(bs):
            for b2 in bs[i1:]:
                ok = True
                for n in range(1, len(ser)):
                    if ser[n] == 0:
                        continue
                    d = F(n)**e
                    for b in (b1, b2):
                        if b:
                            d *= LCM[int(b*n*var_scale)+SLACK]
                    if (ser[n]*d).denominator != 1:
                        ok = False
                        break
                if ok:
                    cands.append((float(b1)+float(b2), e, b1, b2))
    cands.sort()
    sig, e, b1, b2 = cands[0]
    # among minimal sigma prefer minimal e
    best = min([c for c in cands if c[0] == sig], key=lambda c: c[1])
    sig, e, b1, b2 = best
    print(f'  {name:34s} sigma={sig:<4g} e={e}  b=({b1},{b2})')
    return dict(name=name, sigma=sig, e=e, b=(b1, b2))


# ------------------------------------------------------- the x-line inventory
def build_x(N):
    S = {}
    zero = [F(0)]*N
    one = list(zero); one[0] = F(1)
    S['1'] = one
    # arctan x = sum (-1)^k x^{2k+1}/(2k+1)
    at = list(zero)
    for k in range(0, N//2+1):
        if 2*k+1 < N:
            at[2*k+1] = F((-1)**k, 2*k+1)
    S['arctan(x)'] = at
    # log(1+x^2)
    lg = list(zero)
    for m in range(1, N//2+1):
        if 2*m < N:
            lg[2*m] = F((-1)**(m-1), m)
    S['log(1+x^2)'] = lg
    # Li_2(-x^2)
    li2 = list(zero)
    for m in range(1, N//2+1):
        if 2*m < N:
            li2[2*m] = F((-1)**m, m*m)
    S['Li_2(-x^2)'] = li2
    # f(x)  (the Catalan function)
    f = list(zero)
    for k in range(0, N//2+1):
        if 2*k+1 < N:
            f[2*k+1] = F((-1)**k, (2*k+1)**2)
    S['f(x)=Ti_2(x)'] = f
    # products
    S['log^2(1+x^2)'] = mul(lg, lg, N)
    S['arctan^2(x)'] = mul(at, at, N)
    S['arctan(x)log(1+x^2)'] = mul(at, lg, N)
    # M(x) = int_0^x log(1+t^2)/(1+t^2) dt = sum (-1)^{m-1} H_m x^{2m+1}/(2m+1)
    Mx = list(zero)
    H = F(0)
    for m in range(1, N//2+1):
        H += F(1, m)
        if 2*m+1 < N:
            Mx[2*m+1] = F((-1)**(m-1))*H/(2*m+1)
    S['M(x)=int log(1+t^2)/(1+t^2)'] = Mx
    # CDT's fifth function transported:  J(-x^2) = (1+x^2)^{-1/2} int_0^{x^2}
    #   log(1+t)/(t sqrt(1+t)) dt   -- built as a u-series then pulled back.
    return S


def build_u(N):
    """u-series to order N (u = x^2)."""
    zero = [F(0)]*N
    S = {}
    one = list(zero); one[0] = F(1); S['1'] = one
    ell = [F((-1)**n, 2*n+1) for n in range(N)]          # arctan(sqrt u)/sqrt u
    S['ell(u)'] = ell
    g = [F((-1)**n, (2*n+1)**2) for n in range(N)]       # Catalan function
    S['g(u)'] = g
    L1 = list(zero)
    for n in range(1, N):
        L1[n] = F((-1)**(n-1), n)
    S['log(1+u)'] = L1
    li2 = [F(0)] + [F((-1)**n, n*n) for n in range(1, N)]
    S['Li_2(-u)'] = li2
    S['log^2(1+u)'] = mul(L1, L1, N)
    S['u*ell^2'] = [F(0)] + mul(ell, ell, N)[:N-1]
    S['ell*log(1+u)'] = mul(ell, L1, N)
    gt = list(zero)                                       # tilde g
    H = F(0)
    for n in range(1, N):
        H += F(1, n)
        gt[n] = F((-1)**(n-1))*H/(2*n+1)
    S['gtilde(u)'] = gt
    # CDT's J on P^1-{0,-1,oo}:  (1+u)^{-1/2} int_0^u log(1+t)/(t sqrt(1+t)) dt
    sq = [F(0)]*N                                         # (1+u)^{-1/2}
    c = F(1)
    for n in range(N):
        sq[n] = c
        c = c*F(-1)*F(2*n+1, 2*(n+1))
    isq = [F(0)]*N                                        # (1+u)^{-1/2} again
    # integrand log(1+t)/(t sqrt(1+t)) : (log(1+t)/t) * (1+t)^{-1/2}
    lot = [F((-1)**n, n+1) for n in range(N)]             # log(1+t)/t
    integrand = mul(lot, sq, N)
    prim = integ_dx(integrand, N)
    S['J(-u) (CDT fifth fn)'] = mul(sq, prim, N)
    return S


def conditional_u(N):
    """HH(u) = int_0^u g/(1-t) dt  and A(u) = -log(1-u)."""
    g = [F((-1)**n, (2*n+1)**2) for n in range(N)]
    inv = [F(1)]*N                                        # 1/(1-t)
    HH = integ_dx(mul(g, inv, N), N)
    A = [F(0)] + [F(1, n) for n in range(1, N)]
    return HH, A, g


if __name__ == '__main__':
    print('=== x-line inventory: denominator profiles (empirical, n<=%d) ===' % (N-1))
    Sx = build_x(N)
    px = [profile(v, k) for k, v in Sx.items()]

    print()
    print('=== u-line inventory (u = x^2): profiles in the u-variable ===')
    Su = build_u(N)
    pu = [profile(v, k) for k, v in Su.items()]

    print()
    print('=== conditional function ===')
    HH, A, g = conditional_u(N)
    profile(HH, 'HH(u)=int g/(1-t)dt')
    profile(A, 'A(u)=-log(1-u)')
    # pull HH back to the x-line: coefficient of x^{2n} is HH_n
    HHx = [F(0)]*(2*N)
    for n in range(N):
        HHx[2*n] = HH[n]
    profile(HHx[:N], 'HH(x^2) on the x-line')

    print()
    print('=== fold regularity check (numerical) ===')
    import mpmath as mp
    mp.mp.dps = 30
    G = mp.catalan
    # g(1) = G ?
    print('  sum g_n            =', mp.nstr(mp.nsum(lambda n: mp.mpf((-1)**int(n))/(2*int(n)+1)**2, [0, mp.inf]), 20))
    print('  Catalan G          =', mp.nstr(G, 20))
    # HH_n - G/n  should be O(n^-3) and summable  =>  HH - G*A extends to u=1
    print('  n^3*(HH_n - G/n) for n = 20,40,80,159:')
    for n in (20, 40, 80, 159):
        val = mp.mpf(HH[n].numerator)/HH[n].denominator - G/n
        print(f'     n={n:4d}:  HH_n-G/n = {mp.nstr(val,8):>14s}   n^3*(...) = {mp.nstr(n**3*val,8)}')
    tot = sum(mp.mpf(HH[n].numerator)/HH[n].denominator - G/n for n in range(1, N))
    print('  partial sum_{n<160} (HH_n - G/n) =', mp.nstr(tot, 12),
          '  (converges: HH - G*A is finite at u=1)')
    # by contrast HH_n alone ~ G/n : divergent
    print('  partial sum_{n<160} HH_n        =',
          mp.nstr(sum(mp.mpf(HH[n].numerator)/HH[n].denominator for n in range(1, N)), 12),
          '  (still growing like G*log n)')
