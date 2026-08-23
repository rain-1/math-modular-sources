"""EMN projector: core series/arithmetic engine.

Notation (the note's S, H renamed to avoid clashes):
  ang   = the trigonometric uniformiser  S  (z = 1 - cos ang)
  Hfun  = the note's H(z) = z*Phi(z),  Hfun(1-cos ang) = script-H(ang)
  We never use PARI builtin names.

Master ODE:   z(2-z) F'' + (1-z) F' = Psi(z),    Psi = d^2F/d ang^2 expressed in z.
For Hfun: Psi = 1/(2(1-z)).
Series recursion from  n(2n-1) f_n = [z^{n-1}] Psi  +  (n-1)^2 f_{n-1}.
"""
from fractions import Fraction as Fr
import math, json, os

# ---------------------------------------------------------------- series solve
def solve_from_psi(psi_coeffs, N, f0=Fr(0)):
    """psi_coeffs[j] = [z^j] Psi(z), j=0..N.  Returns f[0..N] with f[0]=f0."""
    f = [Fr(0)]*(N+1)
    f[0] = f0
    prev = Fr(0)
    for n in range(1, N+1):
        f[n] = (psi_coeffs[n-1] + Fr((n-1)**2)*prev)/Fr(n*(2*n-1))
        prev = f[n]
    return f

def series_H(N):
    """The EMN generator H(z): Psi = 1/(2(1-z))."""
    psi = [Fr(1,2)]*(N+1)
    return solve_from_psi(psi, N)

# ------------------------------------------------------- rational-function tools
def ratseries(num, den, N):
    """Taylor coefficients of num(z)/den(z) to order N; num,den lists of Fr."""
    assert den[0] != 0
    c = [Fr(0)]*(N+1)
    for n in range(N+1):
        s = num[n] if n < len(num) else Fr(0)
        for k in range(1, min(n, len(den)-1)+1):
            s -= den[k]*c[n-k]
        c[n] = s/den[0]
    return c

def polymul(a, b):
    r = [Fr(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] += x*y
    return r

def chebT_in_z(a):
    """T_a(1-z) as a polynomial in z (list of Fr, low->high)."""
    x = [Fr(1), Fr(-1)]          # 1 - z
    T0, T1 = [Fr(1)], x[:]
    if a == 0: return T0
    if a == 1: return T1
    for _ in range(2, a+1):
        T2 = [2*c for c in polymul(x, T1)]
        for i, c in enumerate(T0):
            T2[i] -= c
        T0, T1 = T1, T2
    return T1

def psi_basis(a, kind):
    """Psi for the symmetric pair  h_{a,b} = script-H(a*ang+b)+script-H(a*ang-b),
       kind='0'   : b=0     -> Psi = a^2 / T_a(1-z)          (single script-H, x2)
       kind='pi3' : b=pi/3  -> Psi = (a^2/2) T_a/(T_a^2-3/4)
       Returns (num, den) polynomials in z.
       NOTE for kind='0' the pair is 2*script-H(a*ang); we use the *single*
       function script-H(a*ang), Psi = (a^2/2)/T_a(1-z)."""
    T = chebT_in_z(a)
    if kind == '0':
        return ([Fr(a*a, 2)], T)
    elif kind == 'pi3':
        num = [Fr(a*a, 2)*c for c in T]
        den = polymul(T, T)
        den[0] -= Fr(3, 4)
        return (num, den)
    raise ValueError(kind)

# --------------------------------------------------------------- measurements
def den_profile(coeffs, ns):
    out = {}
    for n in ns:
        d = coeffs[n].denominator
        out[n] = (math.log(d)/n if d > 1 else 0.0)
    return out

def padic_slope(coeffs, p, ns):
    def vp(fr):
        num, den = fr.numerator, fr.denominator
        if num == 0: return None
        v = 0
        while num % p == 0: num//=p; v+=1
        while den % p == 0: den//=p; v-=1
        return v
    return {n: (vp(coeffs[n]), (vp(coeffs[n])/n if vp(coeffs[n]) is not None else None)) for n in ns}

if __name__ == '__main__':
    N = 2000
    H = series_H(N)
    print("h_1..h_6:", [str(H[n]) for n in range(1,7)])
    ns = [50,100,200,400,800,1200,1600,2000]
    print("log den(h_n)/n :", {n: round(v,5) for n,v in den_profile(H, ns).items()},
          "  2log2 =", round(2*math.log(2),5))
    for p in (2,3,5):
        pr = padic_slope(H, p, ns)
        print(f"v_{p}(h_n): ", {n: pr[n][0] for n in ns})
