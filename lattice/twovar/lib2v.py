"""Shared library: one-variable rows, companions, and two-variable decouplings."""
from fractions import Fraction as F
from math import comb, gcd
from functools import lru_cache

# ---------- one-variable rows and companions ----------

def apery3(N):
    """Apery zeta(3) numbers a_n = sum_k C(n,k)^2 C(n+k,k)^2 and companion b_n
    (normalisation b_0=0, b_1=6, i.e. b_n/a_n -> zeta(3))."""
    a = [sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)) for n in range(N+1)]
    b = [F(0), F(6)]
    # (n+1)^3 u_{n+1} = (34n^3+51n^2+27n+5) u_n - n^3 u_{n-1}
    for n in range(1, N):
        b.append(((34*n**3+51*n**2+27*n+5)*b[n] - n**3*b[n-1])/F((n+1)**3))
    return a, b[:N+1]

def apery2(N):
    """zeta(2) row a_n = sum_k C(n,k)^2 C(n+k,k); companion b_n/a_n -> zeta(2).
    (n+1)^2 u_{n+1} = (11n^2+11n+3) u_n + n^2 u_{n-1}."""
    a = [sum(comb(n,k)**2*comb(n+k,k) for k in range(n+1)) for n in range(N+1)]
    b = [F(0), F(5)]
    for n in range(1, N):
        b.append(((11*n**2+11*n+3)*b[n] + n**2*b[n-1])/F((n+1)**2))
    return a, b[:N+1]

def rowE(N):
    """Zagier row E (12,4,32): a_n = sum_k C(n,k) C(2k,k) C(2n-2k,n-k);
    (n+1)^2 u_{n+1} = (12n^2+12n+4) u_n - 32 n^2 u_{n-1}. b_n/a_n -> G/2 (Catalan)."""
    a = [sum(comb(n,k)*comb(2*k,k)*comb(2*n-2*k,n-k) for k in range(n+1)) for n in range(N+1)]
    b = [F(0), F(1)]
    for n in range(1, N):
        b.append(((12*n**2+12*n+4)*b[n] - 32*n**2*b[n-1])/F((n+1)**2))
    return a, b[:N+1]

def cooper10(N):
    """Cooper s_10: a_n = sum_k C(n,k)^4;  (a,b,c,d)=(6,2,-64,4) in the
    4-term normalisation (n+1)^3 u_{n+1} = (2n+1)(a n^2+a n+b) u_n
      + n(c n^2 + d) u_{n-1}."""
    a = [sum(comb(n,k)**4 for k in range(n+1)) for n in range(N+1)]
    b = [F(0), F(1)]
    for n in range(1, N):
        b.append((2*(2*n+1)*(3*n*n+3*n+1)*b[n] + 4*n*(16*n*n-1)*b[n-1])/F((n+1)**3))
    return a, b[:N+1]

# ---------- two-variable decouplings ----------
# each returns c(a,b) as an exact integer

def z3_D1(a,b):    # sum_k C(a,k)C(b,k)C(a+k,k)C(b+k,k)      -- "balanced"
    return sum(comb(a,k)*comb(b,k)*comb(a+k,k)*comb(b+k,k) for k in range(min(a,b)+1))

def z3_D2(a,b):    # sum_k C(a,k)^2 C(b+k,k)^2
    return sum(comb(a,k)**2*comb(b+k,k)**2 for k in range(a+1))

def z3_D3(a,b):    # sum_k C(a,k)^2 C(a+k,k) C(b+k,k)
    return sum(comb(a,k)**2*comb(a+k,k)*comb(b+k,k) for k in range(a+1))

def z3_D4(a,b):    # sum_k C(a,k)C(b,k)C(a+k,k)^2
    return sum(comb(a,k)*comb(b,k)*comb(a+k,k)**2 for k in range(min(a,b)+1))

def z3_D5(a,b):    # sum_k C(2k,k)^2 C(a+k,2k) C(b+k,2k)   == D1 (identity check)
    return sum(comb(2*k,k)**2*comb(a+k,2*k)*comb(b+k,2*k) for k in range(min(a,b)+1))

def z2_D1(a,b):    # sum_k C(a,k) C(2k,k) C(b+k,2k)   -- diagonal = zeta(2) row
    return sum(comb(a,k)*comb(2*k,k)*comb(b+k,2*k) for k in range(min(a,b)+1))

def z2_D2(a,b):    # sum_k C(a,k)C(b,k)C(a+k,k)
    return sum(comb(a,k)*comb(b,k)*comb(a+k,k) for k in range(min(a,b)+1))

def z2_D3(a,b):    # sum_k C(a,k)^2 C(b+k,k)
    return sum(comb(a,k)**2*comb(b+k,k) for k in range(a+1))

def E_D1(a,b):     # sum_k C(a,k) C(2k,k) C(2b-2k,b-k)
    return sum(comb(a,k)*comb(2*k,k)*comb(2*b-2*k,b-k) for k in range(min(a,b)+1))

def E_D2(a,b):     # sum_k C(2k,k)^2 C(a-k+?)  -- via a_n = sum C(n,k)C(2k,k)C(2n-2k,n-k)
    # decouple the two "n"s in C(2n-2k,n-k) is impossible; use C(a,k) and C(2b-2k,b-k)
    return None

def s10_D1(a,b):   # sum_k C(a,k)^2 C(b,k)^2
    return sum(comb(a,k)**2*comb(b,k)**2 for k in range(min(a,b)+1))

def s10_D2(a,b):   # sum_k C(a,k)^3 C(b,k)
    return sum(comb(a,k)**3*comb(b,k) for k in range(min(a,b)+1))

def lcmrange(n):
    r = 1
    for i in range(1, n+1):
        r = r*i//gcd(r,i)
    return r
