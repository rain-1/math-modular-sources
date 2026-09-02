"""K = Q(sqrt5) as pairs (u,v) of Fractions meaning u + v*sqrt5.
Z[phi] = { (a+b sqrt5)/2 : a,b in Z, a = b mod 2 }."""
from fractions import Fraction as F
from math import gcd
def K(u=0, v=0): return (F(u), F(v))
def kadd(a, b): return (a[0]+b[0], a[1]+b[1])
def ksub(a, b): return (a[0]-b[0], a[1]-b[1])
def kmul(a, b): return (a[0]*b[0]+5*a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def kscal(c, a): return (c*a[0], c*a[1])
def kinv(a):
    n = a[0]*a[0]-5*a[1]*a[1]
    return (a[0]/n, -a[1]/n)
ONE = K(1); ZERO = K(0)
PHI = (F(1, 2), F(1, 2))
PHI5 = (F(11, 2), F(5, 2))          # phi^5 = (11+5 sqrt5)/2
S    = (F(-11, 2), F(-5, 2))        # s = -phi^5
SINV = (F(11, 2), F(-5, 2))         # s^{-1} = -phi^{-5} = (11-5 sqrt5)/2
def isZphi(a):
    u, v = a
    if (2*u).denominator != 1 or (2*v).denominator != 1: return False
    return ((2*u)-(2*v)) % 2 == 0
def kden(a):
    """least positive integer d with d*a in Z[phi];  a = (U+V sqrt5)/2, U=2u, V=2v."""
    from math import lcm
    U, V = 2*a[0], 2*a[1]
    d0 = lcm(U.denominator, V.denominator)
    return d0 if (d0*U - d0*V) % 2 == 0 else 2*d0
