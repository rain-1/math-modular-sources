"""
chi3_forms.py -- exact linear forms in 1 and L(2,chi_{-3}) from rational functions
with double poles at thirds.

Framework (see ONE_CLASS_TWO_WORLDS.md Sec 2):
  R(t) = N(t) / prod_i (t-p_i)^2 ,  p_i in {-(j+1/3)} u {-(j+2/3)},  deg N <= deg D - 2.
  S = sum_{t>=0} R(t)
    = 4*(SA1+SA2)*zeta(2) + (9/2)*(SA1-SA2)*L + (pi/sqrt3)*SB1 + rational
  where SA1 = sum of double-pole coefficients over the 1/3-class, etc.
  Antisymmetry R(-t-(m+1)) = -R(t) forces SA2=-SA1 and SB1=SB2=0, hence
    S = 9*SA1*L(2,chi_{-3}) - (rational).
"""
from fractions import Fraction as F

# ---------- tiny exact polynomial layer (coeff lists, low degree first) ----------
def pmul(a, b):
    r = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i+j] += x*y
    return r

def peval(p, x):
    s = F(0)
    for c in reversed(p):
        s = s*x + c
    return s

def pder(p):
    return [p[i]*i for i in range(1, len(p))] or [F(0)]

def linear(a, b):   # a*t + b
    return [F(b), F(a)]

# ---------- the linear form ----------
def linform(poles, N):
    """poles: list of Fractions p_i (the actual pole locations, R ~ A/(t-p)^2).
       N: numerator polynomial (coeff list).
       Returns dict with SA1,SA2,SB1,SB2 and the rational constant, so that
       sum_{t>=0} R(t) = 4(SA1+SA2) zeta2 + (9/2)(SA1-SA2) L + (pi/sqrt3) SB1 + Crat."""
    n = len(poles)
    Np = pder(N)
    A = [F(0)]*n
    B = [F(0)]*n
    for i, p in enumerate(poles):
        pr = F(1)
        for j, q in enumerate(poles):
            if j != i:
                pr *= (p-q)**2
        Ni = peval(N, p)
        A[i] = Ni/pr
        s = F(0)
        for j, q in enumerate(poles):
            if j != i:
                s += F(1)/(p-q)
        B[i] = (peval(Np, p) - 2*s*Ni)/pr
    assert sum(B) == 0, ("sum B != 0", sum(B))
    SA1 = SA2 = SB1 = SB2 = F(0)
    Crat = F(0)
    for i, p in enumerate(poles):
        x = -p                      # pole of R at t = p, i.e. R ~ A/(t+x)^2 with x=-p
        # x = j + 1/3  or  j + 2/3 ,  j = integer >= 0
        j = int(x - (x - int(x)))
        frac = x - j
        assert frac in (F(1,3), F(2,3)), (x, frac)
        # zeta(2,x) = zeta(2,frac) - sum_{l<j} 1/(l+frac)^2
        # -psi(x)   = -psi(frac) - sum_{l<j} 1/(l+frac)
        tail2 = sum(F(1)/(l+frac)**2 for l in range(j))
        tail1 = sum(F(1)/(l+frac) for l in range(j))
        if frac == F(1,3):
            SA1 += A[i]; SB1 += B[i]
        else:
            SA2 += A[i]; SB2 += B[i]
        Crat += -A[i]*tail2 - B[i]*tail1
    return dict(SA1=SA1, SA2=SA2, SB1=SB1, SB2=SB2, Crat=Crat,
                zeta2=4*(SA1+SA2), L=F(9,2)*(SA1-SA2), pi_over_sqrt3=SB1)

def poles_thirds(m):
    """the 2m+2 poles t = -(j+1/3), -(j+2/3), j=0..m"""
    return [F(-(3*j+1), 3) for j in range(m+1)] + [F(-(3*j+2), 3) for j in range(m+1)]
