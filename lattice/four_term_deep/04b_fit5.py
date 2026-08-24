#!/usr/bin/env python3
"""Validation of the five-term framework: build a five-term row as the signed
binomial transform (generic cusp move) of a known four-term row, fit the
recurrence, and check it against the class parametrisation of 04_fscan5.c."""
from fractions import Fraction as Fr
import itertools, sys

def four_term(P,Q,R,N):
    u=[Fr(0),Fr(0),Fr(1)]
    for n in range(N):
        u.append((P(n)*u[n+2]-Q(n)*u[n+1]+R(n)*u[n])/Fr((n+1)**2))
    return u[2:]

def binom_transform(u,nu):
    from math import comb
    return [sum(Fr(comb(n,m))*Fr(-nu)**(n-m)*u[m] for m in range(n+1)) for n in range(len(u))]

def fit5(v,N=26):
    # (n+1)^2 v_{n+1} = P v_n - Q v_{n-1} + R v_{n-2} - T v_{n-3}, 12 unknowns
    import fractions
    rows=[];rhs=[]
    for n in range(3,N):
        r=[]
        for (poly,sgn,sh) in [(0,1,0),(1,-1,1),(2,1,2),(3,-1,3)]:
            for p in (n*n,n,1): r.append(Fr(sgn)*Fr(p)*v[n-sh])
        rows.append(r); rhs.append(Fr((n+1)**2)*v[n+1])
    # gaussian elimination
    m=len(rows); ncol=12
    A=[rows[i][:]+[rhs[i]] for i in range(m)]
    piv=[]; r=0
    for cidx in range(ncol):
        pr=None
        for i in range(r,m):
            if A[i][cidx]!=0: pr=i;break
        if pr is None: continue
        A[r],A[pr]=A[pr],A[r]
        pvt=A[r][cidx]
        A[r]=[x/pvt for x in A[r]]
        for i in range(m):
            if i!=r and A[i][cidx]!=0:
                fac=A[i][cidx]; A[i]=[A[i][k]-fac*A[r][k] for k in range(ncol+1)]
        piv.append(cidx); r+=1
    free=[c for c in range(ncol) if c not in piv]
    sol=[Fr(0)]*ncol
    for i,cidx in enumerate(piv): sol[cidx]=A[i][ncol]
    return sol, free

# --- the elliptic-K3 four-term row -------------------------------------------
P=lambda n:11*n*n+11*n+4; Q=lambda n:37*n*n+3; R=lambda n:27*n*n-27*n+6
u=four_term(P,Q,R,30)
print("u:", [int(x) for x in u[:9]])
for nu in (1,2,-1):
    v=binom_transform(u,nu)
    sol,free=fit5(v)
    print(f"nu={nu}  free={free}")
    print("   fitted (a,b,c, d,e,f, g,h,j, k,l,m) =", [str(x) for x in sol])
    a,b,c,d,e,f,g,h,j,k,l,m = sol
    if k!=0:
        print("   rho from b/a:", 1-b/a if a else None, " e/-2d:", -e/(2*d) if d else None,
              " h: -(1+3rho)g ->", -h/g if g else None, "  l/k:", -l/k if k else None)
        print("   s1+s2 = -l/k =", -l/k, "  s1*s2 = m/k =", m/k, " (need s1+s2 = 2+4rho)")
    print("   v:", [int(x) for x in v[:8]])
