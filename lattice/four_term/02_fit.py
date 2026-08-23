#!/usr/bin/env python3
"""Fit  (n+1)^2 u_{n+1} = P(n)u_n - Q(n)u_{n-1} + R(n)u_{n-2}  to a sequence,
and construction experiments producing four-term rows from three-term ones."""
from fractions import Fraction as F
import sympy as sp

def fit(u, N=None, degP=2, degQ=2, degR=2):
    """u: dict/list with u[-2]=u[-1]=0, u[0..N]. Returns (P,Q,R) coeff lists or None."""
    if N is None: N = len(u)-1
    U = lambda k: (0 if k < 0 else u[k])
    rows = []
    for k in range(0, N):
        r = [F(k)**i for i in range(degP+1)] * 0
        r = []
        for i in range(degP+1): r.append(F(k)**i * U(k))
        for i in range(degQ+1): r.append(-F(k)**i * U(k-1))
        for i in range(degR+1): r.append(F(k)**i * U(k-2))
        r.append(-F(k+1)**2 * U(k+1))
        rows.append(r)
    M = sp.Matrix([[sp.Rational(x) for x in r] for r in rows])
    ns = M.nullspace()
    out = []
    for v in ns:
        v = list(v)
        if v[-1] == 0: continue
        v = [x/v[-1] for x in v]
        out.append(v[:-1])
    return out

def order2(P, Q, N):
    u = {-2:0,-1:0,0:F(1)}
    for k in range(N):
        u[k+1] = (P(k)*u[k] - Q(k)*u[k-1])/F((k+1)**2)
    return u

def order3(P,Q,R,N):
    u = {-2:0,-1:0,0:F(1)}
    for k in range(N):
        u[k+1] = (P(k)*u[k] - Q(k)*u[k-1] + R(k)*u[k-2])/F((k+1)**2)
    return u

def binom(n,k):
    from math import comb
    return comb(n,k) if 0<=k<=n else 0

if __name__ == '__main__':
    # Zagier A : P = 7n^2+7n+2, Q = -8n^2
    P = lambda k: 7*k*k+7*k+2 ; Q = lambda k: -8*k*k
    N = 26
    u = order2(P,Q,N)
    print("Zagier A u_n:", [u[k] for k in range(8)])
    print("fit as 4-term (should be R=0 family):")
    for v in fit(u,N): print("   ", v)

    # (1) gauge  y -> y/(1-nu t):  v_n = sum_{k<=n} nu^{n-k} u_k
    for nu in (1,2,3,-1):
        v = {-2:0,-1:0}
        s = F(0)
        for k in range(N+1):
            s = s*nu + u[k]; v[k] = s
        sols = fit(v,N)
        print(f"\n[y/(1-{nu}t)] v =", [v[k] for k in range(6)], "-> #sols", len(sols))
        for x in sols: print("     ", x)

    # (2) signed binomial transform with general alpha (integer alpha=1): 
    for lam in (1,2,-1):
        v = {-2:0,-1:0}
        for k in range(N+1):
            v[k] = sum(F(binom(k,m))*F(-lam)**(k-m)*u[m] for m in range(k+1))
        sols = fit(v,N)
        print(f"\n[binom lam={lam} alpha=1] v =", [v[k] for k in range(6)], "-> #sols", len(sols))
        for x in sols: print("     ", x)
