#!/usr/bin/env python3
"""Theorem R2 (descent).  L_1 = th^2 + p*th + r,  L_w = Sym^w L_1 (order w+1, monic).
Compute L_w(g^w)/g^w in terms of u = th g / g and psi = (L_1 g)/g, and check:
  * the psi-free part vanishes identically (that IS the definition of Sym^w);
  * the coefficient of th^{w-1} psi is the constant w;
  * every other linear-in-psi coefficient is a weighted-homogeneous polynomial in
    u,p,r and their th-derivatives of POSITIVE weight (hence in x Q[[x]] when p,r in xQ[[x]]).
Weights: wt(u)=wt(p)=1, wt(r)=wt(psi)=2, wt(th)=+1.
"""
import sympy as sp

MAXW = 5
x = sp.Symbol('x')

def build(w, with_psi, with_p=True):
    """return (Qlist, symbols) where Q_j = th^j(g^w)/g^w."""
    # symbols: u^{(i)}, p^{(i)}, r^{(i)}, psi^{(i)}
    NS = w+3
    U = sp.symbols('u0:%d'%NS); P = sp.symbols('p0:%d'%NS)
    R = sp.symbols('r0:%d'%NS); PS = sp.symbols('s0:%d'%NS)
    def D(e):
        out = 0
        for arr in (U,P,R,PS):
            for i in range(NS-1):
                out += sp.diff(e, arr[i])*arr[i+1]
        # th u = psi - u^2 - p u - r
        out = out.subs(U[1], PS[0]-U[0]**2-(P[0]*U[0] if with_p else 0)-R[0], simultaneous=True)
        return sp.expand(out)
    # need th of u1 etc: define recursively by differentiating the Riccati relation
    # easier: substitute u_{i+1} -> th(u_i) computed from the relation, top-down
    # Build explicit expressions for th^i u in terms of u0,p*,r*,psi*
    thu = [U[0]]
    cur = PS[0]-U[0]**2-(P[0]*U[0] if with_p else 0)-R[0]
    def Dfull(e):
        out = 0
        for arr,base in ((P,None),(R,None),(PS,None)):
            for i in range(NS-1):
                out += sp.diff(e, arr[i])*arr[i+1]
        out += sp.diff(e, U[0])*(PS[0]-U[0]**2-(P[0]*U[0] if with_p else 0)-R[0])
        return sp.expand(out)
    Q=[sp.Integer(1)]
    for j in range(1, w+3):
        Q.append(sp.expand(Dfull(Q[-1]) + w*U[0]*Q[-1]))
    if not with_psi:
        Q=[sp.expand(q.subs({s:0 for s in PS})) for q in Q]
    return Q, (U,P,R,PS)

def wt(mon, U,P,R,PS):
    """weight of a monomial expression term"""
    tot=0
    for i,s in enumerate(U): tot += sp.degree(mon, s)*(1+i)
    for i,s in enumerate(P): tot += sp.degree(mon, s)*(1+i)
    for i,s in enumerate(R): tot += sp.degree(mon, s)*(2+i)
    for i,s in enumerate(PS): tot += sp.degree(mon, s)*(2+i)
    return tot

for w in range(2, MAXW+1):
    Q,(U,P,R,PS) = build(w, True)
    Q0 = [sp.expand(q.subs({s:0 for s in PS})) for q in Q]
    # solve for c_j, j=0..w, with N = th^{w+1} + sum_{j<=w} c_j th^j, from sum c_j Q0_j = 0
    cs = sp.symbols('c0:%d'%(w+1))
    expr = sp.expand(Q0[w+1] + sum(cs[j]*Q0[j] for j in range(w+1)))
    poly = sp.Poly(expr, U[0])
    eqs = poly.all_coeffs()
    sol = sp.solve(eqs, cs, dict=True)
    assert len(sol)==1, (w,sol)
    sol=sol[0]
    print("="*70)
    print("w =", w, "   Sym^w(th^2 + p th + r) = th^%d"%(w+1), end="")
    for j in range(w,-1,-1):
        print("  + (%s) th^%d"%(sp.factor(sp.simplify(sol[cs[j]])), j), end="")
    print()
    # now the full expression
    full = sp.expand(Q[w+1] + sum(sol[cs[j]]*Q[j] for j in range(w+1)))
    full = sp.expand(sp.simplify(full))
    # psi-free part must vanish
    zero = sp.expand(full.subs({s:0 for s in PS}))
    print("   psi-free part == 0 :", sp.simplify(zero)==0)
    # split into linear-in-psi and higher
    pol = sp.Poly(full, *PS)
    lin = {}; higher = 0
    for mono, coeff in pol.terms():
        d = sum(mono)
        if d==1:
            i = mono.index(1); lin[i]=lin.get(i,0)+coeff
        elif d>=2:
            higher += coeff*sp.prod([PS[i]**mono[i] for i in range(len(PS))])
    for i in sorted(lin):
        c = sp.expand(lin[i])
        terms = c.as_ordered_terms()
        wts = set()
        for t in terms:
            wts.add(wt(t,U,P,R,PS))
        print("   coeff of th^%d psi : %s      [weights %s]"%(i, sp.factor(c), sorted(wts)))
    print("   quadratic+ part in psi:", sp.factor(sp.expand(higher)) if higher!=0 else 0)
