import sympy as sp
from sympy import Rational as R, symbols, S
g,c = symbols('g c', positive=True)
import mpmath as mp
mp.mp.dps=40

# 256*W(g):  A-piece (0<g<=1/2): 16 g^2 ; B-piece (1/2<=g<1): 16 (1-g)^4/g^2
PA = 16*g**2
PB = 16*(1-g)**4/g**2

def sup_abs_exact(phi_expr, label, prec=25):
    """phi_expr : polynomial in g (already phi(g^2)). returns (sup, arg) numerically, exact crit pts."""
    out=[]
    for (Pp, lo, hi, nm) in [(PA, S(0), R(1,2),'A'), (PB, R(1,2), S(1),'B')]:
        F = sp.together(Pp*phi_expr)
        d = sp.numer(sp.together(sp.diff(F,g)))
        pol = sp.Poly(sp.expand(d), g)
        pol = sp.Poly(sp.factor_list(pol.as_expr())[0]*sp.prod([f for f,_ in sp.factor_list(pol.as_expr())[1]]), g) if pol.degree()>0 else pol
        rts = [sp.nsimplify(r) for r in sp.real_roots(pol)]
        cands=[lo,hi]+[r for r in rts if lo<=sp.N(r,30)<=hi]
        for cc in cands:
            v = abs(sp.N(F.subs(g,cc), prec))
            out.append((v, sp.N(cc,prec), nm))
    out.sort(reverse=True)
    return out

print("="*72); print("E7 recheck: phi(h)=(1-4h)^2 -> phi(g^2)=(1-4g^2)^2")
o=sup_abs_exact((1-4*g**2)**2,'E7')
for v,a,nm in o[:4]: print(f"   |F|={sp.N(v,18)} at g={sp.N(a,18)} piece {nm}")
print("   16/27 =", sp.N(R(16,27),18), "   1/(2*sqrt(3)) =", sp.N(1/(2*sp.sqrt(3)),18))

print(); print("="*72); print("E8a: phi(h)=1-c h.  Structure:")
print("   A-piece: F_A = 16 g^2 (1-c g^2); interior crit g*=1/sqrt(2c) (in (0,1/2) iff c>2), F_A(g*) = 4/c")
print("            at g=1/2: F_A = 4-c  (negative for c>4)")
print("   B-piece: F_B = 16 (1-g)^4 (1-c g^2)/g^2  (negative on [1/2,1) once c>4)")
def Mof(cv):
    cv=mp.mpf(cv)
    phi = 1 - sp.nsimplify(sp.Float(str(cv),40))*g**2
    o=sup_abs_exact(phi,'x',30)
    return o[0]
for cv in ['4.0','4.4','4.6','4.62','4.624','4.63','4.8','5.2']:
    o=Mof(cv); print(f"   c={cv:>6}: sup|.| = {sp.N(o[0],14)} at g={sp.N(o[1],12)} piece {o[2]};  4/c={sp.N(4/sp.Float(cv),14)}")
# minimax: 4/c  vs  max_{g in [1/2,1)} |F_B|
cs = symbols('cs', positive=True)
FB = 16*(1-g)**4*(cs*g**2-1)/g**2
dFB = sp.numer(sp.together(sp.diff(FB,g)))
print("   dFB numerator:", sp.factor(sp.expand(dFB)))
def maxFB(cv):
    cv=sp.nsimplify(sp.Float(str(cv),40))
    p=sp.Poly(sp.expand(dFB.subs(cs,cv)),g)
    p=sp.Poly(sp.prod([f for f,_ in sp.factor_list(p.as_expr())[1]]),g)
    rts=[r for r in sp.real_roots(p) if R(1,2)<=sp.N(r,30)<=1]
    vals=[sp.N(FB.subs({cs:cv,g:r}),30) for r in rts]+[sp.N(FB.subs({cs:cv,g:R(1,2)}),30)]
    return max(vals), rts
lo,hi=mp.mpf(4),mp.mpf(6)
for _ in range(200):
    mid=(lo+hi)/2
    v,_=maxFB(mid)
    if float(v) < float(4/mid): lo=mid
    else: hi=mid
    if hi-lo<mp.mpf('1e-30'): break
copt=(lo+hi)/2
v,rts=maxFB(copt)
print("   equioscillation c* =", mp.nstr(copt,20))
print("   M(c*) = 4/c* =", mp.nstr(4/copt,20), " ; max|F_B| =", sp.N(v,20), " at g=",[sp.N(r,15) for r in rts])
# exact: solve  max_g F_B(g,c) = 4/c  -> resultant
print()
print("   exact characterisation: c* is a root of Res_g( d/dg F_B , F_B - 4/c )")
res = sp.resultant(sp.expand(dFB), sp.expand(sp.numer(sp.together(FB - 4/cs))), g)
res = sp.factor(sp.expand(res))
print("   resultant factors:", res)
for f,_ in sp.factor_list(res)[1]:
    p=sp.Poly(f,cs)
    if p.degree()>0:
        rr=[sp.N(r,20) for r in sp.real_roots(p)]
        print("     factor", f, " real roots:", rr)

print(); print("="*72); print("E8b: q(h)=1-10h+23h^2 ; 256 max_g W(g)|q(g^2)|")
o=sup_abs_exact(1-10*g**2+23*g**4,'q',30)
for v,a,nm in o[:5]: print(f"   |F|={sp.N(v,18)} at g={sp.N(a,18)} piece {nm}")

print(); print("="*72); print("optimal QUADRATIC in h (minimax) -- compare with 16/27 and with q")
def supq(u,v):
    phi = 1 + sp.nsimplify(sp.Float(u,25))*g**2 + sp.nsimplify(sp.Float(v,25))*g**4
    return float(sup_abs_exact(phi,'',25)[0][0])
cur=[-10.0,23.0]; curv=supq(*cur); step=1.0
for it in range(400):
    improved=False
    for du,dv in [(1,0),(-1,0),(0,1),(0,-1),(1,1),(-1,-1),(1,-1),(-1,1)]:
        nu,nv=cur[0]+du*step, cur[1]+dv*step
        val=supq(nu,nv)
        if val<curv-1e-16: cur=[nu,nv]; curv=val; improved=True; break
    if not improved:
        step/=2
        if step<1e-8: break
print("   minimax quadratic phi(h)=1+u h+v h^2 : u=%.10f v=%.10f  value=%.10f"%(cur[0],cur[1],curv))
print("   q(h)=1-10h+23h^2 gives %.10f ; (1-4h)^2 = 1-8h+16h^2 gives %.10f"%(supq(-10,23), supq(-8,16)))
print("   best LINEAR-in-h  (1-c h): %.10f at c*=4.6240231347"%(4/4.6240231346590769))
