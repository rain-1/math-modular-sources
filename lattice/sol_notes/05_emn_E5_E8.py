import sympy as sp
from sympy import Rational as R, symbols, solve, simplify, nsimplify, sqrt, diff, S, Poly, real_roots, Interval, maximum, minimum
x,y,g,c,h = symbols('x y g c h', real=True)

print("="*72); print("E5: max on Delta of F = x^4 y^4/(1-x^2-y^2)^2")
gg = 1-x**2-y**2
F  = x**4*y**4/gg**2
# interior critical points
cps = sp.solve([sp.diff(F,x), sp.diff(F,y)], [x,y], dict=True)
vals=[]
for s in cps:
    if s.get(x) is None or s.get(y) is None: continue
    xv,yv = sp.nsimplify(s[x]), sp.nsimplify(s[y])
    if xv.is_real and yv.is_real and xv>=0 and yv>=0 and xv+yv<1:
        vals.append(((xv,yv), sp.simplify(F.subs({x:xv,y:yv}))))
print("  interior critical points with 0<=x,y, x+y<1:", vals)
# boundary x+y=1: g = 1-x^2-y^2 = 2xy ; F = (xy)^4/(2xy)^2 = (xy)^2/4, xy in (0,1/4]
print("  on edge x+y=1: g=2xy, F=(xy)^2/4, max at xy=1/4 -> F=", R(1,4)**2/4, "= 1/64")
# edges x=0 or y=0: F=0
print("  on edges x=0,y=0: F=0")
print("  value at x=y=1/2:", sp.simplify(F.subs({x:R(1,2),y:R(1,2)})))
# global: parametrize by p=xy on each level. do it fully:
print("  => global max =", max([v for _,v in vals]+[R(1,64)]))

print()
print("="*72); print("E6: for fixed g, max over Delta of x^4y^4  =?= min(g,1-g)^4/16")
# On Delta with g fixed: x^2+y^2 = 1-g.  xy maximised at x=y => xy=(1-g)/2, provided x+y<=1
# x=y=sqrt((1-g)/2) needs 2 sqrt((1-g)/2)<=1 <=> 2(1-g)<=1 <=> g>=1/2
# if g<1/2 the constraint x+y<=1 binds: x+y=1 and x^2+y^2=1-g => 2xy = 1-(1-g)=g => xy=g/2
print("  regime g>=1/2 : interior, x=y=sqrt((1-g)/2), xy=(1-g)/2, (xy)^4=(1-g)^4/16")
print("  regime g<=1/2 : boundary x+y=1, 2xy=g, xy=g/2, (xy)^4=g^4/16")
print("  min(g,1-g)^4/16 matches both regimes:", "g<=1/2 -> min=g" , "g>=1/2 -> min=1-g")
# numeric confirm
import mpmath as mpm
mpm.mp.dps=30
def maxxy(gv):
    gv=mpm.mpf(gv); best=mpm.mpf(0); bx=None
    # param by x, y determined by x^2+y^2=1-gv
    r2=1-gv
    N=2000000
    # analytic candidates only + scan
    import math
    for i in range(1,20001):
        xv=mpm.mpf(i)/20000*mpm.sqrt(r2)
        s=r2-xv**2
        if s<0: continue
        yv=mpm.sqrt(s)
        if xv+yv<=1+mpm.mpf('1e-30') and xv>=0:
            if xv*yv>best: best=xv*yv; bx=(xv,yv)
    return best,bx
for gv in ['0.1','0.3','0.5','0.7','0.9']:
    b,bx = maxxy(gv)
    gvv=mpm.mpf(gv)
    pred = min(gvv,1-gvv)/2
    print(f"   g={gv}: scanned max xy={mpm.nstr(b,12)}  predicted min(g,1-g)/2={mpm.nstr(pred,12)}  at {(mpm.nstr(bx[0],8),mpm.nstr(bx[1],8))}")

print()
print("="*72); print("E7: W(g)=min(g,1-g)^4/(16 g^2);  256*max_g W(g)(1-4g^2)^2 =?= 16/27")
Wlo = g**4/(16*g**2)          # g<=1/2 : g^2/16
Whi = (1-g)**4/(16*g**2)      # g>=1/2
for name,Wx,dom in [("g<=1/2",Wlo,(S(0),R(1,2))),("g>=1/2",Whi,(R(1,2),S(1)))]:
    E = sp.simplify(256*Wx*(1-4*g**2)**2)
    d = sp.simplify(sp.diff(E,g))
    crit = [r for r in sp.solve(sp.numer(sp.together(d)),g) if r.is_real and dom[0]<=r<=dom[1]]
    print(f"  {name}: 256*W*(1-4g^2)^2 = {sp.factor(E)}")
    print(f"     critical pts in domain: {crit}")
    cand = crit+[dom[0],dom[1]]
    best=None
    for cc in cand:
        v=sp.simplify(E.subs(g,cc))
        if v.is_real and (best is None or v>best[1]): best=(cc,v)
    print(f"     max on this piece: g={best[0]}  value={sp.nsimplify(best[1])} = {sp.N(best[1],15)}")
print("  16/27 =", sp.N(R(16,27),15))
