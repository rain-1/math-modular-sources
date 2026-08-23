import importlib.util, sympy
from sympy import binomial, Rational, simplify, nsimplify, Integer
spec=importlib.util.spec_from_file_location("core","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_core.py")
core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
from mpmath import mp, quad
mp.dps=50

# independent numeric check of a negative-t case
f = lambda x,y: x**4*y**4*(1-x**2-y**2)
print("check I(4,-2): exact", core.I(4,-2), " numeric", quad(lambda x: quad(lambda y: f(x,y),[0,1-x]),[0,1]))
print()
print("="*70); print("E2: b(m,t) table.  claim: b = 4^-m C(m,m/2) C(m,t)")
print(f"{'m':>3}{'t':>4}   {'computed b(m,t)':>22}   {'claimed':>22}   ratio")
bad=[]
for m in [0,2,4,6,8]:
    for t in range(0,m+1):
        rat,gc,pic = core.split(core.I(m,t))
        assert pic==0
        claim = Rational(1,4**m)*binomial(m,m//2)*binomial(m,t)
        r = simplify(gc/claim) if claim!=0 else None
        print(f"{m:>3}{t:>4}   {str(gc):>22}   {str(claim):>22}   {r}")
        if gc!=claim: bad.append((m,t,gc,claim))
print("mismatches:", len(bad), " -- all ratios are (-1)^t:",
      all(simplify(g/c)==Integer(-1)**t for (m,t,g,c) in bad))
print()
print("="*70); print("E3 recurrence residuals")
print("  I(m,t) =?= -1/(2^{t+1} t (t-1) C(2(m-t), m-t)) + (m-1)^2/(4t(t-1)) I(m-2,t-2)")
for (m,t) in [(4,2),(6,2),(6,4),(8,2),(8,4),(8,6),(10,2),(10,6),(10,10),(12,8)]:
    lhs = core.I(m,t)
    rhs = -Rational(1,1)/(Integer(2)**(t+1)*t*(t-1)*binomial(2*(m-t), m-t)) + Rational((m-1)**2,4*t*(t-1))*core.I(m-2,t-2)
    res = sympy.expand(lhs-rhs)
    print(f"  (m,t)=({m},{t}) residual = {sympy.simplify(res)}")
print()
print("="*70); print("E4  b(4n,2n) = C(4n,2n)^2/4^{4n} ~ 1/(2 pi n)")
for n in [1,2,3,4]:
    rat,gc,pic = core.split(core.I(4*n,2*n))
    cl = Rational(1,4**(4*n))*binomial(4*n,2*n)**2
    print(f"  n={n}: computed b(4n,2n)={gc}, C(4n,2n)^2/4^(4n)={cl}, equal={gc==cl}")
mp.dps=30
from mpmath import binomial as B, mpf, pi as MPI
print("  asymptotic ratio b*(2 pi n):")
for n in [1,2,5,10,20,50,100,200]:
    b = B(4*n,2*n)**2/mpf(4)**(4*n)
    print(f"    n={n:>4}  b={mp.nstr(b,10)}  b*2*pi*n={mp.nstr(b*2*MPI*n,12)}")
