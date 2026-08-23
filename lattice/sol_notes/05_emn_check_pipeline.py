import sys; sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/sol_notes')
import importlib.util
spec=importlib.util.spec_from_file_location("core","/home/ubuntu/code/math-modular-sources/lattice/sol_notes/05_emn_core.py")
core=importlib.util.module_from_spec(spec); spec.loader.exec_module(core)
from mpmath import mp, mpf, quad, log, catalan, pi as mpi, binomial as mpbin
mp.dps=60
CAT=+catalan

def num_I(m,t):
    """direct high-precision numeric of I(m,t) via the (v,w) reduction, done as a clean 1D w-integral."""
    def inner(w):
        A = 1+w**2
        # int_0^1 v^m/(2-Av)^{t+1} dv  = A^{-(m+1)} int_{2-A}^{2}(2-u)^m u^{-t-1} du
        lo = 2-A
        tot = mp.mpf(0)
        for j in range(m+1):
            c = mpbin(m,j)*mp.mpf(2)**(m-j)*(-1)**j
            if j==t:
                tot += c*log(2/lo)
            else:
                tot += c*(mp.mpf(2)**(j-t)-lo**(j-t))/(j-t)
        return A**(-(m+1))*tot*(1-w**2)**m
    val = quad(inner,[0,1])
    return mp.mpf(2)**(t+1)/2/mp.mpf(4)**m*val

def num_direct(m,t):
    """fully independent: original 2-D integral over Delta in (x,y)."""
    f = lambda x,y: x**m*y**m/(1-x**2-y**2)**(t+1)
    return quad(lambda x: quad(lambda y: f(x,y),[0,1-x]),[0,1])

def ev(expr):
    return mp.mpf(0)+ core.__dict__ and None

import sympy
def sy_num(expr):
    return mp.mpf(str(sympy.N(expr.subs(core.G, sympy.nsimplify(sympy.Symbol('C'))) if False else expr.subs(core.G, sympy.Catalan), 55)))

print("=== E0 ===")
e = core.I(0,0)
print("exact I(0,0) =", e)
print("numeric(reduction) =", num_I(0,0))
print("numeric(direct 2D) =", num_direct(0,0))
print("Catalan            =", CAT)

print()
print("=== cross-check exact vs numeric for various (m,t) ===")
for (m,t) in [(0,0),(2,0),(2,1),(2,2),(4,0),(4,1),(4,2),(4,3),(4,4),(6,3),(8,5),(4,-2),(6,-4)]:
    e = core.I(m,t)
    rat,gc,pic = core.split(e)
    ex = sy_num(e)
    if t>=0:
        nu = num_I(m,t)
        print(f"I({m},{t}): exact={rat} + ({gc})G   pi-part={pic}   |exact-num|={abs(ex-nu)}")
    else:
        print(f"I({m},{t}): exact={rat} + ({gc})G   pi-part={pic}   (t<0, no numeric)")
