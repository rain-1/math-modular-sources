"""Fourier-Laplace transform of L_f = (1-t^2)d^4 - (6t+2)d^3 - 6d^2 - 1
   under t -> -D_z, d_t -> z  (a Weyl-algebra automorphism).
   Claim: FL(L_f) = z^4 D^2 + 2 z^3 D + (-z^4 + 2z^3 + 1),
   whose normal form (W = zY) is W'' = (1 - 2/z - 1/z^4) W,
   i.e. the z -> 1/z image of W'' = (1/z^4 - 2/z^3 - 1) W.  Both irregular at 0 and oo."""
import sympy as sp
z, t = sp.symbols('z t')
D = sp.Symbol('D', commutative=False)
Z = sp.Symbol('Z', commutative=False)
# work with an explicit action on test functions instead of noncommutative algebra
u = sp.Function('u')
def apply_FL(expr_of_u):
    return sp.expand(expr_of_u)
# L_f with all t's on the LEFT:  d^4 - t^2 d^4 - 6 t d^3 - 2 d^3 - 6 d^2 - 1
# FL:  t -> -D,  d -> z   (order preserved: left factors stay left)
#      d^4 -> z^4 ; t^2 d^4 -> (-D)^2 z^4 = D^2 z^4 ; t d^3 -> -D z^3 ; etc.
f = sp.Function('f')
E = ( z**4*f(z)
      - sp.diff(z**4*f(z), z, 2)
      + 6*sp.diff(z**3*f(z), z)
      - 2*z**3*f(z) - 6*z**2*f(z) - f(z) )
E = sp.expand(sp.simplify(E))
print("FL(L_f) applied to f:", sp.collect(E, f(z)))
print("times -1          :", sp.expand(-E))
# check the claimed form
claim = z**4*sp.diff(f(z),z,2) + 2*z**3*sp.diff(f(z),z) + (-z**4+2*z**3+1)*f(z)
print("matches z^4 f'' + 2z^3 f' + (-z^4+2z^3+1) f ?", sp.simplify(sp.expand(-E - claim)) == 0)
# normal form W = z Y
W = sp.Function('W')
nf = sp.expand(sp.simplify((z**4*sp.diff(W(z)/z,z,2) + 2*z**3*sp.diff(W(z)/z,z) + (-z**4+2*z**3+1)*(W(z)/z))*z/z**3))
print("normal form (W=zY), divided by z^3 :", sp.simplify(nf), "  ->  W'' = (1 - 2/z - 1/z^4) W")
