"""
05_emn_core.py  -- exact machinery for
    I(m,t) = int_Delta x^m y^m / (1-x^2-y^2)^{t+1} dx dy,   Delta = {x,y>=0, x+y<=1}

DERIVATION (all exact, no numerics):
  s=x+y, d=x-y, dxdy=(1/2)ds dd, 1-x^2-y^2 = (2-s^2-d^2)/2, xy=(s^2-d^2)/4.
  d=s w  =>
     I = 2^{t+1} int_0^1 int_0^1 s * (s^2(1-w^2)/4)^m / (2-s^2(1+w^2))^{t+1} dw ds
  v=s^2, A=1+w^2:
     inner_v = int_0^1 v^m/(2-Av)^{t+1} dv = A^{-(m+1)} int_{2-A}^{2} (2-u)^m u^{-t-1} du
             = A^{-(m+1)} sum_j C(m,j) 2^{m-j} (-1)^j * int_{1-w^2}^{2} u^{j-t-1} du
     the j=t term (only when 0<=t<=m) gives log(2/(1-w^2)); all others (2^{j-t}-(1-w^2)^{j-t})/(j-t).
  Hence
     I(m,t) = 2^{t+1}/(2*4^m) * [ sum_{j != t} C(m,j)2^{m-j}(-1)^j/(j-t) * (2^{j-t} X(0) - X(j-t))
                                  + [t in 0..m] * C(m,t) 2^{m-t} (-1)^t * J_m ]
  with   X(p) = int_0^1 (1-w^2)^{m+p}/(1+w^2)^{m+1} dw          (rational + rational*pi)
         J_m  = int_0^1 (1-w^2)^m log(2/(1-w^2))/(1+w^2)^{m+1} dw
  w=tan(theta), phi=2theta:  (1-w^2)^m/(1+w^2)^{m+1} dw = cos^m(phi) dphi/2, 1-w^2 = 2cos phi/(1+cos phi)
     J_m = (1/2) K_m,  K_m = int_0^{pi/2} cos^m phi * log((1+cos phi)/cos phi) dphi
  log((1+c)/c) = int_0^1 dlam/(c+lam):
     K_m = sum_{i=0}^{m-1} (-1)^i W(m-1-i)/(i+1) + (-1)^m S_m,  S_m = int_0^1 lam^m R(lam) dlam
     R(lam) = int_0^{pi/2} dphi/(cos phi + lam) = (2/sin psi) artanh(tan(psi/2)),  lam = cos psi
     => S_m = -2 int_0^{pi/4} sin^m(2u) log tan u du
  For m=2r even, sin^{2r}(2u) = 4^{-r}[C(2r,r) + 2 sum_{k=1}^r (-1)^k C(2r,r-k) cos(4ku)], and with
     A_k := int_0^{pi/4} cos(4ku) log tan u du,  log tan u = -2 sum_{j odd} cos(2ju)/j:
     A_0 = -G ;  A_k = -(1/(2k)) sum_{i=0}^{k-1} (-1)^i/(2i+1)   (RATIONAL for k>=1)
  So ALL of Catalan's constant enters through the single term A_0.
"""
from sympy import Rational, binomial, pi, symbols, nsimplify, Integer, expand, simplify, factorint, lcm, S as Sym

G = symbols('G')          # Catalan
PI = pi

# ---------- Z(k) = int_0^1 dw/(1+w^2)^k  (k any integer) ----------
_Zc = {}
def Z(k):
    if k in _Zc: return _Zc[k]
    if k <= 0:
        n = -k                       # int_0^1 (1+w^2)^n dw
        v = sum(binomial(n, j)*Rational(1, 2*j+1) for j in range(n+1))
    elif k == 1:
        v = PI/4
    else:
        j = k-1                      # Z(j+1) = 1/(2j*2^j) + (2j-1)/(2j) Z(j)
        v = Rational(1, 2*j*2**j) + Rational(2*j-1, 2*j)*Z(j)
    v = expand(v); _Zc[k] = v; return v

# ---------- X(m,p) = int_0^1 (1-w^2)^{m+p}/(1+w^2)^{m+1} dw ----------
_Xc = {}
def X(m, p):
    key = (m, p)
    if key in _Xc: return _Xc[key]
    q = m + p
    assert q >= 0, (m, p)
    # (1-w^2)^q = sum_i C(q,i) 2^{q-i} (-1)^i (1+w^2)^i
    v = sum(binomial(q, i)*Integer(2)**(q-i)*(-1)**i*Z(m+1-i) for i in range(q+1))
    v = expand(v); _Xc[key] = v; return v

# ---------- W(k) = int_0^{pi/2} cos^k ----------
_Wc = {0: PI/2, 1: Integer(1)}
def W(k):
    if k in _Wc: return _Wc[k]
    v = expand(Rational(k-1, k)*W(k-2)); _Wc[k] = v; return v

# ---------- A_k ----------
def A(k):
    if k == 0: return -G
    return -Rational(1, 2*k)*sum(Rational((-1)**i, 2*i+1) for i in range(k))

# ---------- S_m (m even) ----------
_Sc = {}
def Sm(m):
    assert m % 2 == 0
    if m in _Sc: return _Sc[m]
    r = m//2
    v = -Rational(2,1)/Integer(4)**r * (binomial(2*r, r)*A(0)
        + 2*sum((-1)**k*binomial(2*r, r-k)*A(k) for k in range(1, r+1)))
    v = expand(v); _Sc[m] = v; return v

# ---------- J_m = K_m/2 ----------
_Jc = {}
def J(m):
    assert m % 2 == 0
    if m in _Jc: return _Jc[m]
    K = sum(Rational((-1)**i, i+1)*W(m-1-i) for i in range(m)) + Sm(m)
    v = expand(K/2); _Jc[m] = v; return v

# ---------- I(m,t) ----------
_Ic = {}
def I(m, t):
    """exact value of the integral, as a sympy expression in G (and pi, which should cancel)."""
    key = (m, t)
    if key in _Ic: return _Ic[key]
    tot = Sym.Zero
    for j in range(m+1):
        if j == t: continue
        c = binomial(m, j)*Integer(2)**(m-j)*Integer(-1)**j*Rational(1, j-t)
        tot += c*(Integer(2)**(j-t)*X(m, 0) - X(m, j-t))
    if 0 <= t <= m:
        tot += binomial(m, t)*Integer(2)**(m-t)*Integer(-1)**t*J(m)
    v = expand(Integer(2)**(t+1)*Rational(1, 2)*Rational(1, 4**m)*tot)
    v = expand(v); _Ic[key] = v; return v

def split(expr):
    """return (rational_part, G_coeff, pi_part_coeff) ; pi_part should be 0."""
    e = expand(expr)
    gc = e.coeff(G, 1)
    rest = expand(e - gc*G)
    pic = expand(rest).coeff(PI, 1)
    rat = expand(rest - pic*PI)
    return simplify(rat), simplify(gc), simplify(pic)
