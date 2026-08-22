#!/usr/bin/env python3
"""Local exponent data of a second-order Apery-like row.

Row:   (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1},   P = a n^2 + b n + c,
                                                    Q = d n^2 + e n + f.
Operator (theta = t d/dt):   L = theta^2 - t P(theta) + t^2 Q(theta+1).
In D = d/dt form:
   L = t^2 R(t) D^2 + t S(t) D + (-c t + (d+e+f) t^2),
   R(t) = 1 - a t + d t^2,   S(t) = 1 - (a+b) t + (3d+e) t^2.
Singular points: 0, roots t_1,t_2 of R, and infinity.
 * exponents at 0        : (0,0)                       [MUM]
 * exponents at t_i      : (0, rho_i),  rho_i = -T(t_i)/(t_i R'(t_i)),
                           T(t) = S(t) - t R'(t) = 1 - b t + (d+e) t^2
 * exponents at infinity : (1-r_1, 1-r_2), r_i the roots of Q
Fuchs relation: rho_1 + rho_2 = r_1 + r_2 = -e/d.
rho_1, rho_2 are the roots of a rational quadratic:
   rho_1 + rho_2 = -e/d,  rho_1 rho_2 = -Res(R,T)/(d (a^2-4d)).
"""
import sys
from sympy import Rational, sqrt, symbols, Poly, resultant, simplify, nsimplify, Integer, S

t = symbols('t')

def data(a,b,c,d,e,f):
    a,b,c,d,e,f = map(S, (a,b,c,d,e,f))
    out = {}
    R = 1 - a*t + d*t**2
    T = 1 - b*t + (d+e)*t**2
    out['R'] = R
    if d == 0:
        out['note'] = 'd=0: only 3 singular points (t=0, 1/a, infinity)'
        return out
    disc = a**2 - 4*d
    out['disc'] = disc
    # rho_1, rho_2
    s_rho = -e/d
    if disc == 0:
        out['note'] = 'R has a double root: only 3 singular points'
        return out
    Tp = Poly(T, t)
    prodT = resultant(Poly(R,t), Tp)/d**Tp.degree()   # = T(t_1)T(t_2)
    p_rho = prodT/((S(1)/d)*(-disc))                  # rho_1 rho_2
    out['rho_sum'], out['rho_prod'] = simplify(s_rho), simplify(p_rho)
    dd = simplify(s_rho**2 - 4*p_rho)
    out['rho'] = (simplify((s_rho - sqrt(dd))/2), simplify((s_rho + sqrt(dd))/2))
    # r_1, r_2 (roots of Q)
    ddq = simplify(e**2 - 4*d*f)
    out['r'] = (simplify((-e - sqrt(ddq))/(2*d)), simplify((-e + sqrt(ddq))/(2*d)))
    out['delta_inf'] = simplify(sqrt(ddq)/d) if d else None
    out['exp_inf'] = (simplify(1-out['r'][0]), simplify(1-out['r'][1]))
    # characteristic roots lambda_i = 1/t_i : roots of x^2 - a x + d
    out['lam'] = (simplify((a - sqrt(disc))/2), simplify((a + sqrt(disc))/2))
    out['t_sing'] = (simplify((a - sqrt(disc))/(2*d)), simplify((a + sqrt(disc))/(2*d)))
    return out

def frac_mod1(x):
    """x mod 1 as a Rational if x is rational, else None"""
    try:
        r = Rational(x)
    except Exception:
        return None
    return r - r.__floor__()

ORDER = {Rational(0):'cusp (I_n / I_n^*)', Rational(1,2):'order 2 (III / III^*)',
         Rational(1,3):'order 3 (II,IV,IV^*,II^*)', Rational(2,3):'order 3 (II,IV,IV^*,II^*)'}

def signature(a,b,c,d,e,f):
    D = data(a,b,c,d,e,f)
    if 'rho' not in D: return D, None
    kinds = []
    ok = True
    for x in D['rho']:
        m = frac_mod1(x)
        if m is None or m not in ORDER: ok = False; m = m
        kinds.append(m)
    m = frac_mod1(abs(D['delta_inf']))
    if m is not None and m > Rational(1,2): m = 1 - m
    kinds.append(m)
    if m is None or m not in ORDER: ok = False
    D['kinds'] = kinds          # [t_1, t_2, infinity]; 0 = cusp, 1/2, 1/3 elliptic order
    D['PSL2Z_admissible'] = ok
    if ok:
        ell = [k for k in kinds if k != 0]
        mu = 6*(2 - sum(Rational(1, 2 if k==Rational(1,2) else 3) for k in ell))
        D['mu'] = mu
        D['cusps'] = 1 + sum(1 for k in kinds if k == 0)   # +1 for the MUM point
        D['elliptic'] = [2 if k==Rational(1,2) else 3 for k in ell]
    return D, ok

if __name__ == '__main__':
    rows = [
      # name, P coeffs (a,b,c), Q coeffs (d,e,f)
      ("Zagier A (7,2,-8)",   (7,7,2),   (-8,0,0)),
      ("Zagier B (9,3,27)",   (9,9,3),   (27,0,0)),
      ("Zagier C (10,3,9)",   (10,10,3), (9,0,0)),
      ("Zagier D (11,5,125)", (11,11,5), (125,0,0)),
      ("Zagier E (12,4,32)",  (12,12,4), (32,0,0)),
      ("Zagier F (17,6,72)",  (17,17,6), (72,0,0)),
      ("hyp (16,4,0)",        (16,16,4), (0,0,0)),
      ("root Apery (136,10,4)",  (136,68,10),  (16,-16,4)),
      ("root T (24,2,4)",        (24,12,2),    (16,-16,4)),
      ("root Domb (20,2,16)",    (20,10,2),    (64,-64,16)),
      ("root AZ(9,3,-27) (72,6,-108)", (72,36,6), (-432,432,-108)),
      ("root AZ(11,5,125)",      (88,44,10),   (2000,-2000,500)),
      ("root AZ(7,3,81)",        (56,28,6),    (1296,-1296,324)),
      ("root Cooper s7",         (26,13,2),    (-27,27,-6)),
      ("root Cooper s10",        (24,12,2),    (-256,256,-60)),
      ("root Cooper s18",        (56,28,6),    (768,-768,180)),
      ("NCS rho=7/6 row",        (234,-39,-78),(-27,63,-30)),
    ]
    for name,(a,b,c),(d,e,f) in rows:
        D,ok = signature(a,b,c,d,e,f)
        print("="*70); print(name, " P=%s n^2+%s n+%s  Q=%s n^2+%s n+%s"%(a,b,c,d,e,f))
        if 'rho' not in D: print("  ", D.get('note')); continue
        print("  rho (finite sing.) =", D['rho'], "  r (roots of Q) =", D['r'])
        print("  exponents at oo =", D['exp_inf'], "  delta_oo =", D['delta_inf'])
        print("  lambda_1,2 =", [complex(x) for x in D["lam"]], " disc =", D['disc'])
        print("  PSL2(Z)-admissible:", D.get('PSL2Z_admissible'), " kinds(t1,t2,oo) =", D.get('kinds'))
        if D.get('PSL2Z_admissible'):
            print("  signature: %d cusps, elliptic orders %s,  mu = %s"%(D['cusps'],D['elliptic'],D['mu']))
