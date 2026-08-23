#!/usr/bin/env python3
"""Kodaira-admissible normalisation classes for the four-term (five singular point) rows.

Class = (rho ; M, j1, j2):
    P(n) = a n^2 + a(1-rho) n + c
    Q(n) = d n^2 - 2 rho d n + f
    R(n) = C (Mn - j1)(Mn - j2),      j1 + j2 = (1+3rho) M,   s_i = j_i/M
exponents:  (0,0) at t=0 ;  (0,rho) at the three roots of 1-at+dt^2-gt^3 ;
            (2-s_1, 2-s_2) at infinity ;  delta_inf = s_2-s_1.
"""
from fractions import Fraction as F

ADM = {F(0), F(1,2), F(1,3), F(2,3)}
def admissible(x):
    return F(x) - int(F(x)) in ADM or (F(x) % 1) in ADM

def kod(diff):
    diff = F(diff) % 1
    return 'cusp' if diff == 0 else ('ord2' if diff == F(1,2) else 'ord3')

def emin(kind): return {'cusp':1,'ord2':3,'ord3':2}[kind]

def classes(rhos, deltas, Mmax=12):
    out = []
    for rho in rhos:
        rho = F(rho)
        if (rho % 1) not in ADM: continue
        for de in deltas:
            de = F(de)
            if (de % 1) not in ADM: continue
            s1 = (1+3*rho-de)/2; s2 = (1+3*rho+de)/2
            M = (s1.denominator*s2.denominator)//__import__('math').gcd(s1.denominator,s2.denominator)
            if M > Mmax: continue
            j1, j2 = int(s1*M), int(s2*M)
            # euler budget for a RATIONAL elliptic surface
            kinds = ['cusp'] + [kod(rho)]*3 + [kod(de)]
            emn = sum(emin(k) for k in kinds)
            out.append(dict(rho=rho, delta=de, M=M, j1=j1, j2=j2,
                            kinds=kinds, emin=emn, ref=(emn<=12)))
    return out

if __name__ == '__main__':
    rhos = [F(0),F(1,3),F(1,2),F(2,3),F(1),F(4,3),F(3,2),F(2),F(-1),F(-1,2),F(-1,3),F(-2,3)]
    deltas = [F(0),F(1,3),F(1,2),F(2,3),F(1),F(4,3),F(3,2),F(5,3),F(2),F(7,3),F(5,2),F(3),F(4)]
    cs = classes(rhos, deltas)
    print(f"{'rho':>6} {'delta':>6} {'M':>3} {'j1':>4} {'j2':>4}  {'R(n)':<26} {'types (0,t_i x3,inf)':<34} emin REF")
    for c in cs:
        M,j1,j2 = c['M'],c['j1'],c['j2']
        Rs = f"C({M}n-{j1})({M}n-{j2})".replace('n-0','n').replace('--','+')
        print(f"{str(c['rho']):>6} {str(c['delta']):>6} {M:>3} {j1:>4} {j2:>4}  {Rs:<26} "
              f"{c['kinds'][0]},{c['kinds'][1]}x3,{c['kinds'][4]:<8} {c['emin']:>4} {'yes' if c['ref'] else 'NO'}")
    print(len(cs), "classes")
