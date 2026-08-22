#!/usr/bin/env python3
"""Independent exact verification of the two new Herfurtner rows."""
from fractions import Fraction as F
from math import gcd
def run(name,P,Q,N=300):
    u=[F(1),F(P(0))]
    for n in range(1,N):
        u.append((P(n)*u[n]-Q(n)*u[n-1])/F((n+1)**2))
    bad=[n for n,x in enumerate(u) if x.denominator!=1]
    b=[F(0),F(1)]
    for n in range(1,80):
        b.append((P(n)*b[n]-Q(n)*b[n-1])/F((n+1)**2))
    D=[1]*82
    for n in range(1,82): D[n]=D[n-1]*n//gcd(D[n-1],n)
    k=0
    while k<=6:
        if all((D[n]**k*b[n]).denominator==1 for n in range(81)): break
        k+=1
    ksharp = not all((D[n]**(k-1)*b[n]).denominator==1 for n in range(81)) if k>0 else True
    print("%s: integral to n=%d: %s ; first 8 = %s ; k=%d (sharp: %s)"%(
        name,N,("YES" if not bad else "NO at n=%d"%bad[0]),[int(x) for x in u[:8]],k,ksharp))
run("row A (117,21,441) on I1 I7 II II",
    lambda n:117*n*n+78*n+21, lambda n:441*(3*n-1)**2)
run("row B (72,6,108) on I3 III III III",
    lambda n:72*n*n+36*n+6, lambda n:108*(4*n-1)*(4*n-3))
# control: they are not rescalings of an integral row with smaller A
for nm,d,P,Q in [("A/3",3,lambda n:39*n*n+26*n+7, lambda n:49*(3*n-1)**2),
                 ("B/2",2,lambda n:36*n*n+18*n+3, lambda n:27*(4*n-1)*(4*n-3)),
                 ("B/3",3,lambda n:24*n*n+12*n+2, lambda n:12*(4*n-1)*(4*n-3)),
                 ("B/6",6,lambda n:12*n*n+6*n+1, lambda n:3*(4*n-1)*(4*n-3))]:
    u=[F(1),F(P(0))]
    ok=True
    for n in range(1,40):
        u.append((P(n)*u[n]-Q(n)*u[n-1])/F((n+1)**2))
        if u[-1].denominator!=1: ok=False; break
    print("  descaled %s integral? %s"%(nm,ok))
