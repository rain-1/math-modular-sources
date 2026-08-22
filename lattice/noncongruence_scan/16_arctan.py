#!/usr/bin/env python3
"""The arctan Pade family at rho = 7/6:
     (n+1)^2 a_{n+1} = (3n-2)[ (al/6)(2n+1) a_n + 3(3n-5) a_{n-1} ],  a_0=1, a_1=-al/3
   i.e. (al,be,ga,de,ep,ze) = (al, -al/6, -al/3, -27, 63, -30).
   Integral exactly for al = 18 mod 36;  lam_{1,2} = (al +- sqrt(al^2+108))/2;
   k = 2 sharp;  and the Apery limit is
        xi = -(1/(2 sqrt 3)) arctan( 6 sqrt 3 / al ).
   (This is the delta = -27 analogue of the Legendre/Pade log family at delta = 16,
    for which xi = -(1/12) log((x+1)/(x-1)) with al = 8x.)                        """
from fractions import Fraction as Fr
from math import gcd, log, sqrt
import mpmath as mp

def rows(al, N):
    be, ga = -al//6, -al//3
    a=[Fr(1),Fr(ga)]; b=[Fr(0),Fr(1)]
    for n in range(1,N):
        P=al*n*n+be*n+ga; Q=-27*n*n+63*n-30
        a.append((P*a[n]-Q*a[n-1])/Fr((n+1)**2))
        b.append((P*b[n]-Q*b[n-1])/Fr((n+1)**2))
    return a,b

def xi_of(al, N=700, dps=200):
    a,b=rows(al,N)
    W=[a[m]*b[m+1]-a[m+1]*b[m] for m in range(N)]
    S=sum(W[m]/(a[m]*a[m+1]) for m in range(N-1))
    mp.mp.dps=dps
    return mp.mpf(S.numerator)/mp.mpf(S.denominator)

if __name__=="__main__":
    print("alpha with an integral row, alpha <= 3000:")
    good=[al for al in range(6,3001,6) if all(x.denominator==1 for x in rows(al,60)[0])]
    print("  ", good[:12], "...   all == 18 mod 36:", all(al%36==18 for al in good))
    mp.mp.dps=190; s3=mp.sqrt(3)
    D=[1]*702
    for n in range(1,702): D[n]=D[n-1]*n//gcd(D[n-1],n)
    print("\n alpha    score       xi                                    "
          "-atan(6r3/al)/(2r3) - xi     k")
    for al in (18,54,90,126,234,306,558,846):
        x=xi_of(al)
        l2=(al-sqrt(al*al+108))/2
        a,b=rows(al,320)
        k=0
        while k<=6 and not all((D[n]**k*b[n]).denominator==1 for n in range(321)): k+=1
        print("%6d  %+8.5f   %s   %10s   %d"%(
            al, -log(abs(l2))-2, mp.nstr(x,30),
            mp.nstr(-mp.atan(6*s3/al)/(2*s3)-x,4), k))
