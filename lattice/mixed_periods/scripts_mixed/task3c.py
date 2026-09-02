from mpmath import mp, mpf, sqrt, log, atan, pi, quad, clsin, pslq, nstr, fabs
from fractions import Fraction
mp.dps = 150
def period(m,i,j,l):
    D = 4*m-1
    def f(u):
        t=(1-u*u)/(4*m); om=(D+u*u)/mpf(4*m)
        return t**i*om**(-j)*log(om)**l
    return quad(f,[0,1])/(2*m)
print("== independence of {1,c_B,c_D} (PSLQ, dps=150, maxcoeff 10^30) ==")
for m in [1,2,3,5,11]:
    D=4*m-1;a=sqrt(D);th=2*atan(1/a)
    cB=th/a; cD=(th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
    r=pslq([mpf(1),cB,cD],maxcoeff=10**30,maxsteps=10**6,tol=mpf(10)**-130)
    print("  m=%2d: relation among {1,c_B,c_D}: %s"%(m,r))
print()
print("== closed form check: c[log(1-t)^2] == (4 - 2D'*c_B + D'*c_D)/m ? ==")
for m in [1,2,3,5,11,17]:
    D=4*m-1;a=sqrt(D);th=2*atan(1/a)
    cB=th/a; cD=(th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
    E0=period(m,0,0,2); pred=(4-2*D*cB+D*cD)/mpf(m)
    print("  m=%2d: E0=%s  pred-E0=%s"%(m,nstr(E0,30),nstr(fabs(E0-pred),3)))
print()
print("== l=2, j>=1 periods (the [n][n/2][n/3] layer): do they leave span{1,cB,cD}? ==")
for m in [1,2,3]:
    D=4*m-1;a=sqrt(D);th=2*atan(1/a)
    cB=th/a; cD=(th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
    print(" m=%d"%m)
    for j in [1,2,3]:
        c=period(m,0,j,2)
        r=pslq([c,mpf(1),cB,cD],maxcoeff=10**10,maxsteps=10**6,tol=mpf(10)**-120)
        print("   (0,%d,2): c=%s   rel in {1,cB,cD}: %s"%(j,nstr(c,30),r))
