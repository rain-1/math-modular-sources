from mpmath import mp, mpf, sqrt, log, atan, pi, quad, clsin, pslq, nstr, fabs, mpmathify
from fractions import Fraction
mp.dps = 150

def period(m,i,j,l):
    D = 4*m-1
    def f(u):
        t = (1-u*u)/(4*m); om = (D+u*u)/mpf(4*m)
        return t**i * om**(-j) * log(om)**l
    return quad(f,[0,1])/(2*m)

def rat(x, tol=mpf(10)**-100):
    from mpmath import pslq as P
    r = P([x, mpf(1)], maxcoeff=10**12, maxsteps=10**6, tol=tol)
    if r and r[0]!=0: return Fraction(-r[1], r[0])
    return None

for m in [1,2,3,5,11]:
    D=4*m-1; a=sqrt(D); th=2*atan(1/a)
    cB=th/a; cD=(th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
    print("="*70); print("m=%d, D'=%d"%(m,D))
    print("-- l=0 : express c(i,j,0) in {1, c_B} --")
    for i in range(3):
        for j in range(4):
            c=period(m,i,j,0)
            r=pslq([c, mpf(1), cB], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-110)
            if r and r[0]!=0:
                q1=Fraction(-r[1],r[0]); q2=Fraction(-r[2],r[0])
                print("   (i,j)=(%d,%d): c = %s + (%s)*c_B"%(i,j,q1,q2))
            else:
                print("   (i,j)=(%d,%d): NO RELATION in {1,c_B};  c=%s"%(i,j,nstr(c,30)))
    print("-- l=1 : express c(i,j,1) in {1, c_B, c_D} --")
    for i in range(3):
        for j in range(4):
            c=period(m,i,j,1)
            r=pslq([c, mpf(1), cB, cD], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-110)
            if r and r[0]!=0:
                q=[Fraction(-r[k],r[0]) for k in (1,2,3)]
                print("   (i,j)=(%d,%d): c = %s + (%s)*c_B + (%s)*c_D"%(i,j,q[0],q[1],q[2]))
            else:
                print("   (i,j)=(%d,%d): NO RELATION;  c=%s"%(i,j,nstr(c,30)))
    print("-- l=2 : c(i,0,2) for i=0,1,2 : test span {1,c_B,c_D} --")
    for i in range(3):
        c=period(m,i,0,2)
        r=pslq([c, mpf(1), cB, cD], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-110)
        print("   (i,j)=(%d,0): c=%s"%(i,nstr(c,40)), " rel in {1,cB,cD}:", r)
    # are the three l=2 periods related to each other mod {1,cB,cD}?
    E=[period(m,i,0,2) for i in range(3)]
    r=pslq([E[0],E[1],mpf(1),cB,cD], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-110)
    print("   rel(E0,E1,1,cB,cD):",r)
    r=pslq([E[0],E[2],mpf(1),cB,cD], maxcoeff=10**8, maxsteps=10**6, tol=mpf(10)**-110)
    print("   rel(E0,E2,1,cB,cD):",r)
