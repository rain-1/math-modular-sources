from mpmath import mp, mpf, sqrt, log, atan, pi, quad, clsin, pslq, nstr, fabs
mp.dps = 120

def period(m,i,j,l):
    D = 4*m-1
    def f(u):
        t = (1-u*u)/(4*m)
        om = (D+u*u)/mpf(4*m)   # 1-t
        return t**i * om**(-j) * log(om)**l
    return quad(f,[0,1])/(2*m)

print("== Task 3: closed-form checks (dps=120) ==")
for m in [1,2,3,5,11]:
    D=4*m-1; a=sqrt(D); th=2*atan(1/a)
    cB=period(m,0,1,0); cD=period(m,0,1,1)
    cBc = th/a
    cDc = (th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
    print("m=%d D'=%d"%(m,D))
    print("   c_B      = %s"%nstr(cB,55))
    print("   theta/a  = %s   |diff|=%s"%(nstr(cBc,55), nstr(fabs(cB-cBc),3)))
    print("   c_D      = %s"%nstr(cD,55))
    print("   closed   = %s   |diff|=%s"%(nstr(cDc,55), nstr(fabs(cD-cDc),3)))
