from mpmath import mp, mpf, sqrt, log, pi, quad, clsin, atan, nstr, fabs, cos, sin, acos, mpc, polylog, arg, findroot
mp.dps=60
def mahler_lin(a,b,c):
    """m(a+bx+cy) by Jensen in one variable + 1-D quadrature (exact reduction)."""
    a=mpf(a);b=mpf(b);c=mpf(c)
    f=lambda s: log(max(abs(a+b*mpc(cos(s),sin(s))), c))
    # split where |a+b e^{is}| = c  ->  a^2+b^2+2ab cos s = c^2
    pts=[mpf(0),pi]
    val=(c*c-a*a-b*b)/(2*a*b)
    if -1<val<1: pts=[mpf(0),acos(val),pi]
    return quad(f,pts)/pi
print("m | D' | 2*pi*m(1+x+sqrt(D/m)y) - pi*log(D/m)  vs  Q_m   [independent 1-D torus quadrature]")
for m in [1,2,3,5,7,11]:
    D=4*m-1; a=sqrt(D); th=2*atan(1/a)
    Q=2*clsin(2,pi-th)-th*log(mpf(D)/m)
    t=sqrt(mpf(D)/m)
    M=mahler_lin(1,1,t)
    lhs=2*pi*M-pi*log(mpf(D)/m)
    print("%2d | %3d | %s   Q_m=%s   diff=%s"%(m,D,nstr(lhs,30),nstr(Q,30),nstr(fabs(lhs-Q),3)))
print()
print("integer form: Q_m = pi*[ m(m(1+x)^2 - D'y^2) - log D' ]  (2-D torus integral, 20 digits)")
mp.dps=30
def mahler2(P):
    f=lambda s,t: log(abs(P(mpc(cos(s),sin(s)),mpc(cos(t),sin(t)))))
    return quad(f,[0,2*pi],[0,2*pi])/(2*pi)**2
for m in [1,2,3]:
    D=4*m-1; a=sqrt(D); th=2*atan(1/a); Q=2*clsin(2,pi-th)-th*log(mpf(D)/m)
    P=lambda X,Y,m=m,D=D: m*(1+X)**2 - D*Y**2
    M2=mahler2(P)
    print("  m=%d: m(%d(1+x)^2-%dy^2)=%s   pi*(that-log D')-Q_m = %s"%(m,m,D,nstr(M2,22),nstr(pi*(M2-log(mpf(D)))-Q,3)))
