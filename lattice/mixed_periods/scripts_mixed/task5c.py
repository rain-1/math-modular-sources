import mpmath as mp
WD='/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/mixed/'
out=[]
def P(s):
    print(s); out.append(str(s))

def DBW(z):
    z=mp.mpc(z)
    return mp.im(mp.polylog(2,z)) + mp.arg(1-z)*mp.log(abs(z))
def CM(a,b,c):
    a=mp.mpf(a); b=mp.mpf(b); c=mp.mpf(c)
    if a<b+c and b<a+c and c<a+b:
        al=mp.acos((b*b+c*c-a*a)/(2*b*c)); be=mp.acos((a*a+c*c-b*b)/(2*a*c)); ga=mp.acos((a*a+b*b-c*c)/(2*a*b))
        return al*mp.log(a)+be*mp.log(b)+ga*mp.log(c)+DBW((b/a)*mp.exp(1j*ga))
    return mp.pi*mp.log(max(a,b,c))
def mahler_num(a,b,c,dps=40):
    with mp.workdps(dps+15):
        a=mp.mpf(a);b=mp.mpf(b);c=mp.mpf(c)
        def g(s):
            A=mp.sqrt(a*a+b*b+2*a*b*mp.cos(s)); return mp.log(A) if A>c else mp.log(c)
        t=(c*c-a*a-b*b)/(2*a*b); pts=[mp.mpf(0),mp.pi]
        if -1<t<1: pts=[mp.mpf(0),mp.acos(t),mp.pi]
        return +(mp.quad(g,pts)/mp.pi)
def dig(x,y):
    d=abs(x-y); s=max(abs(x),abs(y),mp.mpf(1))
    return mp.inf if d==0 else -mp.log10(d/s)

mp.mp.dps=120
P("")
P("="*78)
P("PART B -- structure of Q_m as a Mahler measure")
P("="*78)

P("")
P("B0. The m=1 control (D=3).")
m=1; D=3; A=mp.sqrt(D); th=2*mp.atan(1/A)
Q1=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/m)
P("   theta = 2*atan(1/sqrt3) = %s   ; pi/3 = %s" % (mp.nstr(th,30), mp.nstr(mp.pi/3,30)))
P("   theta - pi/3 = %s  (theta = pi/3 confirmed to %.0f digits)"%(mp.nstr(th-mp.pi/3,5), float(dig(th,mp.pi/3))))
P("   Q_1 = 2*Cl2(2pi/3) - (pi/3)*log 3 = %s" % mp.nstr(Q1,40))
M1=CM(1,1,1)
P("   M = pi*m(1+x+y)                  = %s" % mp.nstr(M1,40))
P("   Cl2(pi/3)                        = %s  (agree %.0f digits)"%(mp.nstr(mp.clsin(2,mp.pi/3),40),float(dig(M1,mp.clsin(2,mp.pi/3)))))
L2=(mp.zeta(2,mp.mpf(1)/3)-mp.zeta(2,mp.mpf(2)/3))/9
P("   (3sqrt3/4)*L(2,chi_-3)           = %s  (agree %.0f digits)"%(mp.nstr(3*mp.sqrt(3)/4*L2,40),float(dig(M1,3*mp.sqrt(3)/4*L2))))
r=mp.pslq([M1,Q1,th,A,mp.mpf(1)],maxcoeff=10**6,maxsteps=10**6,tol=mp.mpf(10)**-45)
P("   PSLQ [M, Q_1, theta, sqrt3, 1], maxcoeff 1e6, tol 1e-45 : %s"%r)
r=mp.pslq([M1,Q1,th,A,mp.mpf(1)],maxcoeff=10**4,maxsteps=10**6,tol=mp.mpf(10)**-45)
P("   PSLQ [M, Q_1, theta, sqrt3, 1], maxcoeff 1e4            : %s"%r)
P("   -> NO relation in that basis, because Q_1 contains theta*log 3 and log 3 is")
P("      not in the span of {M, theta, sqrt3, 1}.  Adjoin pi*log 3:")
r=mp.pslq([M1,Q1,mp.pi*mp.log(3)],maxcoeff=10**4,maxsteps=10**6,tol=mp.mpf(10)**-80)
P("   PSLQ [M, Q_1, pi*log3], tol 1e-80 : %s"%r)
lhs=4*M1-3*Q1-mp.pi*mp.log(3)
P("   EXACT RELATION (m=1):   4*M - 3*Q_1 - pi*log 3 = 0")
P("      residual = %s   (i.e. 0 to %.0f digits at dps=120)"%(mp.nstr(lhs,5),float(-mp.log10(abs(lhs)) if lhs!=0 else 200)))
P("      equivalently  Q_1 = (4/3)Cl2(pi/3) - (pi/3)log3 = (4/3)*pi*m(1+x+y) - (pi/3)*log 3,")
P("      and  pi*m(1+x+y) = (3/4)Q_1 + (pi/4)log 3 = (3sqrt3/4)L(2,chi_-3).")

P("")
P("B1. General identity found analytically and verified numerically.")
P("   Cassaigne-Maillot for the ISOSCELES triple (1,1,t) with t=sqrt(D/m) gives")
P("   the apex angle gamma = pi - theta  (since cos gamma = 1 - t^2/2 = -(2m-1)/(2m)")
P("   and cos theta = (D-1)/(D+1) = (2m-1)/(2m)), hence")
P("        pi*m(1+x+sqrt(D/m)*y) = (1/2)(pi-theta)*log(D/m) + Cl2(pi-theta),")
P("   and therefore")
P("        Q_m = 2*pi*m(1 + x + sqrt(D/m)*y) - pi*log(D/m).")
P("")
P("   %-4s %-5s %-22s %-42s %s"%("m","D","sqrt(D/m)","Q_m","residual of 2*pi*m(1,1,t)-pi*log(D/m)-Q_m"))
for m in [1,2,3,5,7,11,25]:
    D=4*m-1; A=mp.sqrt(D); th=2*mp.atan(1/A)
    Q=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/m)
    t=mp.sqrt(mp.mpf(D)/m)
    Mt=CM(1,1,t)
    res=2*Mt-mp.pi*mp.log(mp.mpf(D)/m)-Q
    P("   %-4d %-5d %-22s %-42s %s"%(m,D,mp.nstr(t,15),mp.nstr(Q,30),mp.nstr(res,4)))
# independent check against raw numerical Mahler measure for m=2
for m in [2,3,5]:
    D=4*m-1; t=mp.sqrt(mp.mpf(D)/m)
    mn=mahler_num(1,1,t,dps=40); cm=CM(1,1,t)/mp.pi
    P("   independent torus-integral check m(1+x+%s y): CM vs quadrature agree %.0f digits"%(mp.nstr(t,10),float(dig(mn,cm))))

P("")
P("B2. Scale-equivalent integral-square form: sides (sqrt(m), sqrt(m), sqrt(D)).")
for m in [2,3,5]:
    D=4*m-1
    Mv=CM(mp.sqrt(m),mp.sqrt(m),mp.sqrt(D))
    A=mp.sqrt(D); th=2*mp.atan(1/A)
    Q=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/m)
    res=2*Mv-mp.pi*mp.log(mp.mpf(D))-Q   # since m(lam P)=log lam + m(P), lam=sqrt(m)
    P("   m=%d: 2*pi*m(sqrt%d + sqrt%d x + sqrt%d y) - pi*log(%d) - Q_%d = %s"%(m,m,m,D,D,m,mp.nstr(res,4)))

open(WD+'task5_results.txt','a').write("\n".join(out)+"\n")
