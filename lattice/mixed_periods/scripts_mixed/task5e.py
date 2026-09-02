import mpmath as mp
WD='/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/mixed/'
out=[]
def P(s):
    print(s); out.append(str(s))
def DBW(z):
    z=mp.mpc(z); return mp.im(mp.polylog(2,z))+mp.arg(1-z)*mp.log(abs(z))
def CM(a,b,c):
    a=mp.mpf(a);b=mp.mpf(b);c=mp.mpf(c)
    if a<b+c and b<a+c and c<a+b:
        al=mp.acos((b*b+c*c-a*a)/(2*b*c));be=mp.acos((a*a+c*c-b*b)/(2*a*c));ga=mp.acos((a*a+b*b-c*c)/(2*a*b))
        return al*mp.log(a)+be*mp.log(b)+ga*mp.log(c)+DBW((b/a)*mp.exp(1j*ga))
    return mp.pi*mp.log(max(a,b,c))

def mahler_Pm(m,D,dps=50):
    """m( m*(1+x)^2 - D*y^2 ) by Jensen in y then 1-D quadrature in x."""
    with mp.workdps(dps+15):
        m=mp.mpf(m); D=mp.mpf(D); r=mp.sqrt(m/D)
        # roots in y have modulus r*|1+x|, two of them; leading coeff -D
        def g(s):
            u=r*mp.sqrt(2+2*mp.cos(s))     # r*|1+e^{is}|
            return 2*mp.log(u) if u>1 else mp.mpf(0)
        # split where r*|1+e^{is}| = 1  <=> cos s = 1/(2r^2)-1
        t=1/(2*r*r)-1; pts=[mp.mpf(0),mp.pi]
        if -1<t<1: pts=[mp.mpf(0),mp.acos(t),mp.pi]
        return +(mp.log(D)+mp.quad(g,pts)/mp.pi)

mp.mp.dps=120
P("")
P("B3. Integer-coefficient reformulation.")
P("   Since m(1+x-t y)=m(1+x+t y), 2*m(1+x+t y)=m((1+x)^2-t^2 y^2); with t^2=D/m this gives")
P("        Q_m = pi*[ m( m*(1+x)^2 - D*y^2 ) - log D ]     (integer polynomial!)")
P("   %-4s %-5s %-46s %s"%("m","D","m(m(1+x)^2 - D y^2)  [torus quadrature]","pi*(that) - pi*log D - Q_m"))
for mm in [1,2,3,5,7,11]:
    D=4*mm-1; A=mp.sqrt(D); th=2*mp.atan(1/A)
    Q=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/mm)
    val=mahler_Pm(mm,D,dps=50)
    res=mp.pi*val-mp.pi*mp.log(D)-Q
    P("   %-4d %-5d %-46s %s"%(mm,D,mp.nstr(val,35),mp.nstr(res,4)))
P("   (quadrature carried at ~50 working digits, so residuals ~1e-50 confirm the identity")
P("    to that accuracy; the CM-based check in B1 confirms it to 120 digits.)")
open(WD+'task5_results.txt','a').write("\n".join(out)+"\n")
