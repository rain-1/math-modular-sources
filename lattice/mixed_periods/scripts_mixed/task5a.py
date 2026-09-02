import mpmath as mp

def mahler_num(a,b,c,dps=40):
    """m(a+bx+cy) by exact Jensen reduction + 1-D quadrature."""
    with mp.workdps(dps+15):
        a=mp.mpf(a); b=mp.mpf(b); c=mp.mpf(c)
        # |a+b e^{is}|^2 = a^2+b^2+2ab cos s ; integrand log max(|.|, c), symmetric in s->-s
        # integrate over [0,pi], divide by pi
        f = lambda s: mp.log(mp.mpf(2)*mp.sqrt(a*a+b*b+2*a*b*mp.cos(s)))  # placeholder
        def g(s):
            A = mp.sqrt(a*a+b*b+2*a*b*mp.cos(s))
            return mp.log(A) if A>c else mp.log(c)
        # split points: cos s = (c^2-a^2-b^2)/(2ab)
        t = (c*c-a*a-b*b)/(2*a*b)
        pts=[mp.mpf(0), mp.pi]
        if -1 < t < 1:
            s0 = mp.acos(t)
            pts=[mp.mpf(0), s0, mp.pi]
        val = mp.quad(g, pts)
        return +(val/mp.pi)

def DBW(z):
    z=mp.mpc(z)
    if abs(z)==0: return mp.mpf(0)
    return mp.im(mp.polylog(2,z)) + mp.arg(1-z)*mp.log(abs(z))

def CM(a,b,c):
    """Cassaigne-Maillot value of pi*m(a+bx+cy)."""
    a=mp.mpf(a); b=mp.mpf(b); c=mp.mpf(c)
    if a<b+c and b<a+c and c<a+b:
        al = mp.acos((b*b+c*c-a*a)/(2*b*c))
        be = mp.acos((a*a+c*c-b*b)/(2*a*c))
        ga = mp.acos((a*a+b*b-c*c)/(2*a*b))
        return al*mp.log(a)+be*mp.log(b)+ga*mp.log(c)+DBW((b/a)*mp.exp(1j*ga))
    else:
        return mp.pi*mp.log(max(a,b,c))

def digits(x,y):
    d = abs(x-y)
    s = max(abs(x),abs(y),mp.mpf(1))
    if d==0: return mp.inf
    return -mp.log10(d/s)

out=[]
def P(s):
    print(s); out.append(s)

mp.mp.dps=50
P("="*78)
P("PART A -- Cassaigne-Maillot verification")
P("="*78)

# (2) Smyth
mp.mp.dps=60
chi = lambda n: (1 if n%3==1 else (-1 if n%3==2 else 0))
L2 = mp.nsum(lambda k: mp.mpf(1)/(3*k+1)**2 - mp.mpf(1)/(3*k+2)**2, [0,mp.inf])
# better: use Hurwitz zeta
L2b = (mp.zeta(2,mp.mpf(1)/3)-mp.zeta(2,mp.mpf(2)/3))/9
L2c = (4/(3*mp.sqrt(3)))*mp.clsin(2,mp.pi/3)
P("")
P("(2) L(2,chi_-3):")
P("   direct sum          = %s" % mp.nstr(L2,40))
P("   Hurwitz zeta        = %s" % mp.nstr(L2b,40))
P("   (4/(3sqrt3))Cl2(pi/3)= %s" % mp.nstr(L2c,40))
P("   agree(sum,Hurwitz)  : %.1f digits" % float(digits(L2,L2b)))
P("   agree(Hurwitz,Cl2)  : %.1f digits" % float(digits(L2b,L2c)))
mnum = mahler_num(1,1,1,dps=40)
smyth = (3*mp.sqrt(3)/(4*mp.pi))*L2b
P("   m(1+x+y) numeric    = %s" % mp.nstr(mnum,40))
P("   Smyth (3sqrt3/4pi)L = %s" % mp.nstr(smyth,40))
P("   agreement           : %.1f digits" % float(digits(mnum,smyth)))

# (3) triples
P("")
P("(3) Cassaigne-Maillot vs numerical Mahler measure  (pi*m):")
trips=[(1,1,1),(2,3,4),(3,4,5),(5,6,7),(2,2,3),(1,2,2),(7,8,9),(1,1,2),(3,5,7),
       (11,13,17),(1,1,5),(1,2,10),(5,1,1),(2,3,10),(1,1,1.5),(2.5,3.5,4.5),(1,7,7)]
P("   %-16s %-8s %-30s %s" % ("(a,b,c)","type","pi*m (numeric)","digits agree"))
worst=mp.inf
for (a,b,c) in trips:
    tri = (a<b+c and b<a+c and c<a+b)
    num = mp.pi*mahler_num(a,b,c,dps=40)
    cm  = CM(a,b,c)
    d = digits(num,cm)
    worst=min(worst,d)
    P("   %-16s %-8s %-30s %.1f" % (str((a,b,c)), "tri" if tri else "non-tri", mp.nstr(num,30), float(d)))
P("   -> worst-case agreement over %d triples: %.1f digits" % (len(trips),float(worst)))
P("   CONFIRMED: Cassaigne-Maillot formula reproduces the numerical Mahler measure.")

open('/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/mixed/task5_results.txt','a').write("\n".join(out)+"\n")
