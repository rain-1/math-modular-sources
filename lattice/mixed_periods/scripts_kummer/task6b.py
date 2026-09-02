from ident2 import *
import sys, itertools
setup(150)
CL=im(polylog(2,exp(mpc(0,1)*pi/3))); G=mp.catalan
def CMparts(a,b,c):
    a,b,c=mpf(a),mpf(b),mpf(c)
    cl=lambda x: mpf(-1) if x<-1 else (mpf(1) if x>1 else x)
    al=mp.acos(cl((b*b+c*c-a*a)/(2*b*c))); be=mp.acos(cl((a*a+c*c-b*b)/(2*a*c))); ga=mp.acos(cl((a*a+b*b-c*c)/(2*a*b)))
    z=(b/a)*exp(mpc(0,1)*ga)
    return z, re(BW(z)), al*log(a)+be*log(b)+ga*log(c)
print("TASK 6 (revised)\n(a) Bloch part vs pi*m(1+x+y)=Cl2(pi/3), allowing sqrt3-rational coefficients:")
for (k,N,fam) in [(3,9,'M'),(3,63,'P'),(3,126,'M')]:
    h=Host(k,N,fam)
    for a in range(1,k):
        B=parts(h,a)['Bloch']
        r=rat_ratio(B,CL,10**10); r3=rat_ratio(B,sqrt(3)*CL,10**10)
        print(f"   k=3 N={N} a={a}: Bloch = "+(f"{r} * pi*m(1+x+y)" if r else (f"{r3} * sqrt3 * pi*m(1+x+y)   <<< HIT" if r3 else "not a (sqrt3-)rational multiple of pi*m(1+x+y)")))
print("\n(b) triangle with vertices 1, n*zeta_l, n*zeta_j  ->  Cassaigne-Maillot point of the surviving Bloch points:")
for (k,N,fam) in [(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P')]:
    h=Host(k,N,fam); al=h.alphas(); n=mpf(h.D)**(mpf(1)/k)
    print(f"  --- k={k} N={N} {fam} D={h.D}, n={int(mp.nint(n))} ---")
    seen=set()
    for j in range(k):
        for l in range(k):
            if j==l: continue
            z=(1-al[l])/(al[j]-al[l]); d=BW(z)
            if abs(d)<mpf(10)**-60: continue
            A=abs(al[l]-al[j]); Bb=abs(al[l]-1); C=abs(al[j]-1)
            if A>Bb+C or Bb>A+C or C>A+Bb: continue
            zc,dc,lg = CMparts(A,Bb,C)
            ok=None
            for s,w in [(1,zc),(-1,1-zc),(1,mp.conj(zc)),(-1,1-mp.conj(zc)),(1,1/zc),(-1,1/(1-zc))]:
                if abs(w-z)<mpf(10)**-60: ok=s;break
            key=(mp.nstr(A,12),mp.nstr(Bb,12),mp.nstr(C,12))
            if key in seen: continue
            seen.add(key)
            print(f"     triangle sides ({mp.nstr(A,10)}, {mp.nstr(Bb,10)}, {mp.nstr(C,10)}):  CM point {mp.nstr(zc,10)},"
                  f"  D(z_period)={mp.nstr(d,14)}  D(z_CM)={mp.nstr(dc,14)}   match: {ok if ok else 'NO'}"
                  f"   pi*m = {mp.nstr(dc+lg,16)}")
    sys.stdout.flush()
