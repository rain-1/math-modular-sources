from ident import *
setup(140)
import sys
G = mp.catalan
CL = im(polylog(2, exp(mpc(0,1)*pi/3)))
HOSTS = [(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P'),(5,1025,'M'),(5,7775,'P'),(6,15624,'P')]
TOL = mpf(10)**-110
def ratmul(x,y,mc=10**7):
    if abs(y)<mpf(10)**-60: return None
    r=pslq([x,y],tol=TOL,maxcoeff=mc,maxsteps=100000)
    if r is None or r[0]==0: return None
    return Rational(-int(r[1]),int(r[0]))
for (k,N,fam) in HOSTS:
    h=Host(k,N,fam); n=int(mp.nint(mpf(h.D)**(mpf(1)/k)))
    REF = CL if k==3 else (G if k==4 else None); RNAME="Cl2(pi/3)" if k==3 else ("G" if k==4 else None)
    print(f"\n########## k={k}  H=(1{'-' if fam=='M' else '+'}{N}x)^(-1/{k})   D={h.D}={n}^{k}   primes={primeset(h)} ##########")
    for a in range(1,k):
        P=parts(h,a); pts=P['pts']; lams=P['lams']
        reps=[]; refc=mpf(0); coefs=[]
        for z,l in zip(pts,lams):
            c=-im(l); d=BW(z)
            if abs(d)<mpf(10)**-60: continue
            if abs(c)<mpf(10)**-60: continue
            done=False
            if REF is not None:
                q=ratmul(d,REF)
                if q is not None: refc+=c*mpf(q.p)/q.q; done=True
            if not done:
                for i,(zr,dr) in enumerate(reps):
                    q=ratmul(d,dr)
                    if q is not None: coefs[i]+=c*mpf(q.p)/q.q; done=True; break
                if not done: reps.append((z,d)); coefs.append(c)
        Bloch=P['Bloch']; cd=re(h.cD(a))
        chk = refc*(REF if REF is not None else 0) + sum(cc*dd for cc,(zz,dd) in zip(coefs,reps))
        Elem = cd - Bloch
        print(f"\n  --- sector a={a} ---   c^({a})[log(1-t)/(1-t)] = {mp.nstr(cd,40)}")
        print(f"      Bloch(D) part = {mp.nstr(Bloch,30)}   (regroup check {mp.nstr(abs(chk-Bloch),3)})")
        if REF is not None:
            q=ratmul(refc,mpf(1),10**9)
            print(f"      coefficient of {RNAME}: {mp.nstr(refc,25)}   ~ {q if q else '(not rational)'}"
                  + (f"  -> {RNAME} coefficient nonzero" if abs(refc)>mpf(10)**-60 else "   *** ZERO ***"))
        for cc,(zz,dd) in zip(coefs,reps):
            if abs(cc)>mpf(10)**-60:
                qq=ratmul(cc,mpf(1),10**9); qs=ratmul(cc, sqrt(3),10**9)
                print(f"      surviving non-{RNAME} point z={mp.nstr(zz,14)}  D(z)={mp.nstr(dd,18)}  coeff={mp.nstr(cc,20)}"
                      + (f" = {qq}" if qq else (f" = {qs}*sqrt3" if qs else "")))
        for sq in ([None,3] if k in (3,6) else [None]):
            names,vals,_=elem_basis(h,with_sqrt=sq)
            res=show_pslq(Elem,names,vals,tol_digits=95,maxcoeff=10**7)
            if res:
                terms,resid,r=res
                print(f"      elementary remainder = {mp.nstr(Elem,30)}")
                print(f"         = "+" + ".join(f"({t[0]})*{t[1]}" for t in terms)+f"     [residual {mp.nstr(fabs(resid),3)}]")
                break
        else:
            print(f"      elementary remainder = {mp.nstr(Elem,30)}   *** NO relation found in the elementary basis ***")
        sys.stdout.flush()
