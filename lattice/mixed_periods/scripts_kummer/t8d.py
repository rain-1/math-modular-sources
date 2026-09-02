from ident2 import *
import sys
setup(160)
G=mp.catalan; CL=im(polylog(2,exp(mpc(0,1)*pi/3)))
HOSTS=[(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P')]
def fmt(terms): return " + ".join(f"({t[0]})*{t[1]}" for t in terms)
for (k,N,fam) in HOSTS:
    h=Host(k,N,fam); n=int(mp.nint(mpf(h.D)**(mpf(1)/k)))
    REF, RN = (CL,"Cl2(pi/3)") if k==3 else (G,"G")
    sq = 3 if k==3 else None
    names,vals,P,AN = full_basis(h,sqrtd=sq)
    print(f"\n########## k={k}  H=(1{'-' if fam=='M' else '+'}{N}x)^(-1/{k})  D={h.D}={n}^{k} ##########")
    print(f"  log-primes {P};  arg generators {AN};  elementary basis size {len(names)}")
    sys.stdout.flush()
    for a in range(1,k):
        Pt=parts(h,a); reps=[]; refc=mpf(0); coefs=[]
        for z,l in zip(Pt['pts'],Pt['lams']):
            c=-im(l); d=BW(z)
            if abs(d)<mpf(10)**-60 or abs(c)<mpf(10)**-60: continue
            q=rat_ratio(d,REF)
            if q is not None: refc+=c*mpf(q.numerator)/q.denominator; continue
            hit=False
            for i,(zr,dr) in enumerate(reps):
                q2=rat_ratio(d,dr)
                if q2 is not None: coefs[i]+=c*mpf(q2.numerator)/q2.denominator; hit=True; break
            if not hit: reps.append((z,d)); coefs.append(c)
        Bl=Pt['Bloch']; cd=re(h.cD(a)); Elem=cd-Bl
        print(f"\n  --- a={a} ---  c^({a})[log(1-t)/(1-t)] = {mp.nstr(cd,45)}")
        rc = as_rat(refc,10**7) if k==4 else as_rat(refc/sqrt(3),10**7)
        lbl = RN if k==4 else f"sqrt3*{RN}"
        print(f"      coeff of {lbl}: {rc if rc is not None else mp.nstr(refc,25)+' (irrational?)'}")
        for cc,(zz,dd) in zip(coefs,reps):
            if abs(cc)>mpf(10)**-60:
                r1=as_rat(cc,10**7); r2=as_rat(cc/sqrt(3),10**7)
                print(f"      + [{r1 if r1 is not None else (str(r2)+'*sqrt3' if r2 is not None else mp.nstr(cc,20))}] * D({mp.nstr(zz,14)})   (D={mp.nstr(dd,18)})")
        res=find(Elem,names,vals,digits=110,maxcoeff=10**6)
        if res: print(f"      + elementary: {fmt(res[0])}    [residual {mp.nstr(fabs(res[1]),3)}]")
        else:   print(f"      + elementary remainder {mp.nstr(Elem,30)}  *** NOT in the basis ***")
        sys.stdout.flush()
