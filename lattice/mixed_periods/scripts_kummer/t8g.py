"""Final Task 8: exact real closed form
   c^{(a)}[log(1-t)/(1-t)] = sum_z [ beta_z Re Li2(z) + gamma_z D(z) ] + elementary """
from ident2 import *
import sys
setup(300)
G=mp.catalan; CL=im(polylog(2,exp(mpc(0,1)*pi/3)))
def orb_key(z, reps, tol):
    """match z to a representative modulo z->1-z (sign -1) and z->conj (Re:+1, D:-1)"""
    for i,p in enumerate(reps):
        for zz,sR,sD in [(z,1,1),(mp.conj(z),1,-1),(1-z,-1,-1),(1-mp.conj(z),-1,1)]:
            if abs(zz-p)<tol: return i,sR,sD
    return None,None,None
HOSTS=[(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P')]
for (k,N,fam) in HOSTS:
    h=Host(k,N,fam); n=int(mp.nint(mpf(h.D)**(mpf(1)/k)))
    sq = 3 if k==3 else None
    names,vals,P,AN=full_basis(h,sqrtd=sq)
    print(f"\n##### k={k}  H=(1{'-' if fam=='M' else '+'}{N}x)^(-1/{k})  D={h.D}={n}^{k}  logs{P} args{AN} #####")
    for a in range(1,k):
        Pt=parts(h,a); tol=mpf(10)**-200
        reps=[]; bet=[]; gam=[]
        for z,l in zip(Pt['pts'],Pt['lams']):
            i,sR,sD = orb_key(z,reps,tol)
            if i is None:
                reps.append(z); bet.append(mpf(0)); gam.append(mpf(0)); i=len(reps)-1; sR=sD=1
            bet[i]+= sR*re(l); gam[i]+= -sD*im(l)
        # elementary leftovers from folding: use exact reflection/conjugation identities
        # Re Li2(1-z) = pi^2/6 - log|z|log|1-z| + arg z arg(1-z) - Re Li2(z)
        # D(1-z) = -D(z);  Re Li2(zbar)=Re Li2(z); D(zbar)=-D(z)
        corr=mpf(0)
        for z,l in zip(Pt['pts'],Pt['lams']):
            i,sR,sD = orb_key(z,reps,tol)
            if sR==-1:  # z folded as 1-w with w=reps[i]; ReLi2(z)= [that identity] - ReLi2(w)
                corr += re(l)*( pi**2/6 - log(abs(1-z))*log(abs(z)) + arg(1-z)*arg(z) )
        rebuilt = sum(b*re(polylog(2,zz)) for b,zz in zip(bet,reps)) \
                + sum(g*BW(zz) for g,zz in zip(gam,reps))
        Elem = re(h.cD(a)) - rebuilt   # includes the folding correction 'corr'
        print(f"\n  --- a={a} ---   c^({a})_D = {mp.nstr(re(h.cD(a)),40)}")
        for b,g,z in zip(bet,gam,reps):
            if abs(b)>mpf(10)**-100 or abs(g)>mpf(10)**-100:
                rb=as_rat(b,10**7); rb3=as_rat(b/sqrt(3),10**7)
                rg=as_rat(g,10**7); rg3=as_rat(g/sqrt(3),10**7)
                bs = str(rb) if rb is not None else (str(rb3)+"*sqrt3" if rb3 is not None else mp.nstr(b,15))
                gs = str(rg) if rg is not None else (str(rg3)+"*sqrt3" if rg3 is not None else mp.nstr(g,15))
                print(f"      z={mp.nstr(z,16):38s} ReLi2 coeff {bs:16s} D coeff {gs:16s}  D(z)={mp.nstr(BW(z),16)}")
        d = 300-2*len(names)
        r=find(Elem,names,vals,digits=d,maxcoeff=10**6)
        print(f"      elementary part {mp.nstr(Elem,25)}: "+(" + ".join(f"({t[0]})*{t[1]}" for t in r[0]) if r else f"NO RELATION (basis {len(names)}, {d} digits, maxcoeff 1e6)"))
        sys.stdout.flush()
