from ident2 import *
import sys
setup(320)
G=mp.catalan; CL=im(polylog(2,exp(mpc(0,1)*pi/3)))
def fmt(t): return " + ".join(f"({x[0]})*{x[1]}" for x in t)
print("(iv) Elem = c_D - Bloch  against the FULL log+arg product basis, 300 digits, maxcoeff 1e10\n")
for (k,N,fam) in [(3,9,'M'),(3,63,'P'),(4,80,'P'),(4,624,'P')]:
    h=Host(k,N,fam)
    sq = 3 if k==3 else None
    names,vals,P,AN=full_basis(h,sqrtd=sq)
    print(f"  k={k} N={N} {fam}: logs {P}, args {AN}, basis size {len(names)}")
    for a in range(1,k):
        Pt=parts(h,a); El=re(h.cD(a))-Pt['Bloch']
        r=find(El,names,vals,digits=300-2*len(names),maxcoeff=10**10)
        print(f"     a={a} Elem: "+(fmt(r[0])[:300] if r else "NO RELATION"))
        r2=find(Pt['LiRe'],names,vals,digits=300-2*len(names),maxcoeff=10**10)
        print(f"     a={a} LiRe: "+(fmt(r2[0])[:300] if r2 else "NO RELATION"))
    sys.stdout.flush()
print("\n(v) relations among the surviving Bloch points and G / Cl2(pi/3), 300 digits")
for (k,N,fam) in [(4,80,'P'),(4,624,'P'),(3,63,'P'),(3,126,'M')]:
    h=Host(k,N,fam); REF,RN=(CL,"Cl2(pi/3)") if k==3 else (G,"G")
    Pt=parts(h,1); seen=[]
    for z in Pt['pts']:
        d=BW(z)
        if abs(d)<mpf(10)**-100: continue
        if any(abs(abs(d)-abs(x[1]))<mpf(10)**-100 for x in seen): continue
        seen.append((z,d))
    v=[x[1] for x in seen]+[REF]
    r=pslq(v,tol=mpf(10)**-280,maxcoeff=10**8,maxsteps=500000)
    print(f"  k={k} N={N}: {len(seen)} independent |D| values; PSLQ with {RN}: {r}")
    for (z,d) in seen: print(f"       D({mp.nstr(z,12)}) = {mp.nstr(d,25)}   / {RN} = {mp.nstr(d/REF,20)}")
