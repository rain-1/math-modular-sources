from ident2 import *
import sys
setup(320)
G=mp.catalan; CL=im(polylog(2,exp(mpc(0,1)*pi/3)))
print("DECISIVE TESTS at dps=320\n")
# ---- (i) P0 cancellation theorem, numerically ----
print("(i) aggregate Li2-coefficient of each cyclotomic P0 point (should be exactly 0):")
for (k,N,fam) in [(3,9,'M'),(3,63,'P'),(4,80,'P'),(4,624,'P'),(5,1025,'M'),(6,15624,'P')]:
    h=Host(k,N,fam)
    for a in range(1,k):
        al=h.alphas(); n=mpf(h.D)**(mpf(1)/k)
        agg={}
        for j in range(k):
            for l in range(k):
                if j==l: continue
                z=-al[l]/(al[j]-al[l]); c=(mpf(h.sigma)/h.D)*al[j]**a
                key=None
                for kk in agg:
                    if abs(agg[kk][0]-z)<mpf(10)**-200: key=kk;break
                if key is None: agg[len(agg)]=[z,c]
                else: agg[key][1]+=c
        mx=max(abs(v[1]) for v in agg.values())
        print(f"    k={k} N={N} {fam} a={a}: max |aggregate coeff over P0 points| = {mp.nstr(mx,3)}")
sys.stdout.flush()
# ---- (ii) is the Bloch part a rational multiple of G / Cl2(pi/3)? ----
print("\n(ii) Bloch part vs the reference constant (300-digit PSLQ, 2-term):")
for (k,N,fam) in [(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,80,'P'),(4,624,'P')]:
    h=Host(k,N,fam); REF,RN=(CL,"Cl2(pi/3)") if k==3 else (G,"G")
    for a in range(1,k):
        P=parts(h,a); B=P['Bloch']
        r=rat_ratio(B,REF,10**12); r3=rat_ratio(B,sqrt(3)*REF,10**12)
        s = f"= {r} * {RN}" if r else (f"= {r3} * sqrt3 * {RN}" if r3 else f"NOT a (rational or sqrt3-rational) multiple of {RN}")
        print(f"    k={k} N={N} {fam} a={a}: Bloch = {mp.nstr(B,25)}  {s}")
sys.stdout.flush()
# ---- (iii) is Re-Li2 part elementary?  tight basis, 300 digits ----
print("\n(iii) is the Re Li2 part elementary?  (tight basis: pi^2, log p log q, [x sqrt3])")
for (k,N,fam) in [(3,9,'M'),(4,80,'P')]:
    h=Host(k,N,fam)
    LB,P=logbasis(h)
    for sq in ([None,3] if k==3 else [None]):
        names=["pi^2"];vals=[pi**2]
        for i in range(len(LB)):
            for j in range(i,len(LB)):
                names.append(LB[i][0]+"*"+LB[j][0]); vals.append(LB[i][1]*LB[j][1])
        if sq:
            b=list(zip(names,vals))
            for nm,v in b: names.append(f"sqrt3*{nm}"); vals.append(sqrt(3)*v)
        for a in range(1,k):
            Pt=parts(h,a)
            for lbl,tgt in [("Re Li2 part",Pt['LiRe']),("Elem = c_D - Bloch",re(h.cD(a))-Pt['Bloch'])]:
                res=find(tgt,names,vals,digits=260,maxcoeff=10**12)
                print(f"    k={k} N={N} a={a} sqrt3={bool(sq)} {lbl}: "
                      +(" + ".join(f"({t[0]})*{t[1]}" for t in res[0]) if res else "NO RELATION (basis size %d, 260 digits, maxcoeff 1e12)"%len(names)))
        sys.stdout.flush()
