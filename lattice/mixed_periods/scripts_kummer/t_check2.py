from kummer import *
setup(120)
tests = [(2,4,'M'),(2,4,'P'),(3,27,'M'),(3,27,'P'),(3,9,'M'),(3,63,'P'),(3,126,'M'),(4,256,'M'),(4,80,'P'),(4,624,'P'),(5,3125,'M'),(5,1025,'P'),(5,7775,'M'),(6,15624,'P'),(6,46656,'M')]
worst={'ut':0,'cf':0}
for (k,N,fam) in tests:
    h = Host(k,N,fam)
    for a in range(1,k):
        for kern in ['B','D','L']:
            u = h.per_u(a,kern); t = h.per_t(a,kern)
            cf = {'B':h.cB,'D':h.cD,'L':h.cL}[kern](a)
            e1 = abs(u-t)/max(abs(u),1e-99); e2 = abs(u-cf)/max(abs(u),1e-99)
            worst['ut']=max(worst['ut'],float(e1)); worst['cf']=max(worst['cf'],float(e2))
            flag = "" if e2 < mpf(10)**-90 else "   <<< CLOSED FORM MISMATCH"
            print(f"k={k} N={N} {fam} a={a} {kern}: val={mp.nstr(u,30)} relerr(t-quad)={mp.nstr(e1,3)} relerr(closed)={mp.nstr(e2,3)}{flag}")
print("worst rel err u-vs-t:",worst['ut']," u-vs-closedform:",worst['cf'])
