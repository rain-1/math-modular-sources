from ident2 import *
import sys
setup(160)
print("TASK 5: period span per host and sector (PSLQ at 140 digits, maxcoeff 1e12)\n")
KERNS=['B','D','L','tB','tD']
for k in [3,4]:
    for (N,fam) in ([(9,'M'),(27,'M'),(27,'P'),(63,'P')] if k==3 else [(8,'M'),(80,'P'),(256,'M'),(256,'P')]):
        h=Host(k,N,fam)
        print(f"--- k={k} N={N} {fam} D={h.D} ---")
        allper={}
        for a in range(1,k):
            v=[h.per_u(a,x) for x in KERNS]
            allper[a]=v
            r=pslq([mpf(1)]+v,tol=mpf(10)**-140,maxcoeff=10**12,maxsteps=500000)
            print(f"   a={a}: 6-term PSLQ {{1,B,D,L,tB,tD}} -> {r}")
            # is every period in Q + Q*B + Q*D ?
            for i,nm in enumerate(KERNS):
                if nm in ('B','D'): continue
                rr=pslq([v[i],mpf(1),v[0],v[1]],tol=mpf(10)**-140,maxcoeff=10**12,maxsteps=500000)
                print(f"        {nm} in Q+Q*B+Q*D ? {rr}")
        # cross-sector
        for a in range(1,k):
            for b in range(a+1,k):
                for i,nm in enumerate(KERNS[:2]):
                    rr=pslq([allper[b][i],mpf(1),allper[a][i]],tol=mpf(10)**-140,maxcoeff=10**12,maxsteps=500000)
                    print(f"   cross-sector {nm}: c^({b}) in Q + Q c^({a}) ? {rr}")
        sys.stdout.flush()
