import mpmath as mp, itertools, math, time, sys, pickle

WD='/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/mixed/'

def DBW(z):
    z=mp.mpc(z)
    return mp.im(mp.polylog(2,z)) + mp.arg(1-z)*mp.log(abs(z))

def CM(a,b,c):
    a=mp.mpf(a); b=mp.mpf(b); c=mp.mpf(c)
    if a<b+c and b<a+c and c<a+b:
        al = mp.acos((b*b+c*c-a*a)/(2*b*c))
        be = mp.acos((a*a+c*c-b*b)/(2*a*c))
        ga = mp.acos((a*a+b*b-c*c)/(2*a*b))
        return al*mp.log(a)+be*mp.log(b)+ga*mp.log(c)+DBW((b/a)*mp.exp(1j*ga))
    return mp.pi*mp.log(max(a,b,c))

def consts(m,dps):
    with mp.workdps(dps):
        D=4*m-1; A=mp.sqrt(D); th=2*mp.atan(1/A)
        Q=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/m)
        return +Q,+th,+A,D

N=int(sys.argv[1]) if len(sys.argv)>1 else 60
MS=[2,3,5]
mp.mp.dps=70

# sanity: no degenerate relation among the basis constants alone
print("=== basis sanity (no relation among Q_m,theta,sqrtD,1 alone) ===")
for m in MS:
    Q,th,A,D=consts(m,70)
    r=mp.pslq([Q,th,A,mp.mpf(1)],maxcoeff=10**6,maxsteps=10**5,tol=mp.mpf(10)**-45)
    print("  m=%d D=%d :"%(m,D), r)
    r2=mp.pslq([th,A,mp.mpf(1)],maxcoeff=10**6,maxsteps=10**5,tol=mp.mpf(10)**-45)
    print("        [theta,sqrtD,1]:", r2)
sys.stdout.flush()

C={m:consts(m,70) for m in MS}
tol=mp.mpf(10)**-45
one=mp.mpf(1)

trips=[]
for a in range(1,N+1):
    for b in range(a,N+1):
        for c in range(b,N+1):
            if math.gcd(math.gcd(a,b),c)==1:
                trips.append((a,b,c))
print("triples:",len(trips)); sys.stdout.flush()

hits=[]; ctrl_hits=[]
t0=time.time()
for i,(a,b,c) in enumerate(trips):
    if i%2000==0:
        print("  %d/%d  %.0fs"%(i,len(trips),time.time()-t0)); sys.stdout.flush()
    M=CM(a,b,c)
    for m in MS:
        Q,th,A,D=C[m]
        r=mp.pslq([M,Q,th,A,one],maxcoeff=10**6,maxsteps=10**5,tol=tol)
        if r is not None and r[0]!=0:
            hits.append((a,b,c,m,tuple(r)))
        rc=mp.pslq([M,th,A,one],maxcoeff=10**6,maxsteps=10**5,tol=tol)
        if rc is not None and rc[0]!=0:
            ctrl_hits.append((a,b,c,m,tuple(rc)))
print("done %.0fs"%(time.time()-t0))
pickle.dump((hits,ctrl_hits,N),open(WD+'task5b_raw.pkl','wb'))
print("raw hits:",len(hits)," ctrl hits:",len(ctrl_hits))
