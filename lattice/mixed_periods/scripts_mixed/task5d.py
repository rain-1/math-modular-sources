import mpmath as mp, math, time, sys, pickle
WD='/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/9a849c0a-95f8-4d19-b342-98033d0d9c03/scratchpad/mixed/'
def DBW(z):
    z=mp.mpc(z); return mp.im(mp.polylog(2,z))+mp.arg(1-z)*mp.log(abs(z))
mp.mp.dps=70
N=60; MS=[2,3,5]; tol=mp.mpf(10)**-45; one=mp.mpf(1)
C={}
for m in MS:
    D=4*m-1; A=mp.sqrt(D); th=2*mp.atan(1/A)
    Q=2*mp.clsin(2,mp.pi-th)-th*mp.log(mp.mpf(D)/m)
    C[m]=(Q,th,A,D,mp.clsin(2,mp.pi-th),mp.pi*mp.log(mp.mpf(D)/m))
trips=[(a,b,c) for a in range(1,N+1) for b in range(a,N+1) for c in range(b,N+1)
       if math.gcd(math.gcd(a,b),c)==1]
cl_hits=[]; aug_hits=[]
t0=time.time()
for i,(a,b,c) in enumerate(trips):
    if i%3000==0: print(i,len(trips),"%.0fs"%(time.time()-t0)); sys.stdout.flush()
    if not (c<a+b): continue
    A_=mp.mpf(a); B_=mp.mpf(b); Cc=mp.mpf(c)
    al=mp.acos((B_*B_+Cc*Cc-A_*A_)/(2*B_*Cc)); be=mp.acos((A_*A_+Cc*Cc-B_*B_)/(2*A_*Cc)); ga=mp.acos((A_*A_+B_*B_-Cc*Cc)/(2*A_*B_))
    Bl=DBW((B_/A_)*mp.exp(1j*ga))
    M=al*mp.log(A_)+be*mp.log(B_)+ga*mp.log(Cc)+Bl
    for m in MS:
        Q,th,Asq,D,K,pl=C[m]
        r=mp.pslq([Bl,K],maxcoeff=10**6,maxsteps=10**5,tol=tol)
        if r is not None and r[0]!=0 and r[1]!=0: cl_hits.append((a,b,c,m,tuple(r)))
        r2=mp.pslq([M,Q,th,Asq,one,pl],maxcoeff=10**4,maxsteps=10**5,tol=tol)
        if r2 is not None and r2[0]!=0 and r2[1]!=0: aug_hits.append((a,b,c,m,tuple(r2)))
print("done %.0fs"%(time.time()-t0))
pickle.dump((cl_hits,aug_hits),open(WD+'task5d_raw.pkl','wb'))
print("clausen-part hits:",len(cl_hits))
for h in cl_hits[:80]: print("  ",h)
print("augmented-basis hits:",len(aug_hits))
for h in aug_hits[:80]: print("  ",h)
