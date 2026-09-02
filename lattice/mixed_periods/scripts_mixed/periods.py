from mpmath import mp, mpf, sqrt, log, atan, pi, quad, clsin, pslq, polylog, zeta, dirichlet, identify, mpc, im
mp.dps = 60
def cl2(t): return clsin(2,t)
for m in [1,2,3,4,5,6,11,17,41]:
    D = 4*m-1; a = sqrt(D); th = 2*atan(1/a)
    # integrals via u-substitution: c_B = 2 int_0^1 du/(D+u^2); c_D = 2 int_0^1 log((D+u^2)/(4m))/(D+u^2) du
    cB = 2*quad(lambda u: 1/(D+u*u), [0,1])
    cD = 2*quad(lambda u: log((D+u*u)/(4*m))/(D+u*u), [0,1])
    # c_C = int_0^{1/4m} log(1-t)/(t sqrt(1-4mt)) dt : t=(1-u^2)/(4m), dt=-u du/(2m), 1/t = 4m/(1-u^2): = int_0^1 log((D+u^2)/(4m)) * 4m/(1-u^2) * 1/u * u/(2m) du = 2 int_0^1 log((D+u^2)/(4m))/(1-u^2) du
    cC = 2*quad(lambda u: log((D+u*u)/(4*m))/(1-u*u), [0,1])
    print("m=%d D=%d  theta=%s"%(m,D,mp.nstr(th,20)))
    print("  cB=",mp.nstr(cB,25), " check (2/a)atan(1/a)=",mp.nstr(2*atan(1/a)/a,25))
    basis_names = ["cD","th/a","th*log2/a","th*logD/a","th*logm/a","Cl2(th)/a","Cl2(pi-th)/a","pi^2/a","pi*th/a","1"]
    basis = [cD, th/a, th*log(2)/a, th*log(D)/a, (th*log(m)/a if m>1 else mpf(0)), cl2(th)/a, cl2(pi-th)/a, pi**2/a, pi*th/a, mpf(1)]
    if m==1: basis[4]=mpf(0)
    idx=[i for i,b in enumerate(basis) if b!=0]
    r = pslq([basis[i] for i in idx], maxcoeff=10**6, maxsteps=10**6)
    print("  PSLQ cD:", None if r is None else {basis_names[idx[i]]:r[i] for i in range(len(r)) if r[i]!=0})
    basisC = [cC, pi**2, th**2, log(2)**2, log(D)**2, log(2)*log(D), (log(m)**2 if m>1 else mpf(0)), cl2(th), cl2(pi-th), th*log(2), th*log(D), mpf(1)]
    namesC = ["cC","pi^2","th^2","log2^2","logD^2","log2logD","logm^2","Cl2(th)","Cl2(pi-th)","th log2","th logD","1"]
    idx=[i for i,b in enumerate(basisC) if b!=0]
    r = pslq([basisC[i] for i in idx], maxcoeff=10**6, maxsteps=10**6)
    print("  PSLQ cC:", None if r is None else {namesC[idx[i]]:r[i] for i in range(len(r)) if r[i]!=0})
    # L-value test: is Cl2(pi-th) or Cl2(th) a rational multiple of sqrt(D) L(2,chi_{-D}) (D prime => fundamental disc -D when D=3 mod 4)?
    if D in (3,7,11,19,23,31,43,67,163):
        L2 = dirichlet(2, [0]+[(-1 if pow(k,(D-1)//2,D)==D-1 else 1) if k%D else 0 for k in range(1,D)]) if False else None
        # build chi_{-D}(k) = kronecker(-D,k)
        from mpmath import kronecker  # not available? fallback
