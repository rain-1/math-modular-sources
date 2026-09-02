from mpmath import mp, mpf, mpc, factorial, binomial, findroot, nstr, besseli
mp.dps = 40
# minimal solution of y_n = n y_{n-1} + s(-1)^n y_{n-2}
def minimal(s, N=300, M=30):
    y = [mpf(0)]*(N+2); y[N]=mpf(1)
    for n in range(N+1,1,-1):
        y[n-2] = (y[n]-n*y[n-1])/(s*(-1)**n)
    return [v/y[0] for v in y[:M]]
# Taylor coefficients of f(t) = phi(1+t) psi(1-t), phi(u)=sum a^j u^j/(j!)^2, psi(v)=sum b^k v^k/(k!)^2 (I0(2 sqrt(a u)) I0(2 sqrt(b v)))
def coeffs_product(a, b, nmax=8, J=60):
    c = [mpc(0)]*(nmax+1)
    for j in range(J):
        for k in range(J):
            w = a**j * b**k / (factorial(j)**2 * factorial(k)**2)
            if abs(w) < mpf(10)**(-mp.dps-5) and j+k>10: continue
            # [t^n] (1+t)^j (1-t)^k
            for n in range(nmax+1):
                s = mpc(0)
                for i in range(max(0,n-k), min(j,n)+1):
                    s += binomial(j,i)*binomial(k,n-i)*(-1)**(n-i)
                c[n] += w*s
    return c
for s in (1,-1):
    y = minimal(s)
    target = [y[n]/factorial(n) for n in range(9)]  # [t^n] f = y_n/n!
    def resid(a,b):
        c = coeffs_product(a,b,nmax=2)
        return [c[1]/c[0]-target[1]/target[0], c[2]/c[0]-target[2]/target[0]]
    for guess in [(mpc(1,0),mpc(1,0)), (mpc(0.5,0.5),mpc(0.5,-0.5)), (mpc(0,1),mpc(0,-1)), (mpc(-1,0),mpc(1,0)), (mpc(0.25,0),mpc(-0.25,0))]:
        try:
            sol = findroot(lambda a,b: resid(a,b), guess, tol=1e-25, maxsteps=60)
            a,b = sol[0], sol[1]
            c = coeffs_product(a,b,nmax=8)
            errs = [abs(c[n]/c[0]-target[n]/target[0]) for n in range(3,9)]
            print("s=%d guess=%s -> a=%s b=%s ; check errs n=3..8:"%(s,guess,nstr(a,15),nstr(b,15)), [nstr(e,3) for e in errs])
        except Exception as e:
            print("s=%d guess=%s failed: %s"%(s,guess,e))
