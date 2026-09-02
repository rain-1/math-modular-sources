"""High precision r_a = y1/y0 and the asymptotic constant L_a of the minimal solution."""
from mpmath import mp, mpf, factorial, nstr, identify, pslq, sqrt, pi, exp, log, besseli, besselj, gamma, mpc
mp.dps = 220

def cf_value(b, aa, N):
    Am2, Am1 = mpf(1), b(0); Bm2, Bm1 = mpf(0), mpf(1)
    for n in range(1, N+1):
        An = b(n)*Am1 + aa(n)*Am2; Bn = b(n)*Bm1 + aa(n)*Bm2
        Am2, Am1, Bm2, Bm1 = Am1, An, Bm1, Bn
    return Am1/Bm1
C = cf_value(lambda n: mpf(n), lambda n: mpf(1) if n==1 else mpf((-1)**n), 500)
D = cf_value(lambda n: mpf(n+1), lambda n: mpf((-1)**n), 500)
rp = 1 - 1/C     # r_{+1} = y1/y0 for a=+1
rm = 1 - D       # r_{-1}
print("r(+1) =", nstr(rp, 60))
print("r(-1) =", nstr(rm, 60))

def minimal(a, N, M):
    y = [mpf(0)]*(N+2); y[N] = mpf(1)
    for n in range(N+1, 1, -1):
        y[n-2] = (y[n] - n*y[n-1])/(a*(-1)**n)
    return [v/y[0] for v in y[:M]]

# consistency: minimal solution ratio must equal r
for a, r in ((1, rp), (-1, rm)):
    y = minimal(a, 1400, 5)
    print("a=%+d  y1/y0 - r = %s"%(a, nstr(y[1]-r, 5)))

# asymptotic constant: (n+1)! |y_n| -> L.  Use Richardson in 1/n on each parity class.
def Lconst(a, N=1400, M=420):
    y = minimal(a, N, M)
    seq = [factorial(n+1)*abs(y[n]) for n in range(M)]
    return seq
for a in (1,-1):
    seq = Lconst(a)
    # Richardson extrapolation in 1/n, even n only
    import mpmath
    for par in (0,1):
        idx = [n for n in range(200, 380) if n%2==par]
        vals = [seq[n] for n in idx]
        # Richardson on the sequence u_n = L + c1/n + c2/n^2 + ...  use polynomial extrapolation in h=1/n
        K = 14
        ns = idx[:K]; us = [seq[n] for n in ns]
        # Neville in h = 1/n -> h=0
        T = [[mpf(0)]*K for _ in range(K)]
        for i in range(K): T[i][0] = us[i]
        h = [mpf(1)/n for n in ns]
        for j in range(1, K):
            for i in range(K-j):
                T[i][j] = (T[i+1][j-1]*(0-h[i]) - T[i][j-1]*(0-h[i+j]))/(h[i+j]-h[i])
        print("a=%+d parity %d  L = %s"%(a, par, nstr(T[0][K-1], 30)))
