"""Basic data: C, D to high precision; minimal solutions of y_n = n y_{n-1} + a(-1)^n y_{n-2}."""
from mpmath import mp, mpf, factorial, nstr, identify, pslq, sqrt, pi, gamma, exp, log
mp.dps = 140

def cf_value(b, a, N):
    Am2, Am1 = mpf(1), b(0)
    Bm2, Bm1 = mpf(0), mpf(1)
    for n in range(1, N+1):
        An = b(n)*Am1 + a(n)*Am2
        Bn = b(n)*Bm1 + a(n)*Bm2
        Am2, Am1, Bm2, Bm1 = Am1, An, Bm1, Bn
    return Am1/Bm1

C = cf_value(lambda n: mpf(n), lambda n: mpf(1) if n==1 else mpf((-1)**n), 300)
D = cf_value(lambda n: mpf(n+1), lambda n: mpf((-1)**n), 300)
Cref = mpf("0.62873660709895479943558790225263029538343274690833134560117673060807682799050577")
Dref = mpf("0.57663338973018439239789174978291392579614943527571083984110419180764835440124402")
print("C  =", nstr(C, 100))
print("C-Cref =", nstr(C-Cref, 5))
print("D  =", nstr(D, 100))
print("D-Dref =", nstr(D-Dref, 5))

def minimal(a, N=900, M=40):
    """minimal solution of y_n = n y_{n-1} + a(-1)^n y_{n-2}, normalised y_0 = 1."""
    y = [mpf(0)]*(N+2); y[N] = mpf(1)
    for n in range(N+1, 1, -1):
        y[n-2] = (y[n] - n*y[n-1])/(a*(-1)**n)
    y0 = y[0]
    return [v/y0 for v in y[:M]]

for a in (1, -1):
    y = minimal(a)
    print()
    print("=== a = %d ==="%a)
    for n in range(0, 10):
        print("  y_%-2d = %-45s  (n+1)! y_n = %s"%(n, nstr(y[n], 35), nstr(factorial(n+1)*y[n], 30)))
    print("  y_1/y_0 =", nstr(y[1]/y[0], 40))
    print("  y_0/y_1 =", nstr(y[0]/y[1], 40))
    print("  1/(1+y_1/y_0)?", nstr(1/(1+y[1]/y[0]),40))
    print("  signs n=0..24:", [1 if y[n]>0 else -1 for n in range(25)])
    print("  (-1)^(n(n+1)/2) :", [(-1)**((n*(n+1))//2) for n in range(25)])
