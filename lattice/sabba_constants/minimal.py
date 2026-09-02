from mpmath import mp, mpf, factorial, identify, pslq, nstr, besseli, hyper
mp.dps = 60
# recurrence y_n = n*y_{n-1} + eps_n*y_{n-2}, eps_n = s*(-1)^n. Minimal solution by backward recurrence (Miller):
# y_{n-2} = (y_n - n*y_{n-1})/eps_n
def minimal(s, N=400, M=60):
    y = [mpf(0)]*(N+2)
    y[N+1] = mpf(0); y[N] = mpf(1)
    for n in range(N+1, 1, -1):
        eps = s*(-1)**n
        y[n-2] = (y[n] - n*y[n-1])/eps
    return y[:M]
for s in (+1, -1):
    y = minimal(s)
    y0 = y[0]
    print("s =", s, " twist eps_n = s*(-1)^n")
    # normalise so that y_0 = 1 and print n! * y_n and (n+1)! y_n
    for n in range(0, 12):
        print(n, nstr(y[n]/y0, 25), " n!y_n/y0 =", nstr(factorial(n)*y[n]/y0, 25), " (n+1)!y_n/y0 =", nstr(factorial(n+1)*y[n]/y0, 25))
    # ratios giving the tails: T(n) = y_{n-1}... check CF values: C = 1/T_+(1), D = T_-(1)
    # Pincherle: for y_n = n y_{n-1} + eps_n y_{n-2}, minimal y: y_{n-1}/y_{n-2} = n... let's just print y_n/y_{n-1}
    for n in range(1, 5):
        print("  y_%d/y_%d ="%(n,n-1), nstr(y[n]/y[n-1], 30))
