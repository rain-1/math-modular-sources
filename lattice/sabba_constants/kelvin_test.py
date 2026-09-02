from mpmath import mp, mpf, besseli, besselj, exp, pi, mpc, pslq, ber, bei, sqrt
mp.dps = 80
# Sabba's constants via convergents. CF: b0 + a1/(b1 + a2/(b2 + ...))
def cf_value(b, a, N):
    Am2, Am1 = mpf(1), b(0)   # A_{-1}, A_0
    Bm2, Bm1 = mpf(0), mpf(1)
    for n in range(1, N+1):
        An = b(n)*Am1 + a(n)*Am2
        Bn = b(n)*Bm1 + a(n)*Bm2
        Am2, Am1, Bm2, Bm1 = Am1, An, Bm1, Bn
    return Am1/Bm1
# C = 1/(1 + 1/(2 - 1/(3 + 1/(4 - ...)))): b0=0,a1=1,b_n=n (n>=1), a_n=(-1)^n (n>=2)
C = cf_value(lambda n: mpf(n), lambda n: mpf(1) if n==1 else mpf((-1)**n), 120)
# D = 1 - 1/(2 + 1/(3 - 1/(4 + ...))): b_n = n+1, a_n = (-1)^n
D = cf_value(lambda n: mpf(n+1), lambda n: mpf((-1)**n), 120)
# Rabinowitz R = 1 + 1/(2 + 1/(3+...)) : b_n=n+1, a_n=1 ; Fourier Fo = 1 - 1/(2 - 1/(3 - ...)): a_n=-1
R = cf_value(lambda n: mpf(n+1), lambda n: mpf(1), 120)
Fo = cf_value(lambda n: mpf(n+1), lambda n: mpf(-1), 120)
print("C =", C); print("D =", D)
print("R - I0(2)/I1(2) =", R - besseli(0,2)/besseli(1,2))
print("Fo - J0(2)/J1(2) =", Fo - besselj(0,2)/besselj(1,2))
# Kelvin guesses
w = exp(mpc(0,1)*pi/4)
I0w, I1w = besseli(0, 2*w), besseli(1, 2*w)
basis_names = ["Re I0", "Im I0", "Re I1", "Im I1"]
basis = [I0w.real, I0w.imag, I1w.real, I1w.imag]
for name, X in [("C", C), ("D", D)]:
    vec = [X*b for b in basis] + basis
    rel = pslq(vec, maxcoeff=10**6, maxsteps=10**5)
    print(name, "PSLQ over Q:", rel)
    # also allow sqrt2
    vec2 = vec + [sqrt(2)*v for v in vec]
    rel2 = pslq(vec2, maxcoeff=10**5, maxsteps=10**5)
    print(name, "PSLQ over Q(sqrt2):", rel2)
print("ber0(2),bei0(2),ber1(2),bei1(2) =", ber(0,2), bei(0,2), ber(1,2), bei(1,2))
