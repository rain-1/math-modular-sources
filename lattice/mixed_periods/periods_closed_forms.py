"""Closed forms of the fold periods on the hosts 1/sqrt(1-4mx) (imaginary family) and
1/sqrt(1+4mx) (real family); PSLQ checks. See consolidation/MIXED_PERIODS_HYPERGEOMETRIC.md."""
from mpmath import mp, mpf, sqrt, log, atan, pi, quad, clsin, pslq, polylog, nstr
mp.dps = 50
print("imaginary family: c_B = theta/sqrt(D'),  c_D = (theta log(D'/m) - 2 Cl2(pi-theta))/sqrt(D')")
for m in [1, 2, 3, 5, 11]:
    D = 4*m-1; a = sqrt(D); th = 2*atan(1/a)
    cB = 2*quad(lambda u: 1/(D+u*u), [0, 1])
    cD = 2*quad(lambda u: log((D+u*u)/(4*m))/(D+u*u), [0, 1])
    print(" m=%2d  |cB-formula|=%s  |cD-formula|=%s" % (m, nstr(abs(cB-th/a), 5), nstr(abs(cD-(th*log(mpf(D)/m)-2*clsin(2, pi-th))/a), 5)))
print("real family: c_B = log((c+1)/(c-1))/c,  c_D = R_m/c")
for m in [1, 2, 3, 6, 12]:
    A2 = 4*m+1; c = sqrt(A2)
    cB = 2*quad(lambda u: 1/(A2-u*u), [0, 1])
    cD = 2*quad(lambda u: log((A2-u*u)/(4*m))/(A2-u*u), [0, 1])
    eps = (c+1)/(c-1)
    R = 2*polylog(2, (c-1)/(2*c)) - pi**2/6 + log(2*c)*log(eps) + log((c+1)/(2*c))*log((c-1)/(2*c)) - log(c+1)**2/2 + log(c-1)**2/2
    print(" m=%2d c=%s |cB-formula|=%s  |cD-formula|=%s" % (m, nstr(c, 8), nstr(abs(cB-log(eps)/c), 5), nstr(abs(cD-R/c), 5)))
