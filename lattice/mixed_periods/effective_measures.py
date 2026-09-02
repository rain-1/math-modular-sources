"""Effective irrationality exponents for the fold periods, from CDT's quantitative bound
(ICM survey, Prop. `theomainreprise`): m <= BC/(log|phi'(0)| - tau - (1-gamma)(2/mu - (1-gamma)/mu^2) log(1/rho)).
Bivalent map phi(z)=8(z+z^3)/(1+z)^4, BC = log 8 + 4G/pi. Inventory for eta=c_D: {1, log(1-x), log^2(1-x), g, g o iota},
types (0,1,3/2,3/2,3/2), gamma=3/5; for eta=c_B: {1, log(1-x), g, g o iota}, types (0,1,1,1), gamma=1/2.
rho = radius with max_{|z|=rho}|phi| = 1/N (convergence radius of H, H_B, H_D). See MIXED_PERIODS_HYPERGEOMETRIC.md s.10."""
from mpmath import mp, mpf, log, findroot, catalan, pi, nstr, sqrt
mp.dps = 30
BC = log(8) + 4*catalan/pi; l8 = log(8)
def rho_of(delta): return findroot(lambda r: 8*r*(1+r*r)/(1-r)**4 - delta, delta/8)
def kappa(m, gnum, sig_cond, sig_pure, L):
    g = mpf(gnum)/m
    rows = sorted(sig_pure + [sig_cond]*(m-gnum))
    tau = sum((2*(i+1)-1)*rows[i] for i in range(m))/mpf(m)**2
    gap = l8 - tau - BC/m
    if gap <= 0: return None
    if gap > L: return mpf('inf')
    return (1-g)/(1-sqrt(1-gap/L))   # solves (1-g)(2/mu-(1-g)/mu^2) L = gap
if __name__ == "__main__":
    for N in [4, 8, 12, 16, 20, 44, 164, 9, 27, 80, 25]:
        r = rho_of(mpf(1)/N); L = log(1/r)
        kD = kappa(5, 3, mpf(3)/2, [mpf(0), mpf(1), mpf(3)/2], L)
        kB = kappa(4, 2, mpf(1), [mpf(0), mpf(1)], L)
        worst = max(kappa(mm, 3, mpf(3)/2, [mpf(0), mpf(1), mpf(3)/2], L) for mm in range(5, 40))
        print("N=%4d  rho=%s  L=%s  kappa(c_D)=%s (max over m>=5: %s)  kappa(c_B)=%s" % (N, nstr(r,6), nstr(L,6), nstr(kD,6), nstr(worst,6), nstr(kB,6)))
