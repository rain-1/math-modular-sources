"""Beyond Koebe: the Bost-Charles numerator of the multivalent (Andre/Kodaira)
modular map, computed numerically.

For a second-order row the module {1, H, theta H, ...} is holomorphic at x = 0
and has its only finite branch point at x_2 = 1/lam_2.  Any phi : D -> C \ {x_2}
with phi(0) = 0 is admissible (D simply connected => phi lifts; the f_i are
holomorphic at 0).  Two natural families:

  * UNIVALENT (Koebe):  phi = Riemann map of C \ [x_2, oo);  |phi'(0)| = 4|x_2|,
    and BC = log|phi'(0)| exactly (Grunsky).  Deficit Delta := BC - log|phi'(0)| = 0.
  * MULTIVALENT (Kodaira):  phi_r(z) = x_2 * lambda(r z), lambda = 16w -128w^2+...
    the modular lambda function in the nome w = e^{i pi tau}.  lambda(D) = C\{0,1}
    and loops in D map only to loops around 0, so phi_r^* f is single-valued.
    |phi_r'(0)| = 16 |x_2| r  ->  16|x_2| as r -> 1 (the Kodaira ceiling), but
    Delta(r) = BC(phi_r) - log|phi_r'(0)| grows.

BC is computed as   BC = log|phi'(0)| + iint log|(phi(z)-phi(w))/((z-w) phi'(0))|
using  iint_{T^2} log|z-w| = 0  exactly; the diagonal value is log|phi'(z)/phi'(0)|.
Trapezoidal quadrature in both angles (spectrally accurate for the smooth part).

Margin of the CDT bound at m functions:  m(log|phi'(0)| - tau) - BC
                                       = (m-1) log|phi'(0)| - m tau - Delta.
So the multivalent map is worth  (m-1) log(4r) - Delta(r)  over Koebe.
"""
import math, cmath, os, sys
import mpmath as mp

mp.mp.dps = 25

def lam(w):
    return (mp.jtheta(2, 0, w)/mp.jtheta(3, 0, w))**4

def lam_prime(w):
    t3 = mp.jtheta(3, 0, w)
    l = lam(w)
    # d lambda/d tau = i pi theta3^4 lambda (1-lambda);  w = e^{i pi tau}
    # => d lambda/dw = (d lambda/d tau)/(i pi w) = theta3^4 lambda(1-lambda)/w
    return t3**4*l*(1-l)/w

def delta_of_r(r, N=512):
    """Delta = BC(phi_r) - log|phi_r'(0)| for phi_r(z) = lambda(r z)  (scale-free)."""
    r = mp.mpf(r)
    zs = [mp.e**(2j*mp.pi*mp.mpf(k)/N) for k in range(N)]
    f  = [lam(r*z) for z in zs]
    fp = [lam_prime(r*z)*r for z in zs]
    d0 = mp.mpf(16)*r                      # |phi_r'(0)|
    tot = mp.mpf(0)
    for i in range(N):
        for j in range(N):
            if i == j:
                v = mp.log(abs(fp[i]/d0))
            else:
                v = mp.log(abs((f[i]-f[j])/((zs[i]-zs[j])*d0)))
            tot += v
    return float(tot/(N*N))

if __name__ == '__main__':
    print(__doc__)
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 384
    print(f"quadrature N = {N} per circle\n")
    print(f"{'r':>7} {'log|phi_r\'(0)|-log|x_2|':>24} {'Delta(r)':>10} "
          f"{'gain m=2':>10} {'gain m=3':>10} {'gain m=4':>10} {'gain m=14':>10}")
    print("  (gain = (m-1)*log(4r) - Delta(r) : how much phi_r beats the Koebe map at m functions)")
    best = {}
    for r in (0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75,0.80):
        D = delta_of_r(r, N)
        L = math.log(16*r)
        row = [( m, (m-1)*math.log(4*r) - D) for m in (2,3,4,14)]
        print(f"{r:>7.2f} {L:>24.5f} {D:>10.5f} " + " ".join(f"{g:>10.5f}" for _, g in row))
        for m, g in row:
            if m not in best or g > best[m][1]: best[m] = (r, g, D, L)
    print()
    for m in (2,3,4,14):
        r, g, D, L = best[m]
        print(f"  best for m={m:>2}: r={r:.2f}, log|phi'(0)| = log|x_2| + {L:.5f}, "
              f"Delta={D:.5f}, gain over Koebe = {g:+.5f}")
