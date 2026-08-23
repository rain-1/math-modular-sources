"""Bost-Charles numerator  BC(phi) = int int_{T^2} log|phi(z)-phi(w)| dmu dmu,
computed as  BC = int int log|(phi(z)-phi(w))/(z-w)|  (since int int log|z-w| = 0),
which removes the diagonal singularity.  Diagonal terms use log|phi'(z)|.

Checks: for phi(z) = rho z (univalent) BC must be log rho.
"""
import mpmath as mp
import cmath, math

mp.mp.dps = 25


def theta2(z, K=200):
    return 2*z**0.25*sum(z**(n*(n+1)) for n in range(K))


def theta3(z, K=200):
    return 1 + 2*sum(z**(n*n) for n in range(1, K))


def lam(z):
    return (theta2(z)/theta3(z))**4


def dlam(z, h=1e-7):
    return (lam(z*(1+h)) - lam(z*(1-h)))/(2*h*z)


def phiA(zeta):
    """mu_4 template on the x-line: x = i sqrt(lambda(zeta^2)) = 4 i zeta + ...
    The branch is fixed by writing lambda(zeta^2) = 16 zeta^2 * (1 + O(zeta^2))
    and taking the principal square root of the bracket (which is near 1)."""
    if zeta == 0:
        return 0j
    z = zeta*zeta
    return 4j*zeta*cmath.sqrt(lam(z)/(16*z))


def phiB(z):
    """u-line lambda template: u = -lambda(z) = -16 z + ..."""
    return -lam(z)


def hfun(Q, K=300):
    """Gamma_0(2) hauptmodul h = -256 Delta(2tau)/Delta(tau), Q = exp(2 pi i tau)."""
    p = 1+0j
    for n in range(1, K):
        p *= ((1-Q**(2*n))/(1-Q**n))**24
    return -256*Q*p


def phiC(Q):
    """v-line: v = u^2/(u+1) with u = -lambda  =>  v = -h."""
    return -hfun(Q)


def BC(phi, r, N=4096):
    """BC integral of z -> phi(r z) on the unit circle, N-point trapezoid."""
    zs = [cmath.exp(2j*math.pi*j/N) for j in range(N)]
    fs = [phi(r*z) for z in zs]
    # derivative of  z -> phi(r z)  at each node, by spectral differentiation
    import numpy as np
    F = np.array(fs)
    c = np.fft.fft(F)/N
    k = np.fft.fftfreq(N, d=1.0/N)
    dF = np.fft.ifft(c*1j*0)  # placeholder
    # d/dtheta f(e^{i th}) = sum k i c_k e^{i k th};  f'(z) = (1/(i z)) df/dth
    dth = np.fft.ifft(c*(1j*k))*N
    Z = np.array(zs)
    fp = dth/(1j*Z)                       # derivative of z -> phi(rz)
    tot = 0.0
    Fa = np.array(fs)
    for j in range(N):
        d = Fa - Fa[j]
        dz = Z - Z[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            q = np.abs(d/dz)
        q[j] = abs(fp[j])
        tot += np.sum(np.log(q))
    return tot/N**2


if __name__ == '__main__':
    import numpy as np
    print('sanity: phi(z) = rho z  =>  BC = log rho')
    for rho in (0.5, 3.0, 161.08):
        print(f'   rho={rho:8.3f}   BC={BC(lambda z: rho*z, 1.0, 1024):+.8f}'
              f'   log rho={math.log(rho):+.8f}')
    print()
    print('=== host A: phi = mu_4 template on |zeta| < r ===')
    print('  (fold preimages of x=+-1 at |zeta| = 0.20788 (one each) and 0.73040)')
    for r in (0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.73):
        b = BC(phiA, r, 2048)
        lp = math.log(4*r)
        print(f'   r={r:5.3f}  log|phi\'(0)|={lp:+.5f}  BC={b:+.5f}  BC-log|phi\'(0)|={b-lp:+.5f}')
    print()
    print('=== host B: phi = -lambda on |z| < r  (u-line) ===')
    print('  (fold preimages of u=1 at |z| = 0.043214 (one) and 0.533488)')
    for r in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.53):
        b = BC(phiB, r, 2048)
        lp = math.log(16*r)
        print(f'   r={r:5.3f}  log|phi\'(0)|={lp:+.5f}  BC={b:+.5f}  BC-log|phi\'(0)|={b-lp:+.5f}')
    print()
    print('=== host C: phi = -h on |Q| < r  (v-line, CDT symmetrisation) ===')
    print('  (fold preimages of v=1/2 at |Q| = 0.0018674 (one) and 0.284648)')
    for r in (0.01, 0.05, 0.1, 0.2, 0.28):
        b = BC(phiC, r, 2048)
        lp = math.log(256*r)
        print(f'   r={r:5.3f}  log|phi\'(0)|={lp:+.5f}  BC={b:+.5f}  BC-log|phi\'(0)|={b-lp:+.5f}')
