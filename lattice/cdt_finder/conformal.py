"""Poincare (conformal) radius |phi'(0)| of punctured planes.

phi ranges over ALL holomorphic maps D -> Omega with phi(0)=t0 (CDT's "common
pullback" need not be injective), so the supremum of |phi'(0)| is 1/lambda_Omega(t0),
the reciprocal of the hyperbolic density (curvature -4 normalisation, so that
r_D(0)=1).  Monotone: Omega subset Omega' => r_Omega <= r_Omega'.

Exact cases implemented:
  * Omega = C \ {a,b}          (thrice-punctured sphere)   -- via the modular
    lambda function; r = |b-a| * |lambda'(tau0)| * 2 Im(tau0).
  * Omega = C \ [a,oo)  (one radial cut, Koebe)             -- r = 4|a|.
For >=3 finite punctures we return the rigorous upper bound min over pairs
(Schwarz-Pick + monotonicity) and the lower bound  max over discs / one-cut.
"""
import mpmath as mp

mp.mp.dps = 40


def _lam(tau):
    q = mp.e**(1j*mp.pi*tau)
    t2 = mp.jtheta(2, 0, q); t3 = mp.jtheta(3, 0, q)
    return (t2/t3)**4


def _lam_prime(tau):
    q = mp.e**(1j*mp.pi*tau)
    t3 = mp.jtheta(3, 0, q)
    l = _lam(tau)
    return 1j*mp.pi*t3**4*l*(1-l)


def _invlambda(z):
    """some tau in H with lambda(tau) = z, via tau = i K(1-z)/K(z)."""
    tau = 1j*mp.ellipk(1-z)/mp.ellipk(z)
    if mp.im(tau) <= 0:
        tau = -mp.conj(tau)
    # Newton polish
    for _ in range(60):
        d = (_lam(tau)-z)/_lam_prime(tau)
        tau = tau - d
        if abs(d) < mp.mpf(10)**(-mp.mp.dps+5):
            break
    return tau


def r_two_punctures(a, b, t0=0):
    """conformal radius at t0 of C \\ {a,b}."""
    a = mp.mpmathify(a); b = mp.mpmathify(b); t0 = mp.mpmathify(t0)
    z0 = (t0-a)/(b-a)                     # 0->? ; a->0, b->1
    tau0 = _invlambda(z0)
    assert abs(_lam(tau0)-z0) < mp.mpf(10)**(-20), (a, b, abs(_lam(tau0)-z0))
    return abs(b-a)*abs(_lam_prime(tau0))*2*mp.im(tau0)


def r_upper_multi(sing, t0=0):
    """rigorous upper bound for C \\ sing (finite set, >=2 points), at t0."""
    sing = [mp.mpmathify(s) for s in sing]
    best = None
    for i in range(len(sing)):
        for j in range(i+1, len(sing)):
            v = r_two_punctures(sing[i], sing[j], t0)
            if best is None or v < best:
                best = v
    return best


def r_lower_onecut(sing, t0=0):
    """If every puncture lies on one ray from t0, Omega contains the cut plane;
    r >= 4*dist.  Otherwise fall back to the inscribed disc  r >= dist."""
    sing = [mp.mpmathify(s) for s in sing]
    d = min(abs(s-t0) for s in sing)
    args = [mp.arg(s-t0) for s in sing]
    collinear = all(abs(mp.e**(1j*(x-args[0]))-1) < mp.mpf('1e-25') for x in args)
    return (4*d, 'one-cut (all singularities on a common ray)') if collinear else (d, 'inscribed disc')


if __name__ == '__main__':
    # unit test 1: Koebe.  C\{1/4} is not hyperbolic; use C\{1/4, R} with R->oo
    # and compare with the cut plane 4*(1/4)=1: two punctures give MORE room.
    print('C\\{-1,1} at 0 :', r_two_punctures(-1, 1))
    # unit test 2: reproduce the archive's level-10 number.
    # they normalised C\{u+,u-} to C\{0,1} with base z0 = 1/2 - i, got r_u=1.30651071...
    z0 = mp.mpf('0.5')-1j
    tau0 = _invlambda(z0)
    print('tau0 =', mp.nstr(tau0, 22), '  (archive: -0.53238870386173469076+0.84650001063226346634i)')
    r01 = abs(_lam_prime(tau0))*2*mp.im(tau0)
    print('r for C\\{0,1} at z0 =', mp.nstr(r01, 22))
