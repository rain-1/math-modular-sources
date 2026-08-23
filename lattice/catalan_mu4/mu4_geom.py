"""Conformal geometry of the mu_4 hosts.

Three nested hosts (base point x = 0 in every case):

  A  x-line :  phi : D -> C \ {i,-i},   phi(0)=0,  phi^{-1}(0)={0}.
               Ceiling = |varphi_4'(0)| for the "mu_4 template"
               varphi_4(zeta) = i*sqrt(lambda(zeta^2)),  the quotient of the
               universal cover of P^1-{0,i,-i,oo} by the parabolic at x=0.
  B  u-line (u=x^2) :  phi : D -> C \ {-1},  phi^{-1}(0)={0}.  Template
               varphi(z) = -lambda(z),  |varphi'(0)| = 16   (CDT's lambda template).
  C  v-line (v=u^2/(u+1)) : the CDT involution quotient,  |varphi'(0)| = 256.

Also: the Schwarz-Pick ceilings (no condition on phi^{-1}(0)) and the
positions of the fold preimages.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'cdt_finder'))
import mpmath as mp
from conformal import r_two_punctures, _lam, _lam_prime, _invlambda

mp.mp.dps = 30


def lam_of_z(z):
    """lambda as a function of the nome z = exp(i pi tau)."""
    t2 = 2*z**mp.mpf('0.25')*sum(z**(n*(n+1)) for n in range(0, 40))
    t3 = 1 + 2*sum(z**(n*n) for n in range(1, 40))
    return (t2/t3)**4


def phi4(zeta, nterms=60):
    """the mu_4 template  x = i sqrt(lambda(zeta^2)),  a power series in zeta."""
    z = zeta**2
    return 1j*mp.sqrt(lam_of_z(z))


if __name__ == '__main__':
    print('=== ceilings ===')
    # (A) exact:  lambda(z) = 16 z + O(z^2)  =>  sqrt(lambda(zeta^2)) = 4 zeta + ...
    for zt in ('1e-4', '1e-5', '1e-6'):
        zt = mp.mpf(zt)
        print(f'   |phi4(zeta)/zeta| at zeta={mp.nstr(zt,3):>7s} : '
              f'{mp.nstr(abs(phi4(zt))/zt, 15)}   (-> 4 exactly)')
    print('   A  x-line   ceiling |phi\'(0)| = 4        log =', mp.nstr(mp.log(4), 10))
    print('   B  u-line   ceiling |phi\'(0)| = 16       log =', mp.nstr(mp.log(16), 10))
    print('   C  v-line   ceiling |phi\'(0)| = 256      log =', mp.nstr(mp.log(256), 10))
    print()
    print('=== Schwarz-Pick ceilings (no phi^{-1}(0)={0} condition) ===')
    r = r_two_punctures(1j, -1j, 0)
    print('   C\\{i,-i} at 0                :', mp.nstr(r, 15), '  log =', mp.nstr(mp.log(r), 10))
    for pair in [(1, -1), (1, 1j), (1j, -1)]:
        rr = r_two_punctures(pair[0], pair[1], 0)
        print(f'   C\\{{{pair[0]},{pair[1]}}} at 0'.ljust(33) + ':', mp.nstr(rr, 15))
    print('   => the 6-point host P^1-{0,+-1,+-i,oo} (phi must avoid all four of')
    print('      +-1,+-i) has Schwarz-Pick ceiling <= min over pairs =',
          mp.nstr(min(r_two_punctures(a, b, 0)
                      for a, b in [(1, -1), (1j, -1j), (1, 1j), (1, -1j), (-1, 1j), (-1, -1j)]), 12))
    print()
    print('=== fold preimages in the mu_4 template (host A) ===')
    print('   fold points x = +-1  <=>  u = 1  <=>  lambda = -1  <=>  tau in Gamma(2).(i-1)')
    print('   |zeta| = |z|^{1/2} = exp(-pi*Im(tau)/2):')
    # Gamma(2) orbit of tau0 = -1+i : Im = 1/((d-c)^2+c^2), c even, d odd
    seen = []
    for c in range(-8, 9, 2):
        for d in range(-9, 10, 2):
            if (c, d) == (0, 0):
                continue
            n = (d-c)**2 + c**2
            if n == 0:
                continue
            seen.append((1/mp.mpf(n), c, d))
    seen.sort(reverse=True)
    shown = set()
    for im, c, d in seen:
        key = mp.nstr(im, 8)
        if key in shown:
            continue
        shown.add(key)
        print(f'     Im(tau)=1/{int(round(1/im)):<4d}  |z|=exp(-pi*Im)={mp.nstr(mp.e**(-mp.pi*im),8):>12s}'
              f'   |zeta|={mp.nstr(mp.e**(-mp.pi*im/2),8):>12s}')
        if len(shown) >= 6:
            break
