"""Final ledger for the mu_4 polylogarithm architectures for Catalan's constant.

For each host: inventory (measured denominator types), tau(b;e), the admissible
concentric contour radius r (fixed by the fold-preimage geometry), the exactly
computed Bost-Charles numerator, entry and signed margin.

    margin := m ( log|phi'(0)| - tau )  -  BC(phi)      ( > 0  ==  contradiction )
"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cdt_finder'))
from fractions import Fraction as F
from mu4_tau import tau_of
from mu4_bc import BC, phiA, phiB, phiC

half = F(1, 2)

HOSTS = {
 'A  x-line  P^1-{0,i,-i,oo}': dict(
    scale=4.0, phi=phiA, rmax=0.73040, folds='x=+-1 at |zeta|=0.20788; next 0.73040',
    pure=[('1', 0, []), ('arctan x', 1, []), ('log(1+x^2)', 1, []),
          ('Li_2(-x^2)', 2, []), ('f(x)=Ti_2(x)', 2, []),
          ('log^2(1+x^2)', 1, [half]), ('M(x)', 1, [half]),
          ('arctan^2 x', 1, [1]), ('arctan x*log(1+x^2)', 1, [1]),
          ('J(-x^2)', 0, [half, half])],
    pure_relax=[('1', 0, []), ('arctan x', 0, [1]), ('log(1+x^2)', 0, [1]),
          ('Li_2(-x^2)', 0, [1, 1]), ('f(x)', 0, [1, 1]),
          ('log^2(1+x^2)', 0, [1, half]), ('M(x)', 0, [1, half]),
          ('arctan^2 x', 0, [1, 1]), ('arctan x*log', 0, [1, 1]),
          ('J(-x^2)', 0, [half, half])],
    cond=(1, [1, 1]), cond_relax=(0, [1, 1])),
 'B  u-line  P^1-{0,-1,oo}': dict(
    scale=16.0, phi=phiB, rmax=0.53349, folds='u=1 at |z|=0.043214; next 0.533488',
    pure=[('1', 0, []), ('ell(u)', 0, [2]), ('g(u)', 0, [2, 2]),
          ('log(1+u)', 1, []), ('Li_2(-u)', 2, []), ('log^2(1+u)', 1, [1]),
          ('u ell^2', 1, [2]), ('ell log(1+u)', 0, [1, 2]),
          ('gtilde(u)', 0, [1, 2]), ('J(-u)', 0, [1, 1])],
    pure_relax=[('1', 0, []), ('ell(u)', 0, [2]), ('g(u)', 0, [2, 2]),
          ('log(1+u)', 0, [1]), ('Li_2(-u)', 0, [1, 1]), ('log^2(1+u)', 0, [1, 1]),
          ('u ell^2', 0, [2, 1]), ('ell log(1+u)', 0, [1, 2]),
          ('gtilde(u)', 0, [1, 2]), ('J(-u)', 0, [1, 1])],
    cond=(1, [2, 2]), cond_relax=(0, [2, 2])),
 'C  v-line  (CDT involution quotient of B)': dict(
    scale=256.0, phi=phiC, rmax=0.284648, folds='v=1/2 at |Q|=0.0018674; next 0.284648',
    pure=[('1', 0, []), ('Sym ell', 0, [4]), ('Sym g', 0, [4, 4]),
          ('Sym log', 1, []), ('Sym Li_2', 2, []), ('Sym log^2', 1, [2]),
          ('Sym u ell^2', 1, [4]), ('Sym ell log', 0, [2, 4]),
          ('Sym gtilde', 0, [2, 4]), ('Sym J', 0, [2, 2])],
    pure_relax=[('1', 0, []), ('Sym ell', 0, [4]), ('Sym g', 0, [4, 4]),
          ('Sym log', 0, [2]), ('Sym Li_2', 0, [2, 2]), ('Sym log^2', 0, [2, 2]),
          ('Sym u ell^2', 0, [4, 2]), ('Sym ell log', 0, [2, 4]),
          ('Sym gtilde', 0, [2, 4]), ('Sym J', 0, [2, 2])],
    cond=(1, [4, 4]), cond_relax=(0, [4, 4])),
}


def run(hostname, H, c, relax, rs):
    pure = H['pure_relax'] if relax else H['pure']
    ce, cl = H['cond_relax'] if relax else H['cond']
    inv = list(pure) + [(f'H_{j}', ce, cl) for j in range(c)]
    d = tau_of(inv)
    m, tau = d['m'], d['tau']
    best = None
    for r in rs:
        lp = math.log(H['scale']*r)
        bc = BC(H['phi'], r, 1024)
        marg = m*(lp - tau) - bc
        if best is None or marg > best[0]:
            best = (marg, r, lp, bc)
    marg, r, lp, bc = best
    print(f'   c={c:2d} {"relax" if relax else "meas ":5s} m={m:3d} '
          f'sigma_m={float(d["sigma_m"]):.1f} tau^f={float(d["tau_flat"]):.4f} '
          f'tau^#={d["tau_sharp"]:.4f} tau={tau:.4f} | ceiling entry='
          f'{math.log(H["scale"])-tau:+.4f} | best r={r:.3f} log|phi\'|={lp:+.4f} '
          f'entry={lp-tau:+.4f} BC={bc:+.4f} margin={marg:+.2f}')


if __name__ == '__main__':
    for name, H in HOSTS.items():
        print(f'=== {name}   ceiling {H["scale"]} (log {math.log(H["scale"]):.4f});'
              f' folds: {H["folds"]} ===')
        rs = [H['rmax']*t for t in (0.2, 0.35, 0.5, 0.65, 0.8, 0.9, 0.97, 1.0)]
        for c in (4, 7, 10):
            for relax in (False, True):
                run(name, H, c, relax, rs)
        print()

    print('=== how good would the pure module have to be?  (host A) ===')
    print('   entry needs tau < log 4 = 1.3863 with sigma_m = 2 :')
    for th in (0.50, 0.554, 0.60, 0.70):
        print(f'     u_1/m = u_2/m = {th:.3f} -> tau^flat = {2-2*th*th:.4f}'
              f'  entry(ceiling) = {math.log(4)-(2-2*th*th):+.4f}')
    print('   i.e. STRICTLY MORE THAN 55.4% of all m functions must carry no lcm')
    print('   layer at all (type n^e only).  The conditional orbit never does.')
    print()
    print('=== reference: Catalan on the level-8 modular host (CDT_FINDER.md) ===')
    print('   unsymmetrised ceiling log(16/4) = %.4f, tau = 2.2756  -> entry %+0.4f'
          % (math.log(4), math.log(4)-2.2756))
    print('   symmetrised   ceiling log(256/4)= %.4f, tau = 4.2355  -> entry %+0.4f'
          % (math.log(64), math.log(64)-4.235459))
