"""tau(b;e) for the mu_4 architectures, and the entry / margin ledger.

Denominator descriptions are the MEASURED ones of mu4_series.py.  Each function
may be described in several ways (a type n^e is a fortiori of type [1..n]^e,
and a layer [1..n/2] is a fortiori [1..n]); we search over the descriptions to
minimise tau, since CDT's array b must be step-shaped column by column.
"""
import sys, os, itertools
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'cdt_finder'))
from fractions import Fraction as F
import math
from cdt_bound import tau_sharp

# ------------------------------------------------------------------ tau^flat
def tau_flat_rows(rows):
    """rows = list of layer-multisets (each a sorted tuple of rates b>0).
    Build the step-shaped array by sorting functions by sigma and padding each
    column upward to its maximal rate; return (sigma_m, tau^flat)."""
    m = len(rows)
    r = max((len(x) for x in rows), default=0)
    if r == 0:
        return F(0), F(0)
    # pad rows to length r with zeros, sort each row descending (largest layer first)
    padded = [tuple(sorted(x, reverse=True)) + (F(0),)*(r-len(x)) for x in rows]
    # order functions by increasing sigma (row sum)
    padded.sort(key=lambda t: sum(t))
    # step shape: column j must be 0..0,b_j..b_j -> b_j := max over the column,
    # and every nonzero entry is raised to b_j; every zero below a nonzero is raised too
    cols = []
    for j in range(r):
        col = [padded[i][j] for i in range(m)]
        bj = max(col)
        if bj == 0:
            continue
        # u_j = number of leading zeros in the (already sorted) column;
        # entries after the first nonzero are raised to b_j
        first = next(i for i in range(m) if col[i] != 0)
        cols.append((first, bj))
    sigma_m = sum(b for _, b in cols)
    tf = sigma_m - sum(F(u)**2*b for u, b in cols)/F(m)**2
    return sigma_m, tf, cols


def tau_of(inv):
    """inv = list of (name, e, layers).  Returns dict."""
    rows = [tuple(F(b) for b in lay) for (_, _, lay) in inv]
    e = [ee for (_, ee, _) in inv]
    m = len(inv)
    sm, tf, cols = tau_flat_rows(rows)
    ts, xi = tau_sharp(m, e)
    return dict(m=m, sigma_m=sm, cols=cols, tau_flat=tf, tau_sharp=ts, xi=xi,
                tau=float(tf)+ts, e=e)


def report(name, inv, ceiling, realised_loss=0.0, BC_shape=6.763):
    d = tau_of(inv)
    m = d['m']
    for lab, lp in (('ceiling', ceiling), ('realised', ceiling+realised_loss)):
        entry = lp - d['tau']
        BC = lp + BC_shape
        marg = m*entry - BC
        print(f'  {name:28s} m={m:3d}  sigma_m={float(d["sigma_m"]):.2f} '
              f'tau^f={float(d["tau_flat"]):.4f} tau^#={d["tau_sharp"]:.4f} '
              f'tau={d["tau"]:.4f} | {lab:8s} log|phi\'|={lp:+.4f} '
              f'entry={entry:+.4f} margin={marg:+.3f}')
    return d


# --------------------------------------------------------------- inventories
half = F(1, 2)

# x-line pure module, MEASURED types (mu4_series.py)
PURE_X = [
    ('1',                    0, []),
    ('arctan x',             1, []),
    ('log(1+x^2)',           1, []),
    ('Li_2(-x^2)',           2, []),
    ('f(x) = Ti_2(x)',       2, []),
    ('log^2(1+x^2)',         1, [half]),
    ('M(x)',                 1, [half]),
    ('arctan^2 x',           1, [1]),
    ('arctan x log(1+x^2)',  1, [1]),
    ('J(-x^2)',              0, [half, half]),
]
# relaxed alternative: describe the n^e functions as [1..n]^e (e -> 0)
PURE_X_RELAX = [
    ('1',                    0, []),
    ('arctan x',             0, [1]),
    ('log(1+x^2)',           0, [1]),
    ('Li_2(-x^2)',           0, [1, 1]),
    ('f(x) = Ti_2(x)',       0, [1, 1]),
    ('log^2(1+x^2)',         0, [1, half]),
    ('M(x)',                 0, [1, half]),
    ('arctan^2 x',           0, [1, 1]),
    ('arctan x log(1+x^2)',  0, [1, 1]),
    ('J(-x^2)',              0, [half, half]),
]
COND_X = ('H_j', 1, [1, 1])           # conditional orbit, x-line, type n[1..n]^2

if __name__ == '__main__':
    print('=== host A: x-line, P^1-{0,i,-i,oo}, mu_4 template ceiling |phi\'(0)| = 4 ===')
    ceilA = math.log(4)
    lossA = math.log(0.7304)          # concentric disc |zeta|<0.7304 (mu4_geom.py)
    for nc in (4, 7, 10):
        inv = PURE_X + [(f'H_{j}', 1, [1, 1]) for j in range(nc)]
        report(f'measured types, c={nc}', inv, ceilA, lossA)
    print()
    for nc in (4, 7, 10):
        inv = PURE_X_RELAX + [(f'H_{j}', 0, [1, 1]) for j in range(nc)]
        report(f'all-lcm types, c={nc}', inv, ceilA, lossA)
    print()
    print('  best conceivable: pure orbit entirely denominator-free (u_j = m/2)')
    for m in (14, 20, 30, 50):
        inv = [('free', 0, [])]*(m//2) + [('H', 0, [1, 1])]*(m - m//2)
        report(f'fantasy u=m/2, m={m}', inv, ceilA, lossA)


# ---------------------------------------------------------------------------
# Like-for-like comparison with Catalan's MODULAR host (Zagier row E, level 8).
# Pure module = the CDT five functions on P^1 - {0, 1/4, oo} (lambda_2 = 4);
# conditional row = the Eichler companion, type [1..n]^2.
LEVEL8_X = [('1', 0, []), ('log(1-4x)', 1, []), ('Li_2(4x)', 2, []),
            ('log^2(1-4x)', 1, [half]), ('J(4x)', 0, [1, 1])]
LEVEL8_X_RELAX = [('1', 0, []), ('log(1-4x)', 0, [1]), ('Li_2(4x)', 0, [1, 1]),
                  ('log^2(1-4x)', 0, [1, half]), ('J(4x)', 0, [1, 1])]
# CDT's own array for the B_i (Lemma bdenominators + (13.0.2)): u_1=1, u_2=3.
LEVEL8_Y = [('B_1', 0, []), ('B_2', 0, [2]), ('B_3', 1, [2]),
            ('B_4', 0, [2, 2]), ('B_5', 0, [2, 2]), ('B_6', 1, [2, 2]),
            ('B_7', 1, [2, 2])]


def level8():
    print('=== like-for-like: Catalan on the level-8 modular host (Zagier E) ===')
    for c in (4, 7, 9):
        inv = LEVEL8_X_RELAX + [(f'H_{j}', 0, [1, 1]) for j in range(c)]
        d = tau_of(inv)
        print(f'   unsymmetrised  c={c:2d} m={d["m"]:3d} sigma_m={float(d["sigma_m"]):.1f} '
              f'tau={d["tau"]:.4f}  ceiling log(16/4)={math.log(4):.4f} '
              f'entry={math.log(4)-d["tau"]:+.4f}')
    for c in (7,):
        inv = LEVEL8_Y + [(f'H_{j}', 1 if j >= 2 else 0, [2, 2]) for j in range(c)]
        d = tau_of(inv)
        print(f'   symmetrised    c={c:2d} m={d["m"]:3d} sigma_m={float(d["sigma_m"]):.1f} '
              f'tau={d["tau"]:.4f}  ceiling log(256/4)={math.log(64):.4f} '
              f'entry={math.log(64)-d["tau"]:+.4f}')
    print('   (CDT-transported inventory: tau=4.2355, entry -0.0766 -- CDT_FINDER.md)')
