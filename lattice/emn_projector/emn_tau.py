"""tau(b;e) and the entry ledger for the EMN host.

Coordinate: Z = z/2.  Singular set {0, 1/2, 1, oo}; fold Z=1 (the cusp z=2 where
H(2)=2G is an honest Apery limit); surviving outer singularity t2 = 1/2.
Width-law ceiling after fold removal: |phi'(0)| <= 16*|t2| = 8.
"""
import sys, os, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'catalan_mu4'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'cdt_finder'))
from fractions import Fraction as F
from mu4_tau import tau_flat_rows, tau_of

half = F(1,2)
CEIL = math.log(8.0)            # width law: log(16*|t2|), t2 = 1/2

# pure module on P^1-{0,1/2,oo} in Z (polylogs of 2Z), all regular at the fold Z=1
PURE = [
    ('1',                0, []),
    ('Li_1(2Z)',         1, []),
    ('Li_2(2Z)',         2, []),
    ('Li_3(2Z)',         3, []),
    ('log^2(1-2Z)',      1, [1]),
    ('Li_2 . log',       1, [1]),
    ('log^3(1-2Z)',      1, [1]),
    ('Li_{2,1}(2Z)',     1, [1]),
    ('J(2Z)  [CDT 5th]', 0, [half, half]),
    ('Li_4(2Z)',         4, []),
]
PURE_RELAX = [(n, 0, [F(1)]*e + list(l)) for (n, e, l) in PURE]
COND = ('B_j', 1, [F(1), F(2)])          # measured shape n[1..n][1..2n], sigma=3
COND_TRUE = ('B_j', 1, [F(1), F(1.95)])  # 'true size' 2.95 packed into two layers

def run(pure, cond, cmax=8, label=''):
    print(f"--- {label}")
    for c in range(1, cmax+1):
        inv = pure + [(f'B_{j}',)+cond[1:] for j in range(c)]
        d = tau_of(inv)
        entry = CEIL - d['tau']
        print(f"   m={d['m']:3d} (pure {len(pure)}, cond {c})  sigma_m={float(d['sigma_m']):.3f} "
              f"tau^f={float(d['tau_flat']):.4f} tau^#={d['tau_sharp']:.4f} tau={d['tau']:.4f}"
              f"   ceiling={CEIL:+.4f}   ENTRY={entry:+.4f}")

if __name__ == '__main__':
    print(f"EMN host, Z=z/2: ceiling log(16*1/2)= {CEIL:.4f}")
    run(PURE, COND, 8, 'measured types, shaped conditional n[1..n][1..2n]')
    run(PURE_RELAX, COND, 8, 'relaxed pure types (n^e -> [1..n]^e)')
    # sensitivity: what would a cheaper conditional row buy?
    print('--- sensitivity: conditional row of total lcm weight k (pure module as above)')
    for k in (2, 2.5, 3, 3.5, 4):
        inv = PURE + [(f'B_{j}', 1, [F(1), F(k-1)]) for j in range(4)]
        d = tau_of(inv)
        print(f"   k={k}: tau={d['tau']:.4f} entry={CEIL-d['tau']:+.4f}")
