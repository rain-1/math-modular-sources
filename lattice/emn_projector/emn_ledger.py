"""Like-for-like entry ledger: the CDT five-function pure module transported to the
EMN host, plus c conditional rows, exactly as in mu4_tau.level8()."""
import sys, math
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/catalan_mu4')
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from fractions import Fraction as F
from mu4_tau import tau_of
half = F(1,2)

PURE5 = [('1',0,[]), ('log(1-2Z)',1,[]), ('Li_2(2Z)',2,[]),
         ('log^2(1-2Z)',1,[half]), ('J(2Z)',0,[1,1])]
PURE5R = [('1',0,[]), ('log(1-2Z)',0,[1]), ('Li_2(2Z)',0,[1,1]),
          ('log^2(1-2Z)',0,[1,half]), ('J(2Z)',0,[1,1])]

def run(label, pure, cond, ceil_, cs=(4,7,9)):
    print(f'--- {label}   ceiling={ceil_:+.4f}')
    for c in cs:
        inv = pure + [(f'B_{j}',)+cond for j in range(c)]
        d = tau_of(inv)
        print(f'    c={c:2d} m={d["m"]:3d} sigma_m={float(d["sigma_m"]):.3f} '
              f'tau^f={float(d["tau_flat"]):.4f} tau^#={d["tau_sharp"]:.4f} '
              f'tau={d["tau"]:.4f}   ENTRY={ceil_-d["tau"]:+.4f}')

C = math.log(8.0)
print('EMN host, Z=z/2, punctures {0,1/2,1,oo}; fold Z=1 (z=2, H(2)=2G); t2=1/2')
run('measured types, conditional n[1..n][1..2n] (k=3)', PURE5, (1,[F(1),F(2)]), C)
run('relaxed pure, conditional n[1..n][1..2n]', PURE5R, (0,[F(1),F(2)]), C)
print()
print('counterfactuals:')
run('if the conditional row were n[1..n]^2 (k=2)', PURE5, (1,[F(1),F(1)]), C)
run('if the fold sat at the outer singularity (ceiling log16, k=2.95)',
    PURE5, (1,[F(1),F(1.95)]), math.log(16.0))
print()
print('benchmarks: CDT L(2,chi_-3) +0.77 | level 8 sym -0.077 | mu_4 host A -0.563 | level16 -1.46')
