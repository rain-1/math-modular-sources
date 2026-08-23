# 10_score.py -- CDT scoring of the chi_{-4} rows found by 06-09.
# Conventions: consolidation/CDT_FINDER.md (tau, ceilings, contour loss),
# consolidation/CDT_NONCONGRUENCE.md (architectures K/D/S, per-function units).
#
# Ceiling = conformal radius at t=0 of the universal (orbifold) cover of
#   Omega = P^1 minus (Sigma minus {fold}), with 0 counted as a puncture.
# Exact values used:
#   P^1 minus {0,s,infty}          : r = 16|s|
#   P^1 minus {0,a,b} (infty free) : r = 16|a b /(b-a)|
#   {0,infty} punctures + one order-2 cone point at c : r = 64|c|
# For >3 remaining points we use monotonicity: r <= min over 3-subsets.
import sys, math, itertools
sys.path.insert(0, '/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from cdt_bound import tau as cdt_tau

M = 14; COLS = [(1,2),(3,2)]
EVEC = [0,0,1,0,0,0,0,0,0,1,1,1,1,1]
TAU_S = cdt_tau(M, COLS, EVEC)['tau']          # k=2, symmetrised (b_j=2)
TAU_D = cdt_tau(3, [(1,1),(1,1)], [0,0,0])['tau']   # m=3 inventory {1,H,theta H}, tau^flat=16/9
LOSS = math.log(0.62922326804)
BC0 = 11.845

INF = float('inf')
def r3(pts):
    """radius at 0 of the universal cover of P^1 minus ({0} union pts), |pts|=2."""
    a, b = pts
    if a == INF: a, b = b, a
    if b == INF: return 16*abs(a)
    return 16*abs(a*b/(b-a))

def rmulti(rest):
    """rigorous upper bound: min over 2-subsets of rest (excluding 0)."""
    pts = [p for p in rest if p != 0]
    best = INF
    for c in itertools.combinations(pts, 2):
        best = min(best, r3(c))
    return best

# rows: name, Sigma_finite_nonzero (list, complex ok), fold, S-architecture data
ROWS = [
 ("Gamma_0(8) po2   Zagier E  [CALIBRATION]", [0.125,0.25], 0.125,
   dict(kind='cone', desc='y=4x^2/(4x-1), free fixed point x=1/2 -> cone-2 at y=1', ceil=64.0)),
 ("Gamma_0(16) po2  F=theta^2 [CALIBRATION]", [0.25,0.5,0.25+0.25j,0.25-0.25j], 0.25,
   dict(kind='punct', desc='y=4x^2/(4x+1): both sigma-fixed points are cusps; Omega_y = P1-{0,-1/2,-1}', ceil=r3((-0.5,-1.0)))),
 ("Gamma_0(12) po1", [1/6,0.25,1/3,0.5], 1/6, None),
 ("Gamma_0(12) po2", [1/6,0.25,1/3,0.5], 1/6, None),
 ("Gamma_0(12) po3", [0.25,0.5,-0.5,1.0], 0.25,
   dict(kind='punct', desc='sigma(t)=1/(4t) swaps fold<->1: image stays a puncture; Omega_v=P1-{0,4/5,1,-1}',
        ceil=rmulti([4/5,1.0,-1.0]))),
 ("Gamma_0(12) po5", [0.25,0.5,-0.5,1.0], 0.25,
   dict(kind='punct', desc='same', ceil=rmulti([4/5,1.0,-1.0]))),
 ("Gamma_0(12) po4  [fold is a +-pair; no Apery limit]", [1/3,-1/3,1.0,-1.0], 1/3, None),
]

print(f"tau_S (CDT k=2 m=14 inventory) = {TAU_S:.6f}")
print(f"tau_D (m=3 unsymmetrised)      = {TAU_D:.6f}")
print()
hdr = (f"{'row':46s}{'l1':>7s}{'l2':>7s}{'#rem':>5s}{'score':>8s}"
       f"{'ceilD':>8s}{'lnD':>8s}{'entD':>8s}{'ceilS':>8s}{'lnS':>8s}{'entS':>8s}{'margS':>9s}")
print(hdr); print('-'*len(hdr))
for name, sing, fold, Sdat in ROWS:
    rest = [p for p in sing if abs(p-fold) > 1e-12]
    dd = []
    for p in rest:
        if not any(abs(p-q) < 1e-12 for q in dd): dd.append(p)
    mods = sorted(set(round(abs(p),12) for p in sing))
    lam1, lam2 = 1/abs(fold), 1/mods[1]
    nrem = len(dd)+2
    ceilD = rmulti(dd+[INF])
    entD = math.log(ceilD) - TAU_D
    if Sdat:
        ceilS = Sdat['ceil']; lnS = math.log(ceilS); entS = lnS - TAU_S
        BC = BC0 + math.log(1/lam2)
        margS = M*(entS+LOSS) - BC
        sS, sE, sM = f"{ceilS:8.3f}", f"{lnS:8.4f}", f"{entS:8.4f}"
        sMa = f"{margS:9.2f}"
    else:
        sS, sE, sM, sMa = "     n/a", "     n/a", "     n/a", "      n/a"
    print(f"{name:46s}{lam1:7.4f}{lam2:7.4f}{nrem:5d}{math.log(1/lam2)-2:8.4f}"
          f"{ceilD:8.3f}{math.log(ceilD):8.4f}{entD:8.4f}{sS}{sE}{sM}{sMa}")
print()
for name, sing, fold, Sdat in ROWS:
    if Sdat: print(f"  {name}: {Sdat['desc']}")
