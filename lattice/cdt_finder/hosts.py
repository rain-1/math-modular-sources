"""CDT finder: score every modular Apery host in the project census against the
Calegari-Dimitrov-Tang arithmetic holonomy bound (arXiv:2408.15403, Thm 7.0.1).

Model (calibrated on CDT's own L(2,chi_{-3}) proof; see CDT_FINDER.md for the
derivation and for every place where an input is estimated rather than computed):

  * host coordinate  x = q + O(q^2),  A_n in Z,  D_n^k B_n in Z  (k = w+1);
    singular set of the Picard-Fuchs operator = {0, 1/lam_1, 1/lam_2, infty}.
  * CDT normalise so that the OUTER singularity sits at s := 1/lam_2 and use the
    "pure" polylogarithm module on P^1 \ {0, s, infty}.  Integrality of the pure
    functions Li_j(x/s) = sum (lam_2 x)^n / n^j forces  s = 1/lam_2  with
    lam_2 a nonzero RATIONAL INTEGER; CDT have lam_2 = 1.
  * normaliser descent:  w(x) = s x/(x - s),  y = x + w(x) = x^2/(x - s),
    fixing x = 0 and exchanging s <-> infty; y-singularities {0, 4s, infty, y(1/lam_1)}.
  * uniformisation:  x = s*lambda(z)  (lambda = 16 z - 128 z^2 + ...), so the
    symmetrised coordinate has  y = -256 s z^2 + ...  and the CEILING on the
    conformal size is  |phi'(0)| <= 256 s = 256/lam_2.
    CDT realise 161.081 = 256 * 0.62922 (Lemma A.4.4).  Because the ambient
    orbifold is the SAME for every host (only the scale s changes), a contour of
    CDT's shape rescales as phi -> s*phi, giving
        log|phi'(0)| = log 161.081 - log lam_2,
        BC-integral  = 11.845     - log lam_2.
  * denominator array: the conditional orbit carries all k LCM layers doubled by
    the symmetrisation (b_j = 2), the pure orbit carries u_j of them; m = 2*p_0.

  ==> margin(m, k, u, lam_2) = m*(log|phi'(0)| - tau(b;e)) - BC.
"""
import json, math, os, sys
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cdt_bound import tau_flat, tau_sharp, tau  # noqa

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))

# --- CDT's realised contour, on the s = 1 (lam_2 = 1) ambient -----------------
LOG_PHI_CDT = math.log(256*5448339453535586608000000000
                       / 8658833407565631122430056127)      # 5.0819077...
BC_CDT      = 11.845          # Bost-Charles integral, CDT (A.5.1)
BC_BEST     = 11.529          # their best convexity-improved numerator (13.621)


def cdt_model(k, p_0, lam2, u=None, e_count=None, bc=BC_CDT):
    """k     = Eichler denominator exponent of the conditional function (w+1);
       p_0   = size of each orbit, m = 2*p_0;
       lam2  = outer characteristic root (|lam2| >= 1, integer);
       u     = list of k integers, u_j = #functions missing the j-th LCM layer
               (default: CDT's pattern rescaled, u_j = (2j-1)*p_0/7 for k=2).
       e_count = sum of the integration vector (default: CDT's 6 scaled)."""
    m = 2*p_0
    if u is None:
        u = [round((2*j - 1)*p_0/7.0) for j in range(1, k + 1)]
    cols = [(uj, 2) for uj in u]
    if e_count is None:
        e_count = round(6*p_0/7.0)
    e = [1]*e_count + [0]*(m - e_count)
    T = tau(m, cols, e)
    L = LOG_PHI_CDT - math.log(abs(lam2))
    N = bc - math.log(abs(lam2))
    den = L - T['tau']
    return dict(m=m, k=k, p0=p_0, lam2=lam2, u=u, sigma_m=float(T['sigma_m']),
                tau_flat=float(T['tau_flat']), tau_sharp=T['tau_sharp'],
                tau=T['tau'], ceiling=math.log(256.0/abs(lam2)),
                logphi=L, bc=N, entry=den,
                bound=(N/den if den > 0 else float('inf')),
                margin=m*den - N)


def is_int(x, tol=1e-6):
    return abs(x - round(x)) < tol and abs(round(x)) >= 1


def load_rows():
    """paper census (hand-curated, exact) + sporadic_scan2 (85 machine rows)."""
    census = [
        # name, level, period, r(order), lam1, lam2, k(denominator exponent)
        ('Zagier A (7,2,-8)',    '6',   'zeta(2)/4',              2,  8.0,   -1.0,  2),
        ('Zagier B (9,3,27)',    '36',  'none (complex fold)',    2,  5.196, None, 2),
        ('Zagier C (10,3,9)',    '6',   'L(2,chi_-3)/2',          2,  9.0,    1.0,  2),
        ('Zagier D (11,3,-1)',   '11',  'zeta(2)/5',              2, 11.09,  -0.0902, 2),
        ('Zagier E (12,4,32)',   '8',   'G = L(2,chi_-4)/2',      2,  8.0,    4.0,  2),
        ('Zagier F (17,6,72)',   '6/8', '5L(2,chi_-3)/8',         2,  9.0,    8.0,  2),
        ('Cooper s_7',           '7',   'zeta(2)/7',              3, 27.0,   -1.0,  2),
        ('Cooper s_10',          '10',  'zeta(2)/5',              3, 16.0,   -4.0,  2),
        ('Cooper s_18',          '18',  'L(2,chi_-3)/2',          3, 16.0,   12.0,  2),
        ('AZ (7,3,81)',          'B-IV','none (complex fold)',    3,  9.0,   None, 3),
        ('AZ (9,3,-27)',         '-',   'L(3,chi_-3)/3',          3, 19.39, -1.392, 3),
        ('Domb (10,4,64)',       '6',   '7 zeta(3)/24',           3, 16.0,   4.0,  3),
        ('eta (11,5,125)',       'ord6','L(3,chi_5)/2',           3, 11.18, None, 3),
        ('T (12,4,16)',          '6',   '7 zeta(3)/32',           3, 23.31, 0.6863, 3),
        ('Apery (17,5,1)',       '5/6', 'zeta(3)/6',              3, 33.97, 0.02944, 3),
        ('zeta(7) level 24',     '24',  '1463 zeta(7)/13824',     7,  6.828, 4.828, 7),
    ]
    out = []
    for nm, lev, per, r, l1, l2, k in census:
        out.append(dict(src='paper census', name=nm, level=lev, period=per,
                        w=r-1, lam1=l1, lam2=l2, k=k))
    p = os.path.join(ROOT, 'lattice/sporadic_scan2/table.json')
    if os.path.exists(p):
        for row in json.load(open(p)):
            if row.get('cplx') or row.get('degen'):
                continue
            l2 = row.get('lam2')
            if l2 is None or row.get('k') in (None, 0):
                continue
            out.append(dict(src='scan2', name=f"scan2 #{row['ti']} (N={row['N']},w={row['w']})",
                            level=str(row['N']), period=(row.get('known') or row.get('limit') or '?'),
                            w=row['w'], lam1=row.get('lam1'), lam2=l2, k=row['k']))
    return out


if __name__ == '__main__':
    print('=== calibration: CDT level-6 weight-3, their own numbers ===')
    r = cdt_model(k=2, p_0=7, lam2=1, u=[1, 3], e_count=6)
    print(f"  m={r['m']}  tau^flat={r['tau_flat']:.6f} (191/49={191/49:.6f})  "
          f"tau^#={r['tau_sharp']:.6f} (27/80=0.3375)  tau={r['tau']:.6f}")
    print(f"  ceiling log(256/lam2)={r['ceiling']:.4f}  realised log|phi'(0)|={r['logphi']:.4f}"
          f"  BC={r['bc']:.3f}  entry={r['entry']:.6f}  bound m<={r['bound']:.4f}  margin={r['margin']:+.4f}")
    assert abs(r['bound'] - 13.9938) < 2e-3, r['bound']
    print('  CALIBRATION OK (CDT: m <= 13.9938 < 14).')

    print()
    print('=== ranked hosts (CDT-identical architecture: p_0=7, u=[1,3,..], e=6) ===')
    rows = load_rows()
    res = []
    for h in rows:
        l2 = h['lam2']
        if l2 is None:
            continue
        k = h['k']
        if k < 1 or k > 8:
            continue
        rational = is_int(l2)                # integral outer root: pure polylogs integral
        lam2eff = abs(l2)
        u = [1, 3] if k == 2 else [round((2*j-1)*7/7.0) for j in range(1, k+1)]
        try:
            m = cdt_model(k=k, p_0=7, lam2=max(lam2eff, 1e-9), u=u, e_count=6)
        except Exception as ex:
            continue
        res.append((m['margin'], h, m, rational))
    res.sort(key=lambda z: -z[0])
    seen = set()
    print(f"{'host':38s} {'lvl':>5s} {'w':>2s} {'k':>2s} {'lam2':>9s} {'int?':>5s} "
          f"{'ceil':>7s} {'logphi':>7s} {'tau':>7s} {'entry':>8s} {'m<=':>8s} {'margin':>9s}  period")
    for mar, h, mm, rat in res:
        key = (h['name'])
        if key in seen:
            continue
        seen.add(key)
        print(f"{h['name'][:38]:38s} {h['level']:>5s} {h['w']:>2d} {h['k']:>2d} "
              f"{h['lam2']:>9.4f} {('yes' if rat else 'NO'):>5s} "
              f"{mm['ceiling']:>7.3f} {mm['logphi']:>7.3f} {mm['tau']:>7.3f} "
              f"{mm['entry']:>8.4f} {mm['bound']:>8.3f} {mar:>+9.3f}  {h['period'][:28]}")
