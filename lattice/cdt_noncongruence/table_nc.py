r"""Final margins table: every second-order Apery row of the project, scored in
four architectures and reduced to one comparable number -- the THRESHOLD ON THE
ELEMENTARY SCORE at which each architecture starts to give a contradiction.

  (E)  elementary (Beukers/Apery):  d_n^k |a_n xi - b_n| -> 0.        score > 0.
  (K)  CDT, univalent Koebe map of C \ [x_2,oo); BC = log(4|x_2|) exactly
       (Grunsky).  Fully rigorous shape: no contour has to be designed.
  (D)  CDT, multivalent Kodaira map phi_r(z) = x_2 lambda(r z); BC computed
       numerically (bc_multivalent.py -> delta_table.json).
  (S)  CDT's own symmetrised architecture: normaliser descent y = x^2/(x-s),
       s = 1/lam_2; ceiling 256 s; denominators doubled (b_j = 2); m = 14,
       u = (1,3), e-count 6; CDT's realised contour loss 0.62922 and
       BC = 11.845 + log s.  Needs s in Q; over a number field the
       CDT_FINDER.md Sec.3 normalisation replaces |lam_2| by |N(lam_2)|^{1/[K:Q]}.

margin (CDT_FINDER units) = m(log|phi'(0)| - tau) - BC ;  it equals
(m-1)*(score - threshold) in (K)/(D) and 13*(score - threshold) in (S).
"""
import json, math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'adelic_holonomy'))
from cdt_bound import tau_flat, tau_sharp
from adelic_bound import gamma_p
from arch_k import inventory, tau_of
from optimum import G, best_inventory, m2
from hosts_nc import HOSTS

LOG_LOSS = math.log(5448339453535586608000000000/8658833407565631122430056127)
BC_CDT   = 11.845
TAU_CDT  = float(tau_flat(14, [(1,2),(3,2)])[1]) + tau_sharp(14, [0,0,1,0,0,0,0,0,0,1,1,1,1,1])[0]
L4 = math.log(4)

def inv_list(m, allow_pure):
    if m == 2: return [m2()]
    if allow_pure: return [best_inventory(m)]
    mm, cols, e = inventory(m-3, 1)
    T = tau_of(mm, cols, e)[2]
    return [(T, 1, m-3, 0, 0)]

def arch_KD(logx2, allow_pure, multivalent, mmax=20):
    """min over m of the score-threshold; returns that m and its margin."""
    best = None
    for m in range(2, mmax+1):
        for T, p, J, _, _ in inv_list(m, allow_pure):
            g, r = (G(m) if multivalent else (0.0, None))
            thresh = (m*T - g)/(m-1) - L4 - 2
            score = logx2 - 2 + 2                     # placeholder; score = logx2-2
            margin = (m-1)*(logx2 + L4) - m*T + g
            logphi = logx2 + L4 + (math.log(4*r) if r else 0.0)
            entry  = logphi - T
            rec = dict(m=m, p=p, J=J, tau=T, r=r, gain=g, thresh=thresh,
                       margin=margin, entry=entry, logphi=logphi, BC=logphi - margin + (m-1)*logphi - 0)
            rec['BC'] = m*entry - margin
            if best is None or thresh < best['thresh']: best = rec
    return best

def arch_S(lam2_norm, gamma=0.0, loss=True):
    X = -math.log(lam2_norm)                 # = log s
    c = math.log(256) + (LOG_LOSS if loss else 0.0)
    margin = 13*X + 14*(c - TAU_CDT + gamma) - BC_CDT
    thresh = (14*(TAU_CDT - gamma) + BC_CDT - 14*c)/13 - 2
    return dict(m=14, tau=TAU_CDT, gamma=gamma, margin=margin, thresh=thresh,
                entry=c + X - TAU_CDT + gamma, logphi=c + X, BC=BC_CDT + X)

def adelic_gamma(lam2):
    """pure orbit of p-adic slope v_p(lam_2) in x, 2 v_p(lam_2) in y (CDT 7+7)."""
    if lam2 is None: return 0.0, {}
    l2 = int(round(abs(lam2))); g = 0.0; d = {}
    if abs(abs(lam2) - l2) > 1e-9 or l2 == 0: return 0.0, {}
    for p in (2,3,5,7,11,13):
        v = 0
        while l2 % p == 0: l2 //= p; v += 1
        if v:
            gp = gamma_p(14, [v]*7 + [0]*7, p); g += gp; d[p] = (v, gp)
    return g, d

if __name__ == '__main__':
    print(__doc__)
    print(f"tau_CDT = {TAU_CDT:.6f} (=16603/3920={16603/3920:.6f});  log(contour loss) = {LOG_LOSS:.6f}")
    print(f"(S) threshold on the score: {arch_S(1.0)['thresh']:+.6f}  "
          f"(at the hard ceiling: {arch_S(1.0, loss=False)['thresh']:+.6f})")
    print(f"    i.e. (S) needs lam2_norm <= {math.exp(-(arch_S(1.0)['thresh']+2)):.6f} "
          f"-- the Apery-perfect / unit-fold condition |N(lam_2)| = 1.")
    print(f"(K) threshold, generic (p=1): {arch_KD(0,False,False)['thresh']:+.5f}   "
          f"with a pure module: {arch_KD(0,True,False)['thresh']:+.5f}")
    print(f"(D) threshold, generic (p=1): {arch_KD(0,False,True)['thresh']:+.5f}   "
          f"with a pure module: {arch_KD(0,True,True)['thresh']:+.5f}\n")

    print(f"{'host':<36s}{'score':>8s}{'l2^n':>7s}  | "
          f"{'(K) m':>6s}{'margin':>9s}{'defic':>8s} | "
          f"{'(D) m':>6s}{'r':>5s}{'margin':>9s}{'defic':>8s} | "
          f"{'(S) margin':>11s}{'defic':>8s}{'+adelic':>9s}")
    rows = []
    for h in HOSTS:
        logx2 = math.log(abs(h['x2'])); sc = h['score']
        ap = h['lam2_rational']
        K = arch_KD(logx2, ap, False); D = arch_KD(logx2, ap, True)
        gam, gd = adelic_gamma(h['lam2'] if ap else None)
        S = arch_S(h['lam2_norm']); Sa = arch_S(h['lam2_norm'], gamma=gam)
        rows.append((h, K, D, S, Sa, gam, gd))
        print(f"{h['name']:<36s}{sc:>+8.4f}{h['lam2_norm']:>7.3f}  | "
              f"{K['m']:>6d}{K['margin']:>+9.3f}{sc-K['thresh']:>+8.4f} | "
              f"{D['m']:>6d}{str(D['r']):>5s}{D['margin']:>+9.3f}{sc-D['thresh']:>+8.4f} | "
              f"{S['margin']:>+11.3f}{sc-S['thresh']:>+8.4f}{Sa['margin']:>+9.3f}")
    print("\n'defic' = score - threshold: positive means that architecture gives a contradiction.")
    print("\nadelic gains (architecture S only; the (K)/(D) inventories have no function of positive slope):")
    for h, K, D, S, Sa, gam, gd in rows:
        if gam:
            print(f"  {h['name']:<36s} " + ", ".join(f"p={p}: v={v}, gamma={g:+.4f}" for p,(v,g) in gd.items())
                  + f"   margin {S['margin']:+.3f} -> {Sa['margin']:+.3f}, deficit {h['score']-S['thresh']:+.4f} -> {h['score']-Sa['thresh']:+.4f}")
    print("\nCALIBRATION  (CDT's own host, Zagier C):", 
          f"(S) margin = {arch_S(1.0)['margin']:+.4f}   [CDT: +0.0053]")

# ---------------------------------------------------------------------------
def final_report():
    print("\n" + "="*118)
    print("FINAL RANKING -- deficit in SCORE UNITS (nats).  'proved' = deficit > 0.")
    print("  (K) univalent Koebe, threshold %+.5f | (D) Kodaira phi_r, threshold %+.5f"
          % (arch_KD(0,False,False)['thresh'], arch_KD(0,False,True)['thresh']))
    print("  (D+) with a pure polylogarithm module (needs lam_2 in Z), threshold %+.5f"
          % arch_KD(0,True,True)['thresh'])
    print("  (S) CDT symmetrised, m=14: threshold %+.5f (realised contour) / %+.5f (hard ceiling)"
          % (arch_S(1)['thresh'], arch_S(1, loss=False)['thresh']))
    print("      -- (S) scores the host by its NORMALISED score -log|N(lam_2)|^{1/[K:Q]} - 2, not by 'score'.\n")
    print(f"{'host':<36s}{'score':>8s}{'sc_S':>8s} |{'(K)':>8s}{'(D)':>8s}{'(D+)':>8s} |"
          f"{'(S)real':>9s}{'(S)ceil':>9s}{'(S)ad-r':>9s}{'(S)ad-c':>9s} | best")
    out = []
    for h in HOSTS:
        logx2 = math.log(abs(h['x2'])); sc = h['score']; ap = h['lam2_rational']
        scS = -math.log(h['lam2_norm']) - 2
        gam, gd = adelic_gamma(h['lam2'] if ap else None)
        dK  = sc - arch_KD(logx2, False, False)['thresh']
        dD  = sc - arch_KD(logx2, False, True)['thresh']
        dDp = (sc - arch_KD(logx2, True, True)['thresh']) if ap else None
        dSr = scS - arch_S(1)['thresh']; dSc = scS - arch_S(1, loss=False)['thresh']
        dSar= scS - arch_S(1, gamma=gam)['thresh']
        dSac= scS - arch_S(1, gamma=gam, loss=False)['thresh']
        cand = [('K',dK),('D',dD)] + ([('D+',dDp)] if dDp is not None else []) + \
               [('S',dSr),('S-ceil',dSc),('S+ad',dSar),('S+ad-ceil',dSac)]
        bn, bv = max(cand, key=lambda z: z[1])
        out.append((h, bn, bv))
        print(f"{h['name']:<36s}{sc:>+8.4f}{scS:>+8.4f} |{dK:>+8.4f}{dD:>+8.4f}"
              f"{('%+8.4f'%dDp) if dDp is not None else '       -':>8s} |"
              f"{dSr:>+9.4f}{dSc:>+9.4f}{dSar:>+9.4f}{dSac:>+9.4f} | {bn} {bv:+.4f}")
    print("\nranked by best deficit:")
    for h, bn, bv in sorted(out, key=lambda z: -z[2]):
        print(f"  {bv:>+8.4f}  [{bn:<9s}] {h['name']:<36s}  {h['cong']:<24s} {h['period'][:40]}")

final_report()
