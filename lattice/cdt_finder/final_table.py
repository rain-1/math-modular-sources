"""Final ranked CDT-finder table (see consolidation/CDT_FINDER.md)."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cdt_bound import tau_flat, tau_sharp

LOSS = math.log(0.6292232680)     # CDT's realised contour loss (Lemma A.4.4)
BC   = 11.845                     # CDT's Bost-Charles integral at scale s = 1

HOSTS = [
 # name, level, w, k, lam2_norm (|N(lam2)|^{1/deg}), rational?, period, note
 ('Zagier C  (10,3,9)',      '6',  1, 2,  1.0, 'Q',       'L(2,chi_-3)/2', 'CDT Theorem A'),
 ('Zagier A  (7,2,-8)',      '6',  1, 2,  1.0, 'Q',       'zeta(2)/4',     'CDT second orientation'),
 ('Cooper s_7',              '7',  2, 2,  1.0, 'Q',       'zeta(2)/7',     'k=2 by free integration'),
 ('Zagier D  (11,3,-1)',     '11', 1, 2,  1.0, 'Q(sqrt5) unit', 'zeta(2)/5', 'Apery-perfect fold'),
 ('Zagier E  (12,4,32)',     '8',  1, 2,  4.0, 'Q',       'G = Catalan',   'nearest miss'),
 ('Cooper s_10',             '10', 2, 2,  4.0, 'Q',       'zeta(2)/5',     ''),
 ('Zagier F  (17,6,72)',     '6/8',1, 2,  8.0, 'Q',       '5L(2,chi_-3)/8',''),
 ('Cooper s_18',             '18', 2, 2, 12.0, 'Q',       'L(2,chi_-3)/2', ''),
 ('scan2 #942/3/4 (N=16)',   '16', 1, 2,  2.0, 'Q',       '?',             ''),
 ('Apery  (17,5,1)',         '5/6',2, 3,  1.0, 'Q(sqrt2) unit', 'zeta(3)/6', 'c=1'),
 ('X_1(5) Sym^2 (Beukers)',  '5',  2, 3,  1.0, 'Q(sqrt5) unit', '8z(3)-5r5L(3,chi5)', 'PRIORITY CASE'),
 ('T  (12,4,16)',            '6',  2, 3,  4.0, 'Q(sqrt2)','7 zeta(3)/32',  ''),
 ('Domb (10,4,64)',          '6',  2, 3,  4.0, 'Q',       '7 zeta(3)/24',  ''),
 ('AZ (9,3,-27)',            '9',  2, 3,  5.196,'Q(sqrt3)','L(3,chi_-3)/3', ''),
 ('zeta(7) level 24',        '24', 6, 7,  4.828,'Q(sqrt2)','1463 z(7)/13824',''),
]

def tau_for(m, k, u):
    sm, tf = tau_flat(m, [(uj, 2) for uj in u])
    ec = round(6*m/14.0)
    ts, _ = tau_sharp(m, [1]*ec+[0]*(m-ec))
    return float(tf)+ts

def score(lam2, k, m, mode):
    s = 1.0/lam2
    if mode == 'cdt':          # CDT-proportional pure inventory  u_j = (2j-1)m/14
        u = [max(0, min(m, round((2*j-1)*m/14.0))) for j in range(1, k+1)]
    else:                      # best case: whole pure orbit denominator-free
        u = [m//2]*k
    T = tau_for(m, k, u)
    ceil_ = math.log(256*s)
    real = ceil_ + LOSS
    N = BC + math.log(s)
    return dict(u=u, tau=T, ceil=ceil_, real=real, entry_c=ceil_-T,
                entry_r=real-T, margin=m*(real-T)-N, N=N)

if __name__ == '__main__':
    for mode, m, lab in (('cdt', 14, 'A. CDT-identical architecture (m=14, p_0=7, CDT pure inventory)'),
                         ('best', 14, 'B. m=14, best conceivable pure inventory (all 7 pure denominator-free)'),
                         ('best', 50, 'C. m=50, best conceivable pure inventory (p_0=25)')):
        print('\n' + lab)
        print(f"{'host':26s}{'lvl':>5s}{'w':>3s}{'k':>3s}{'|lam2|_norm':>12s}{'field':>16s}"
              f"{'tau':>8s}{'ceil':>8s}{'entryC':>8s}{'entryR':>8s}{'margin':>10s}  period")
        rows = []
        for nm, lev, w, k, l2, fld, per, note in HOSTS:
            r = score(l2, k, m, mode)
            rows.append((r['margin'], nm, lev, w, k, l2, fld, per, r, note))
        rows.sort(key=lambda z: -z[0])
        for mar, nm, lev, w, k, l2, fld, per, r, note in rows:
            print(f"{nm:26s}{lev:>5s}{w:>3d}{k:>3d}{l2:>12.4f}{fld:>16s}"
                  f"{r['tau']:>8.3f}{r['ceil']:>8.3f}{r['entry_c']:>+8.3f}{r['entry_r']:>+8.3f}"
                  f"{mar:>+10.3f}  {per}")
