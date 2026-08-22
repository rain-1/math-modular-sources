"""The best RIGOROUS-shape CDT criterion for a second-order Apery row.

margin(m, inventory, r) = (m-1) log|phi'(0)| - m tau + [ -Delta(r) ]
   with log|phi'(0)| = log(16 r |x_2|) = log(4|x_2|) + log(4r)
   => margin = (m-1)(log 4 + 2 + score) - m tau + [(m-1) log(4r) - Delta(r)]
   with score = log(1/|lam_2|) - 2 and log|x_2| = score + 2.
G(m) := max(0, max_r (m-1) log(4r) - Delta(r))   [0 = the Koebe map, Delta = 0]
threshold(m) := smallest score for which margin > 0
              = ( m*tau - G(m) )/(m-1) - log 4 - 2 .
"""
import json, math, sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
from cdt_bound import tau_flat, tau_sharp
from arch_k import inventory, tau_of

DEL = {float(k): v for k, v in json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                          'delta_table.json'))).items()}

def G(m):
    best, argr = 0.0, None          # r -> 0 limit / Koebe map: gain 0 by definition
    for r, d in sorted(DEL.items()):
        g = (m-1)*math.log(4*r) - d
        if g > best: best, argr = g, r
    return best, argr

def best_inventory(m):
    """all (p pure, J integrals) with p + 2 + J = m, p >= 1, J >= 0"""
    out = []
    for p in range(1, m-1):
        J = m - 2 - p
        mm, cols, e = inventory(J, p)
        assert mm == m, (mm, m)
        tf, ts, T = tau_of(mm, cols, e)
        out.append((T, p, J, tf, ts))
    return sorted(out)[0]

def m2():
    tf = tau_flat(2, [(1,1),(1,1)]); ts, _ = tau_sharp(2, [0,0])
    return float(tf[1])+ts, 1, 0, float(tf[1]), ts

if __name__ == '__main__':
    print(__doc__)
    print(f"{'m':>3} {'p':>2} {'J':>2} {'tau':>8} {'G(m)':>8} {'r*':>5} "
          f"{'thresh score (pure>=2 ok)':>26} {'thresh score (p=1 only)':>24}")
    rows = []
    for m in range(2, 21):
        Tb, pb, Jb, _, _ = m2() if m == 2 else best_inventory(m)
        # p = 1 variant (the only one available when lam_2 is not a rational integer)
        if m == 2:
            T1, p1, J1 = Tb, 1, 0
        else:
            mm, cols, e = inventory(m-3, 1); T1 = tau_of(mm, cols, e)[2]
        g, r = G(m)
        th_b = (m*Tb - g)/(m-1) - math.log(4) - 2
        th_1 = (m*T1 - g)/(m-1) - math.log(4) - 2
        rows.append((m, pb, Jb, Tb, g, r, th_b, th_1))
        print(f"{m:>3} {pb:>2} {Jb:>2} {Tb:>8.5f} {g:>8.5f} {str(r):>5} "
              f"{th_b:>26.5f} {th_1:>24.5f}")
    bb = min(rows, key=lambda z: z[6]); b1 = min(rows, key=lambda z: z[7])
    print(f"\nBEST (any inventory)   : m={bb[0]}, p={bb[1]}, J={bb[2]}, tau={bb[3]:.5f}, "
          f"r*={bb[5]}, threshold score = {bb[6]:+.5f}")
    print(f"BEST (p = 1, generic)  : m={b1[0]}, tau={b1[3]:.5f}, r*={b1[5]}, "
          f"threshold score = {b1[7]:+.5f}")
    print(f"For reference: the elementary (Beukers) criterion is  score > 0.")
    print(f"So the CDT route is worth  {-b1[7]:+.5f} nats over the elementary route "
          f"(generic inventory),")
    print(f"                     and   {-bb[6]:+.5f} nats with a pure polylogarithm module "
          f"(needs lam_2 in Z).")
