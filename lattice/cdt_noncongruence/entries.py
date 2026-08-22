"""Entry conditions log|phi'(0)| > tau, in the same presentation as
CDT_FINDER.md Sec.4 (entryC = at the hard ceiling, entryR = with CDT's realised
contour loss 0.62922), for every architecture and every host."""
import math, os, sys
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'cdt_finder'))
from cdt_bound import tau_flat, tau_sharp
from arch_k import inventory, tau_of
from optimum import G
from table_nc import LOG_LOSS, BC_CDT, TAU_CDT, arch_S, adelic_gamma
from hosts_nc import HOSTS

TAU3 = tau_of(*inventory(0, 1))[2]          # m = 3, {1, H, theta H}
TAU2 = float(tau_flat(2, [(1,1),(1,1)])[1])
g3, r3 = G(3)
LOGR3 = math.log(4*r3)                       # log|phi_r'(0)| - log(4|x_2|)
DELTA3 = 2*LOGR3 - g3                        # BC - log|phi'(0)| for the m=3 optimum

if __name__ == '__main__':
    print(f"m=2 inventory {{1,H}}      : tau = {TAU2:.6f}")
    print(f"m=3 inventory {{1,H,thH}}  : tau = {TAU3:.6f} = 16/9")
    print(f"(D) optimum at m=3        : r* = {r3}, log|phi'(0)| = log(4|x_2|) {LOGR3:+.6f}, "
          f"Delta = BC - log|phi'(0)| = {DELTA3:.6f}")
    print(f"(S) CDT                   : tau = {TAU_CDT:.6f}, BC = 11.845 + log s, "
          f"log(loss) = {LOG_LOSS:.6f}\n")
    print(f"{'host':<36s}{'log 4|x2|':>10s}{'(K)entry':>9s}{'(D)entry':>9s}  | "
          f"{'ceil_S':>8s}{'(S)entryC':>10s}{'(S)entryR':>10s}{'gamma':>8s}"
          f"{'(S)adC':>9s}{'(S)adR':>9s}")
    for h in HOSTS:
        L0 = math.log(4*abs(h['x2']))
        eK = L0 - TAU3
        eD = L0 + LOGR3 - TAU3
        ap = h['lam2_rational']
        gam, _ = adelic_gamma(h['lam2'] if ap else None)
        ceil_ = math.log(256/h['lam2_norm'])
        print(f"{h['name']:<36s}{L0:>10.4f}{eK:>+9.4f}{eD:>+9.4f}  | "
              f"{ceil_:>8.4f}{ceil_-TAU_CDT:>+10.4f}{ceil_+LOG_LOSS-TAU_CDT:>+10.4f}"
              f"{gam:>+8.4f}{ceil_-TAU_CDT+gam:>+9.4f}{ceil_+LOG_LOSS-TAU_CDT+gam:>+9.4f}")
    print("\n(K)entry uses the Koebe map (|phi'(0)| = 4|x_2|); (D)entry the Kodaira map phi_{r*}.")
    print("An entry <= 0 means CDT's theorem does not even apply in that architecture.")
