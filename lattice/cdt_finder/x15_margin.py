"""X_1(5) Sym^2 (Beukers 1987 Thm 4 configuration) scored against CDT."""
import math, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cdt_bound import tau_flat, tau_sharp
from conformal import r_two_punctures
import mpmath as mp

phi = (1+math.sqrt(5))/2
t1, t2 = phi**-5, -phi**5
print(f't1 = phi^-5 = {t1:.10f}   t2 = -phi^5 = {t2:.10f}   t1*t2 = {t1*t2:.6f}  |N(t2)|=1')

# --- conformal sizes, per archimedean place ---------------------------------
print('\nconformal sizes |phi\'(0)| (base point t=0):')
for lab, s in (('v1 (sqrt5 -> +2.236)', abs(t2)), ('v2 (sqrt5 -> -2.236)', abs(t1))):
    print(f'  {lab}: outer singularity |s|={s:.6f}')
    print(f'     univalent slit plane C\\(-inf,s]   (Koebe)   4|s|      = {4*s:9.4f}  log={math.log(4*s):+.4f}')
    print(f'     thrice-punctured C\\{{t1,t2}} (exact, lambda)          = '
          f'{float(r_two_punctures(t1,t2)) if s==abs(t2) else float(r_two_punctures(t1,t2)):9.4f}')
    print(f'     multivalent modular ceiling, unsymmetrised 16|s|      = {16*s:9.4f}  log={math.log(16*s):+.4f}')
    print(f'     multivalent modular ceiling, symmetrised  256|s|      = {256*s:9.4f}  log={math.log(256*s):+.4f}')
print(f'\n  exact Poincare radius of C\\{{t1,t2}} at 0 = {float(r_two_punctures(t1,t2)):.6f}'
      f'   log = {math.log(float(r_two_punctures(t1,t2))):.4f}')

def budget(kind):
    if kind == 'sum':      # unnormalised sum over the two real places
        return lambda f: math.log(f*abs(t2)) + math.log(f*abs(t1))
    return lambda f: 0.5*(math.log(f*abs(t2)) + math.log(f*abs(t1)))   # normalised

print('\n--- entry test  log|phi\'(0)| > tau(b;e),  k = 3 (sharp, no free integration) ---')
LOSS = math.log(0.6292232680)     # CDT's realised contour loss factor
for m, u, ecount, lab in ((14, [1, 3, 5], 6, 'CDT-proportional pure inventory'),
                          (14, [7, 7, 7], 6, 'best case: whole pure orbit denominator-free'),
                          (50, [25, 25, 25], 21, 'best case, p_0 = 25'),
                          (60, [30, 30, 30], 26, 'best case, p_0 = 30')):
    sm, tf = tau_flat(m, [(uj, 2) for uj in u])
    ts, _ = tau_sharp(m, [1]*ecount+[0]*(m-ecount))
    T = float(tf)+ts
    for nrm, nlab in ((0.5, 'normalised (average over places, [K:Q]=2)'),
                      (1.0, 'unnormalised (sum over places)')):
        L_ceil = nrm*(math.log(256*abs(t2))+math.log(256*abs(t1)))
        L_real = L_ceil + (2*nrm)*LOSS
        N = 11.845 + math.log(abs(t2)) if nrm == 1.0 else 11.845   # BC integral, same shape
        if nrm == 1.0: N = (11.845+math.log(abs(t2))) + (11.845+math.log(abs(t1)))
        print(f'  m={m:3d} u={str(u):12s} tau={T:.4f} | {nlab[:28]:28s} '
              f'ceil={L_ceil:6.3f} real={L_real:6.3f} entry={L_real-T:+7.3f} '
              f'margin={m*(L_real-T)-N:+9.3f}')
