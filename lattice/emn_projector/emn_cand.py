"""Build the series of an arbitrary projector combination and measure it."""
from fractions import Fraction as Fr
import math
from emn_core import psi_basis, ratseries, polymul, solve_from_psi

def psi_series(comb, N):
    """comb: dict (a,kind)->Fr.  kind in {'0','p3'}.  Returns Taylor coeffs of Psi."""
    tot = [Fr(0)]*(N+1)
    for (a, kind), c in comb.items():
        if c == 0: continue
        num, den = psi_basis(a, '0' if kind == '0' else 'pi3')
        s = ratseries(num, den, N)
        for i in range(N+1): tot[i] += c*s[i]
    return tot

def build(comb, N, f0=Fr(0)):
    return solve_from_psi(psi_series(comb, N), N, f0)

def vp(fr, p):
    if fr == 0: return None
    num, den = fr.numerator, fr.denominator
    v = 0
    while num % p == 0: num//=p; v+=1
    while den % p == 0: den//=p; v-=1
    return v

def measure(name, F, ns, ps=(2,3,5)):
    d = {}
    for n in ns:
        d[n] = math.log(F[n].denominator)/n if F[n] != 0 else 0.0
    sl = {}
    for p in ps:
        vv = [vp(F[n], p) for n in ns]
        sl[p] = vv
    return dict(name=name, logden=d, slopes=sl)

def show(m, ns):
    print(f"{m['name']:44s} logden/n " + " ".join(f"{n}:{m['logden'][n]:.4f}" for n in ns))
    for p, vv in m['slopes'].items():
        print(f"      v_{p} " + " ".join(f"{n}:{v}" for n, v in zip(ns, vv)) +
              f"   slope={(vv[-1]/ns[-1]) if vv[-1] is not None else 0:+.4f}")

if __name__ == '__main__':
    N = 600
    ns = [100, 200, 400, 600]
    cands = {
      'F1 = h(S)  [= H]':                       {(1,'0'): Fr(1)},
      'F3 = h(3S) - 3*hpair(S,pi/3)':           {(3,'0'): Fr(1), (1,'p3'): Fr(-3)},
      'F9 = h(9S)-3*hpair(3S)+9*hpair(S)':      {(9,'0'): Fr(1), (3,'p3'): Fr(-3), (1,'p3'): Fr(9)},
      'K  = h(3S)+3h(S)   [note sec.13]':       {(3,'0'): Fr(1), (1,'0'): Fr(3)},
      'full cubic trace h(3S)+3h(S)-3hpair':    {(3,'0'): Fr(1), (1,'0'): Fr(3), (1,'p3'): Fr(-3)},
      'hpair(S,pi/3) alone':                    {(1,'p3'): Fr(1)},
      'h(2S)':                                  {(2,'0'): Fr(1)},
    }
    series = {}
    for k, c in cands.items():
        F = build(c, N)
        series[k] = F
        show(measure(k, F, ns), ns)
    # linear dependence check: is F3 proportional to F1 mod constants?
    F1 = series['F1 = h(S)  [= H]']; F3 = series['F3 = h(3S) - 3*hpair(S,pi/3)']
    F9 = series['F9 = h(9S)-3*hpair(3S)+9*hpair(S)']
    print("F3/F1 ratios n=1..5:", [str(F3[n]/F1[n]) for n in range(1,6)])
    print("F9/F1 ratios n=1..5:", [str(F9[n]/F1[n]) for n in range(1,6)])
    T = series['full cubic trace h(3S)+3h(S)-3hpair']
    print("full trace coefficients n=1..6 (should vanish):", [str(T[n]) for n in range(1,7)])
