"""Calibration of the CDT holonomy-bound arithmetic on the two examples of
arXiv:2408.15403 that we will imitate for Catalan:

  (i)  Theorem mainA  : 1, zeta(2), L(2,chi_-3)   (their Sec. 13)   -- m=14
  (ii) Theorem logsmain: products of two logarithms (their Sec. 9)  -- m=17
  (iii) their Sec. 7 remark: the SAME 17 functions read on P^1-{0,1,oo} with
        the integrated type n[1..n]^2 (the "unsymmetrised" comparison).

Run:  python3 calib.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', 'cdt_finder'))
from fractions import Fraction
import math
from cdt_bound import tau_flat, tau_sharp, tau, I_uvw

def show(name, m, cols, e, exp_flat=None, exp_sharp=None, exp_tot=None):
    sm, tf = tau_flat(m, cols)
    ts, xi = tau_sharp(m, e)
    tot = float(tf) + ts
    print(f'--- {name}:  m={m}, cols={cols}, sum(e)={sum(e)}, max(e)={max(e)}')
    print(f'    sigma_m    = {sm}')
    print(f'    tau^flat   = {tf} = {float(tf):.6f}' +
          (f'   (CDT: {exp_flat})' if exp_flat else ''))
    print(f'    tau^sharp  = {ts:.8f}  at xi={xi:.4f}' +
          (f'   (CDT: {exp_sharp} = {float(Fraction(exp_sharp)):.8f})' if exp_sharp else ''))
    print(f'    tau        = {tot:.8f}' +
          (f'   (CDT: {exp_tot} = {float(Fraction(exp_tot)):.8f})' if exp_tot else ''))
    return tot

print('=== (i) CDT Theorem A, 1/zeta(2)/L(2,chi_-3), Y_0(2) picture ===')
t1 = show('mainA', 14, [(1,2),(3,2)],
          [0,0,1, 0,0,0,0,0,0, 1,1,1,1,1],
          '191/49', '27/80', '16603/3920')
lp = math.log(256*5448339453535586608000000000/8658833407565631122430056127)
print(f'    log|phi\'(0)| = {lp:.6f};  entry = {lp-t1:+.6f};  BC = 11.845;'
      f'  bound m <= {11.845/(lp-t1):.4f}  (CDT 13.9938);'
      f'  margin {14*(lp-t1)-11.845:+.4f}')

print()
print('=== (ii) CDT Theorem logsmain, 17 functions, y-line, type n[1..2n]^2 ===')
t2 = show('logs (y-line)', 17, [(1,2),(3,2)],
          [0,0,1, 0,0,0,0,0,0, 1,1,1,1,1,1,1,1],
          '1136/289', '78419/242760', '1032659/242760')
print('    CDT quote the number 1032659/242760 = 4.2538... and m <= 16.2.')

print()
print('=== (iii) CDT Sec.7 remark: 17 functions on P^1-{0,1,oo}, type n[1..n]^2 ===')
# b array: columns (u_1,b_1)=(2,1), (u_2,b_2)=(4,1)   [their b array, entries 0/1]
t3 = show('logs (x-line)', 17, [(2,1),(4,1)],
          [0,1,1,1,1,1,1,1, 0, 0,0,0,0,1,1,1,1],
          '558/289', '83711/242760', '552431/242760')
Gp = 0.9163768
lp3 = math.log(16*Gp)
print(f'    |G\'(0)| = {Gp};  log|phi\'(0)| = log16 + log|G\'(0)| = {lp3:.6f}')
print(f'    entry = {lp3-t3:+.6f};  CDT quotient 22.7527  =>  BC = {22.7527*(lp3-t3):.4f}')
print(f'    (CDT: "the quotient comes out to approximately 22.7527, a long distance from 17")')
