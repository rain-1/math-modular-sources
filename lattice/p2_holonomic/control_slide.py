#!/usr/bin/env python3
"""lattice/p2_holonomic/control_slide.py
Control for the E-slide of run_h7/run_h8.  Sliding the 2-adic exponent E by 1
does exactly this to the Hermite data: h22 is unchanged, h11 doubles, and h12
either stays or gains h11 -- one fresh bit.  So theta_j = h12/h11 follows
    theta_{j+1} = theta_j/2  or  theta_j/2 + 1/2,
i.e. theta_j is the *reversed* binary expansion of a random real, and the
balance target R_j doubles.  This script runs that model with fair coins and
measures the persistence of the parity of the balance index, so that the
Catalan persistence can be compared against it."""
import random, sys
from math import log2

def cf_ladder(num, den):
    """partial quotients and convergents of num/den (0<=num<den)"""
    a, b = num, den
    ps, qs = [0, 1], [1, 0]          # p_{-1},p_0 ; q_{-1},q_0  (theta<1)
    A, B = [], []
    p0, p1 = 0, 1; q0, q1 = 1, 0     # convergent i-1, i-2 style
    # standard: theta = [0;a1,a2,...]; convergents p_i/q_i, p_0/q_0 = 0/1
    P = [0]; Q = [1]
    x, y = num, den
    # continued fraction of x/y
    quots = []
    while y:
        q = x // y
        quots.append(q); x, y = y, x - q*y
    # quots[0] = 0 since num<den
    p_m1, q_m1 = 1, 0
    p_0, q_0 = quots[0], 1
    P = [p_0]; Q = [q_0]
    pp, qq = p_m1, q_m1
    pc, qc = p_0, q_0
    for a in quots[1:]:
        pn, qn = a*pc + pp, a*qc + qq
        pp, qq = pc, qc; pc, qc = pn, qn
        P.append(pc); Q.append(qc)
    A = [abs(num*Q[i] - den*P[i]) for i in range(len(Q))]
    return A, Q, quots

def balance_idx(A, B, tnum, tden):
    """argmin_i A_i^2 + (tnum/tden) B_i^2 with exact integers"""
    best = 0; bv = A[0]*A[0]*tden + tnum*B[0]*B[0]
    for i in range(1, len(A)):
        v = A[i]*A[i]*tden + tnum*B[i]*B[i]
        if v < bv: bv = v; best = i
    return best

def run(nchain=400, steps=120, j0=380, seed=20260823):
    rng = random.Random(seed)
    agree = [0]*13; cnt = [0]*13
    evens = 0; tot = 0
    for c in range(nchain):
        h12 = rng.getrandbits(j0) | 1
        pars = []
        for j in range(steps+1):
            den = 1 << (j0 + j)
            num = h12
            g = 1
            while num % 2 == 0 and den % 2 == 0:   # keep it reduced
                num //= 2; den //= 2
            A, B, quots = cf_ladder(num, den)
            i0 = balance_idx(A, B, 1, 1)
            pars.append(i0 & 1)
            h12 = h12 + (rng.getrandbits(1) << (j0 + j))
        evens += sum(1 for p in pars if p == 0); tot += len(pars)
        for lag in range(1, 13):
            agree[lag] += sum(1 for a, b in zip(pars, pars[lag:]) if a == b)
            cnt[lag]   += len(pars) - lag
    print(f"control: {nchain} chains x {steps} steps, theta of bit-length {j0}")
    print(f"  even fraction {evens/tot:.4f}")
    print("  parity agreement at lag: " +
          " ".join(f"{l}:{agree[l]/cnt[l]:.3f}" for l in range(1, 13)))

if __name__ == "__main__":
    run(nchain=int(sys.argv[1]) if len(sys.argv) > 1 else 300,
        steps=int(sys.argv[2]) if len(sys.argv) > 2 else 60,
        j0=int(sys.argv[3]) if len(sys.argv) > 3 else 400)
