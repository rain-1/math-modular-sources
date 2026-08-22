"""Efficient exact generator for the chi_{-3} well-poised hypergeometric rows.

R_m(t) = (2t+b+1) * prod_{j=1..a}(t-j+1)(t+b+j)
         / prod_{j=0..b} (t+j+1/3)^2 (t+j+2/3)^2

sum_{t>=0} R_m(t) = Q L(2,chi_{-3}) - P,  Q,P in Q.
Everything is done in the variable u = 3t, so the poles sit at the integers
u = -(3j+1), -(3j+2); all quantities are rationals with 3-power denominators
times ordinary integers.
"""
from fractions import Fraction as F

def row(m, a, b=None, want_P=True):
    if b is None: b = m
    poles = [F(-(3*j+1),3) for j in range(b+1)] + [F(-(3*j+2),3) for j in range(b+1)]
    n = len(poles)
    A = [None]*n; B = [None]*n
    for i, p in enumerate(poles):
        # value and log-derivative of N at p, in product form
        val = 2*p + b + 1
        dlog = F(2, 2*p+b+1)
        for j in range(1, a+1):
            f1 = p - j + 1; f2 = p + b + j
            val *= f1*f2
            dlog += F(1,f1) + F(1,f2)
        # 1/prod_{i'!=i}(p-p')^2 and its log-derivative
        pr = F(1); s = F(0)
        for i2, q in enumerate(poles):
            if i2 != i:
                d = p - q
                pr *= d*d
                s += F(1,d)
        Ai = val/pr
        A[i] = Ai
        B[i] = Ai*(dlog - 2*s)
    assert sum(B) == 0
    SA1 = sum(A[i] for i in range(b+1))
    SA2 = sum(A[i] for i in range(b+1, n))
    SB1 = sum(B[i] for i in range(b+1))
    assert SA1 + SA2 == 0, "zeta(2) coefficient not zero"
    assert SB1 == 0, "pi/sqrt3 coefficient not zero"
    Q = 9*SA1
    if not want_P:
        return Q, None
    # cumulative tails: tail2[j] = sum_{l<j} 1/(l+f)^2 , tail1[j] = sum_{l<j} 1/(l+f)
    P = F(0)
    for (off, base) in ((0, F(1,3)), (b+1, F(2,3))):
        t1 = F(0); t2 = F(0)
        for j in range(b+1):
            i = off + j
            P += A[i]*t2 + B[i]*t1
            t2 += F(1,(j+base)**2); t1 += F(1,j+base)
    return Q, P            # sum_t R(t) = Q*L - P
