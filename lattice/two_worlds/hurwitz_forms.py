from fractions import Fraction as F
from math import factorial

def Jform(m, s):
    """J_m(s) = int int [x(1-x)y(1-y)]^m (xy)^{s-1}/(1-xy)^{m+1}
       = CA * zeta(2,s) - CD ,  CA,CD in Q.
       R(k) = m! prod_{j=1..m}(k+j) / prod_{i=m..2m}(k+s+i)^2 """
    s = F(s)
    poles = list(range(m, 2*m+1))          # i
    A = []; B = []
    for i in poles:
        k0 = -(s+i)
        pr = 1
        for i2 in poles:
            if i2 != i: pr *= (i2-i)**2
        num = factorial(m)
        for j in range(1, m+1): num *= (k0+j)
        f = F(num, 1)/pr
        df = sum(F(1)/(k0+j) for j in range(1, m+1)) - 2*sum(F(1)/(i2-i) for i2 in poles if i2 != i)
        A.append(f); B.append(f*df)
    assert sum(B) == 0, sum(B)
    CA = sum(A)
    CD = sum(A[t]*sum(F(1)/(l+s)**2 for l in range(poles[t])) for t in range(len(poles))) \
       + sum(B[t]*sum(F(1)/(l+s)    for l in range(poles[t])) for t in range(len(poles)))
    return CA, CD

for m in range(0, 8):
    c1, d1 = Jform(m, F(1,3))
    c2, d2 = Jform(m, F(2,3))
    print("m=%d  C(1/3)=%s  C(2/3)=%s  ratio=%s" % (m, c1, c2, c1/c2 if c2 else None))
