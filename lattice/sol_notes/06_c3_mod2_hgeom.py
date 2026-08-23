from itertools import product
from fractions import Fraction as F
import sympy as sp

mz0 = sp.Matrix([[-1,2,-1,-1],[1,-4,2,2],[1,8,-4,-3],[-6,-19,9,7]])
mzi = sp.Matrix([[-8,-5,2,1],[14,9,-4,-2],[-22,-15,5,3],[45,32,-12,-7]])
mz1 = sp.Matrix([[2,1,0,0],[-3,-2,0,0],[0,0,1,0],[0,0,-3,-1]])
jz  = sp.Matrix([[0,0,0,1],[0,0,1,1],[0,-1,0,0],[-1,-1,0,0]])

print("=== Sp4(F2) = S6 action on the 6 odd theta characteristics ===")
J2 = [[int(jz[i,j])%2 for j in range(4)] for i in range(4)]
vecs = list(product([0,1],repeat=4))
def qform(d, x):
    s = 0
    for i in range(4):
        for j in range(i+1,4):
            s ^= (J2[i][j] & x[i] & x[j])
    for i in range(4):
        s ^= (d[i] & x[i])
    return s
# sanity: q(x+y)+q(x)+q(y) = <x,y>
def sympair(x,y):
    s=0
    for i in range(4):
        for j in range(4):
            s ^= (J2[i][j] & x[i] & y[j])
    return s
ok=True
for d in vecs:
    for x in vecs:
        for y in vecs:
            xy = tuple((x[i]^y[i]) for i in range(4))
            if (qform(d,xy)^qform(d,x)^qform(d,y)) != sympair(x,y): ok=False
print("quadratic refinement axiom holds for all 16 d:", ok)
odd = []
for d in vecs:
    z = sum(1 for x in vecs if qform(d,x)==0)
    if z==6: odd.append(d)
    assert z in (6,10), z
print("number of odd (Arf=1) theta characteristics:", len(odd))

def act(mat, d):
    # (M.q)(x) = q(M^{-1} x)
    Minv = [[int(x)%2 for x in row] for row in (mat.inv()%2).tolist()]
    # find d' with q_{d'}(x) = q_d(Minv x) for all x
    for dd in vecs:
        good = True
        for x in vecs:
            mx = tuple(sum(Minv[i][j]*x[j] for j in range(4))%2 for i in range(4))
            if qform(dd,x) != qform(d,mx): good=False; break
        if good: return dd
    raise RuntimeError("no image")

def cycletype(mat):
    perm = {i: odd.index(act(mat, odd[i])) for i in range(6)}
    seen=set(); ct=[]
    for i in range(6):
        if i in seen: continue
        l=0; j=i
        while j not in seen:
            seen.add(j); j=perm[j]; l+=1
        ct.append(l)
    return tuple(sorted(ct, reverse=True))

for nm, mm in [("mz0 (T_0)",mz0),("mz1 (T_1)",mz1),("mzi (T_inf)",mzi),
               ("mzi^-1",mzi.inv()),("mz0*mz1",mz0*mz1)]:
    print(f"  {nm:12s} order mod 2 -> cycle type on 6 odd thetas: {cycletype(mm)}")
print("note (14.1) claims local cycle types (6), (2,2,2), (3,1,1,1)")

print()
print("=== hypergeometric alpha/beta vs 5F4(1/2,1/6,1/6,5/6,5/6 ; 1,1,2/3,4/3) ===")
a = [F(1,2),F(1,6),F(1,6),F(5,6),F(5,6)]
b = [F(1),F(1),F(2,3),F(4,3),F(1)]     # 4 stated denominators + implicit b_5 = 1
exp0 = sorted(( (1-bj) % 1 for bj in b ))
expinf = sorted(( ai % 1 for ai in a ))
print("numerator params a_i           =", a)
print("denominator params b_j (+1)    =", b)
print("exponents at 0  (= 1-b_j mod1) =", exp0)
print("exponents at oo (= a_i mod 1)  =", expinf)
alpha = sorted([F(0),F(0),F(0),F(1,3),F(2,3)])
beta  = sorted([F(1,6),F(1,6),F(1,2),F(5,6),F(5,6)])
print("note's alpha                   =", alpha, " matches exponents at 0 :", alpha==exp0)
print("note's beta                    =", beta,  " matches exponents at oo:", beta==expinf)
print("alpha, beta disjoint (irreducible):", set(alpha).isdisjoint(set(beta)))
print("sum alpha =", sum(alpha), " sum beta =", sum(beta), " -> T_1 nontrivial eigenvalue exp(2pi i (sum a - sum b)) =",
      "exp(2 pi i *", (sum(beta)-sum(alpha))%1, ") =", "-1" if (sum(beta)-sum(alpha))%1==F(1,2) else "?")

# Hodge numbers by the interlacing / Fedorov rule
h={}
for k,bk in enumerate(beta):
    p = sum(1 for aj in alpha if aj < bk) - k
    h[p]=h.get(p,0)+1
print("Hodge numbers of V by interlacing rule:", [h.get(p,0) for p in sorted(h, reverse=True)], " (note claims (1,3,1))")

print()
print("=== is V = Lambda^2_0 W compatible with the local spectra of section 2? ===")
def lam20(exps):
    s = sorted([ (exps[i]+exps[j])%1 for i in range(len(exps)) for j in range(i+1,len(exps)) ])
    # remove one 0 (the symplectic form)
    s.remove(F(0))
    return sorted(s)
W0   = [F(1,3),F(1,3),F(2,3),F(2,3)]        # omega J2, omega^-1 J2
W1   = [F(1,2),F(1,2),F(0),F(0)]            # (-1,-1,1,1)
Winf = [F(1,2),F(1,2),F(1,6),F(5,6)]        # (-J2, zeta6, zeta6^-1)
for nm, e, target in [("0",W0,alpha),("1",W1,None),("oo",Winf,beta)]:
    L = lam20(e)
    print(f"  Lambda^2_0 W exponents at {nm:2s}: {L}" + (f"   target {target}   match: {L==target}" if target else "   (V's T_1 must be a reflection (0,0,0,0,1/2))"))
print()
print("  KEY: Lambda^2_0 of ANY rank-4 symplectic system always has eigenvalue 1 (exponent 0)")
print("       in every local monodromy, since eigenvalues are {t,1/t,s,1/s} and t*(1/t)=1.")
print("       beta =", beta, "contains no 0  =>  V is NOT Lambda^2_0 of any symplectic rank-4 system.")
print()
print("  Twist fix: let chi be the rank-1 system with monodromy (1,-1,-1) at (0,1,oo):")
for nm, e, target in [("0",W0,alpha),("oo",Winf,beta)]:
    sh = F(0) if nm=="0" else F(1,2)
    L = sorted([(x+sh)%1 for x in lam20(e)])
    print(f"   (Lambda^2_0 W (x) chi) at {nm:2s}: {L}   target {target}   match: {L==target}")
L1 = sorted([(x+F(1,2))%1 for x in lam20(W1)])
print(f"   (Lambda^2_0 W (x) chi) at 1 : {L1}   = reflection (0,0,0,0,1/2)? {L1==sorted([F(0)]*4+[F(1,2)])}")
