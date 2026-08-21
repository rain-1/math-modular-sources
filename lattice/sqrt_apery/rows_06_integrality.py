"""
rows_06_integrality.py -- the 2-adic lemma that explains the whole lambda column,
and the sharpness obstruction to a general 2^n version.

RIVER'S LEMMA [proved].  For G in 1 + t Z[[t]],  4^n [t^n] sqrt(G) in Z,
because binom(1/2,k) = (-1)^(k-1) 2 Cat(k-1)/4^k and only k <= n contribute to [t^n].

GRADED FORM [proved, this file].  Suppose 2^e | G_n for every n >= 1, i.e.
G = 1 + 2^e K with K in t Z[[t]].  Then
    sqrt(G) = sum_k binom(1/2,k) 2^(e k) K^k,
    binom(1/2,k) 2^(e k) = (-1)^(k-1) 2 Cat(k-1) 2^((e-2)k),
and K^k in t^k Z[[t]], so for lambda = max(1, 2^(2-e)):
    lambda^n [t^n] sqrt(G) = sum_{k<=n} (-1)^(k-1) 2 Cat(k-1) lambda^(n-k)
                              * (lambda 2^(e-2))^k [t^n] K^k    in Z.
   e = 0 -> lambda = 4      e = 1 -> lambda = 2      e >= 2 -> lambda = 1.
So the lambda of each row is READ OFF from the 2-divisibility of its parent.

SHARPNESS OF THE LEMMA [proved].  No 2^n version can hold for general
G in 1+tZ[[t]]: for G = 1+t, 2^n[t^n]sqrt(1+t) = (-1)^(n-1) 2 Cat(n-1)/2^n and
Cat(n-1) is odd exactly when n is a power of 2, so v_2 = 1-n < 0 there.
"""
from fractions import Fraction as Fr
from math import comb

N = 400

def parent_AZ(a,b,c,N):
    A=[Fr(1),Fr(b)]
    for n in range(1,N): A.append(((2*n+1)*(a*n*n+a*n+b)*A[n]-c*n**3*A[n-1])/Fr((n+1)**3))
    return [int(x) for x in A]

def parent_cooper(name,N):
    if name=='s7':  A=[Fr(1),Fr(4)]; f=lambda n,u1,u0:(2*n+1)*(13*n*n+13*n+4)*u1+3*n*(9*n*n-1)*u0
    if name=='s10': A=[Fr(1),Fr(2)]; f=lambda n,u1,u0:2*(2*n+1)*(3*n*n+3*n+1)*u1+4*n*(16*n*n-1)*u0
    if name=='s18': A=[Fr(1),Fr(6)]; f=lambda n,u1,u0:2*(2*n+1)*(7*n*n+7*n+3)*u1-12*n*(16*n*n-1)*u0
    for n in range(1,N): A.append(f(n,A[n],A[n-1])/Fr((n+1)**3))
    return [int(x) for x in A]

def v2(m):
    if m==0: return 10**9
    v=0
    while m%2==0: m//=2; v+=1
    return v

def sqrt_series(A):
    c=[Fr(1)]
    for n in range(1,len(A)): c.append((Fr(A[n])-sum(c[k]*c[n-k] for k in range(1,n)))/2)
    return c

ROWS = [("Domb",        parent_AZ(10,4,64,N),  1),
        ("T",           parent_AZ(12,4,16,N),  1),
        ("AZ(9,3,-27)", parent_AZ(9,3,-27,N),  4),
        ("AZ(11,5,125)",parent_AZ(11,5,125,N), 4),
        ("AZ(7,3,81)",  parent_AZ(7,3,81,N),   4),
        ("Cooper s7",   parent_cooper('s7',N), 1),
        ("Cooper s10",  parent_cooper('s10',N),2),
        ("Cooper s18",  parent_cooper('s18',N),2)]

print(f"{'row':14s} {'e = min v_2(A_n), n>=1':24s} {'lambda predicted':17s} {'lambda used':11s} match")
for nm, A, lam in ROWS:
    e = min(v2(A[n]) for n in range(1,N+1))
    e = min(e, 2)                      # only e<=2 matters
    pred = max(1, 2**(2-e))
    print(f"{nm:14s} {e:<24d} {pred:<17d} {lam:<11d} {pred==lam}")

print()
print("sharp 2-adic profile of the raw square root, min_n ( n - v_2(den [t^n]sqrt F) ):")
for nm, A, lam in ROWS:
    c = sqrt_series(A)
    prof = []
    for n in range(N+1):
        d = c[n].denominator; v=0
        while d%2==0: d//=2; v+=1
        prof.append(n-v)
    # for lambda=1 rows the denominators should be trivial
    dens = max((c[n].denominator for n in range(N+1)))
    print(f"  {nm:14s} min(n - v_2(den)) = {min(prof):3d}   max denominator of [t^n]sqrt F = "
          f"{'1 (already integral)' if dens==1 else '2^'+str(max(n-prof[n] for n in range(N+1)))}")

print()
print("obstruction to a general 2^n lemma (G = 1+t):")
for k in range(1,9):
    n=2**k; cat=comb(2*(n-1),n-1)//n
    print(f"   n=2^{k}={n:5d}:  v_2( 2^n [t^n] sqrt(1+t) ) = {1+v2(cat)-n}")
