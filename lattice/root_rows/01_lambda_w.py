#!/usr/bin/env python3
"""Theorem R1: minimal lambda_w with lambda^n [x^n] G^{1/w} in Z for all G in 1+xZ[[x]].

Claim: lambda_w = w * rad(w) = prod_{p|w} p^{v_p(w)+1}.
Graded: lambda = prod_{p|w} p^{max(0, v_p(w)+1-e_p)},  e_p = min_{n>=1} v_p([x^n]G).

Numerical verification of the two ingredients:
  (a) for p | w, v_p(binom(1/w,k)) = -k v_p(w) - (k - s_p(k))/(p-1);
  (b) for p not dividing w, v_p(binom(1/w,k)) >= 0.
plus a direct series test on random G.
"""
from sympy import Rational, binomial, factorint, primefactors, Integer
import random

def vp(x, p):
    x = Rational(x)
    if x == 0: return None
    n, d = x.p, x.q
    v = 0
    while n % p == 0: n//=p; v+=1
    while d % p == 0: d//=p; v-=1
    return v

def sdig(k,p):
    s=0
    while k: s+=k%p; k//=p
    return s

print("=== (a),(b): exact 2-adic/p-adic valuation of binom(1/w,k) ===")
ok=True
for w in range(2,13):
    for p in primefactors(w)+[q for q in [2,3,5,7,11,13] if w%q]:
        a = vp(w,p) if w%p==0 else 0
        for k in range(1,60):
            v = vp(binomial(Rational(1,w),k), p)
            if w%p==0:
                pred = -a*k - (k-sdig(k,p))//(p-1)
                assert (k-sdig(k,p))%(p-1)==0
                if v!=pred: ok=False; print("MISMATCH",w,p,k,v,pred)
            else:
                if v<0: ok=False; print("NEG",w,p,k,v)
print("all valuations match / nonneg:", ok)

print()
print("=== minimal lambda_w ===")
for w in range(2,17):
    lam = 1
    for p,a in factorint(w).items(): lam *= p**(a+1)
    # verify sufficiency: lam^k * binom(1/w,k) in Z for k<=200
    suff = all(vp(binomial(Rational(1,w),k)*Integer(lam)**k, p) >= 0
               for k in range(1,120) for p in primefactors(w))
    # verify minimality: for each p|w, lam/p fails
    minim = True
    for p in primefactors(w):
        lam2 = lam//p
        if all(vp(binomial(Rational(1,w),k)*Integer(lam2)**k, p) >= 0 for k in range(1,400)):
            minim = False
    print(f"w={w:3d}  lambda_w = w*rad(w) = {lam:6d}   sufficient={suff}  minimal={minim}")

print()
print("=== direct series test: lambda^n [x^n] G^{1/w} in Z for random integral G ===")
N=40
random.seed(1)
for w in [2,3,4,5,6,8,9,12]:
    lam=1
    for p,a in factorint(w).items(): lam*=p**(a+1)
    bad=False; badsmall=False
    for trial in range(6):
        S=[0]+[random.randint(-9,9) for _ in range(N)]
        # G^{1/w} by Newton recursion on coefficients
        g=[Rational(1)]+[Rational(0)]*N
        for n in range(1,N+1):
            # (g^w)_n = S_n  -> w*g_n + (terms from lower) = S_n
            # compute conv of g^w up to n with g_n unknown: use  w*g_n*g_0^{w-1}
            gp=[Rational(1)]+[Rational(0)]*N   # g^w with g_n set to 0
            gg=list(g); gg[n]=Rational(0)
            # power
            cur=[Rational(1)]+[Rational(0)]*N
            for _ in range(w):
                new=[Rational(0)]*(N+1)
                for i in range(n+1):
                    if cur[i]==0: continue
                    for j in range(n+1-i):
                        new[i+j]+=cur[i]*gg[j]
                cur=new
            g[n]=(Rational(S[n])-cur[n])/w
        for n in range(N+1):
            if (Integer(lam)**n*g[n]).q!=1: bad=True
        for p in primefactors(w):
            lam2=lam//p
            for n in range(N+1):
                if (Integer(lam2)**n*g[n]).q!=1: badsmall=True
    print(f"w={w:3d} lambda={lam:6d}: integral={not bad}   smaller-lambda-fails={badsmall}")
