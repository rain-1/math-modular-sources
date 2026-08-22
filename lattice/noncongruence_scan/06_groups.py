#!/usr/bin/env python3
"""Complete classification of the genus-zero subgroups Gamma <= PSL_2(Z) of index
<= 12 that have EXACTLY FOUR special points (cusps + elliptic points).

Why index <= 12: a 3-term recurrence means the second-order operator has exactly
four singular points {0,t_1,t_2,inf} on the t-line, i.e. Gamma has exactly four
special points.  Gauss-Bonnet with g=0, r elliptic points of orders e_i and
s cusps, r+s=4:
      area/(pi/3) = 6(2g-2+s+sum(1-1/e_i)) = 6(2 - sum_i 1/e_i)  <= 12,
with equality iff r=0 (four cusps: the Beauville/Zagier case).

Subgroups of index n of PSL_2(Z) = <S,T | S^2 = T^3 = 1> correspond to transitive
pairs (s,t) in S_n with s^2 = t^3 = 1, up to simultaneous conjugacy.  Then
  e_2 = #fix(s), e_3 = #fix(t), cusps = cycles of st,
  g = 1 + n/12 - e_2/4 - e_3/3 - c/2.
Congruence test (Wohlfahrt): N = lcm of cusp widths; Gamma is congruence iff the
coset action factors through PSL_2(Z/N).
"""
import sys, itertools
from sympy.combinatorics import Permutation, PermutationGroup
from math import gcd
from functools import reduce

def lcm(a,b): return a*b//gcd(a,b)

def involutions(n):
    """all s in S_n with s^2=1, as tuples"""
    res=[]
    def rec(pos, used, cur):
        if pos==n: res.append(tuple(cur)); return
        if pos in used: rec(pos+1, used, cur); return
        # fixed
        cur[pos]=pos; used.add(pos); rec(pos+1, used, cur); used.discard(pos)
        for j in range(pos+1, n):
            if j in used: continue
            cur[pos]=j; cur[j]=pos; used.add(pos); used.add(j)
            rec(pos+1, used, cur)
            used.discard(pos); used.discard(j)
    rec(0,set(),[0]*n)
    return res

def order3(n):
    res=[]
    def rec(pos, used, cur):
        if pos==n: res.append(tuple(cur)); return
        if pos in used: rec(pos+1, used, cur); return
        cur[pos]=pos; used.add(pos); rec(pos+1, used, cur); used.discard(pos)
        rest=[j for j in range(pos+1,n) if j not in used]
        for a,b in itertools.permutations(rest,2):
            cur[pos]=a; cur[a]=b; cur[b]=pos
            used.update((pos,a,b)); rec(pos+1, used, cur)
            used.difference_update((pos,a,b))
    rec(0,set(),[0]*n)
    return res

def cycles(p):
    n=len(p); seen=[False]*n; out=[]
    for i in range(n):
        if seen[i]: continue
        c=[]; j=i
        while not seen[j]: seen[j]=True; c.append(j); j=p[j]
        out.append(c)
    return out

def transitive(s,t):
    n=len(s); seen={0}; st=[0]
    while st:
        i=st.pop()
        for p in (s,t):
            if p[i] not in seen: seen.add(p[i]); st.append(p[i])
    return len(seen)==n

# ---------- PSL_2(Z/N) as a permutation group on P^1(Z/N) ------------------
def p1(N):
    pts=[]
    for c in range(N):
        for d in range(N):
            if gcd(gcd(c,d),N)==1: pts.append((c%N,d%N))
    # normalise up to units
    units=[u for u in range(N) if gcd(u,N)==1]
    reps={}; out=[]
    for p in pts:
        key=min(tuple(((u*p[0])%N,(u*p[1])%N)) for u in units)
        if key not in reps: reps[key]=len(out); out.append(key)
    return out, reps, units

def psl_perms(N):
    pts,reps,units=p1(N)
    def act(M,p):
        c,d=p
        return (( M[0][0]*c + M[1][0]*d )%N, ( M[0][1]*c + M[1][1]*d )%N)
    def perm(M):
        img=[]
        for p in pts:
            q=act(M,p)
            key=min(tuple(((u*q[0])%N,(u*q[1])%N)) for u in units)
            img.append(reps[key])
        return img
    S=[[0,-1],[1,0]]; T=[[1,1],[0,1]]
    return perm(S), perm(T), len(pts)

def is_congruence(s,t,N):
    """does the coset action factor through PSL_2(Z/N)?"""
    if N==1: return True
    Sp,Tp,m=psl_perms(N)
    n=len(s)
    G1=PermutationGroup([Permutation(Sp),Permutation(Tp)])
    # graph subgroup inside S_{m+n}
    def joint(a,b):
        return Permutation(list(a)+[m+x for x in b])
    G2=PermutationGroup([joint(Sp,s), joint(Tp,t)])
    return G2.order()==G1.order()

def run(nmax=12):
    print("n  e2 e3 cusps  widths            genus  signature            congruence  Wohlfahrt")
    found=[]
    for n in range(1,nmax+1):
        invs=involutions(n); o3=order3(n)
        seen=set()
        for s in invs:
            for t in o3:
                if not transitive(s,t): continue
                e2=sum(1 for i in range(n) if s[i]==i)
                e3=sum(1 for i in range(n) if t[i]==i)
                st=[s[t[i]] for i in range(n)]
                cy=cycles(st); c=len(cy)
                if 12*(1) + n - 3*e2 - 4*e3 - 6*c != 0:   # g=0  <=> 1+n/12-e2/4-e3/3-c/2 = 0
                    pass
                g2 = 12 + n - 3*e2 - 4*e3 - 6*c          # = 12*g
                if g2 != 0: continue
                if e2+e3+c != 4: continue
                # canonical form up to simultaneous conjugacy
                key=canon(s,t,n)
                if key in seen: continue
                seen.add(key)
                widths=sorted(len(x) for x in cy)
                N=reduce(lcm, widths, 1)
                cong=is_congruence(s,t,N)
                sig="(0;"+",".join(["2"]*e2+["3"]*e3)+";"+str(c)+" cusps)"
                found.append((n,e2,e3,c,widths,sig,cong,N))
                print(f"{n:2d} {e2:3d}{e3:3d}{c:6d}  {str(widths):18s} 0     {sig:20s} {str(cong):10s} {N}")
    print()
    nc=[f for f in found if not f[6]]
    print("NON-CONGRUENCE candidates:", len(nc))
    for f in nc: print("  ", f)
    print("total groups:", len(found))

def canon(s,t,n):
    """canonical label under simultaneous conjugacy (brute force for small n)"""
    best=None
    for perm in itertools.permutations(range(n)):
        inv=[0]*n
        for i,p in enumerate(perm): inv[p]=i
        ss=tuple(perm[s[inv[i]]] for i in range(n))
        tt=tuple(perm[t[inv[i]]] for i in range(n))
        cand=(ss,tt)
        if best is None or cand<best: best=cand
    return best

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 9)
