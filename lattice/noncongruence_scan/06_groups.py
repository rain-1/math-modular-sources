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

def _unused_involutions(n):
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

def _unused_order3(n):
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
    # PSL_2(Z) = < S, L | S^2 = L^3 = 1 >,  S = [[0,-1],[1,0]], L = ST = [[0,-1],[1,1]]
    S=[[0,-1],[1,0]]; L=[[0,-1],[1,1]]
    return perm(S), perm(L), len(pts)

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


def gen_t(n, e3):
    """all t in S_n with t^3=1 and exactly e3 fixed points (n-e3 divisible by 3)"""
    res=[]
    def rec(rem, cur, nfix):
        if not rem:
            if nfix==e3: res.append(tuple(cur))
            return
        i=rem[0]; rest=rem[1:]
        if nfix<e3:
            cur[i]=i; rec(rest, cur, nfix+1)
        for a_i in range(len(rest)):
            for b_i in range(len(rest)):
                if a_i==b_i: continue
                a,b=rest[a_i],rest[b_i]
                cur[i]=a; cur[a]=b; cur[b]=i
                rec([x for x in rest if x not in (a,b)], cur, nfix)
    rec(list(range(n)), [0]*n, 0)
    return res

def std_s(n, e2):
    """standard involution with e2 fixed points"""
    s=list(range(n)); k=(n-e2)//2
    for i in range(k): s[2*i],s[2*i+1]=2*i+1,2*i
    return tuple(s)

def canon(s,t,n):
    """BFS canonical form: relabel by discovery order from each base point"""
    best=None
    for b in range(n):
        lab={b:0}; order=[b]; qi=0
        while qi<len(order):
            x=order[qi]; qi+=1
            for y in (t[x], s[x]):
                if y not in lab: lab[y]=len(order); order.append(y)
        if len(lab)!=n: continue
        ss=[0]*n; tt=[0]*n
        for x in range(n): ss[lab[x]]=lab[s[x]]; tt[lab[x]]=lab[t[x]]
        cand=(tuple(ss),tuple(tt))
        if best is None or cand<best: best=cand
    return best

def run(nmax=12):
    sigs=[]
    for e2 in range(0,5):
        for e3 in range(0,5-e2):
            c=4-e2-e3
            if c<1: continue
            n=3*e2+4*e3+6*c-12
            if n<1 or n>nmax: continue
            if (n-e2)%2 or (n-e3)%3: continue
            sigs.append((n,e2,e3,c))
    sigs.sort()
    found=[]
    print("  n  e2 e3  cusps  widths           signature              congruence  Wohlfahrt level")
    for (n,e2,e3,c) in sigs:
        s=std_s(n,e2); seen=set()
        for t in gen_t(n,e3):
            if not transitive(s,t): continue
            st=[s[t[i]] for i in range(n)]
            cy=cycles(st)
            if len(cy)!=c: continue
            key=canon(s,t,n)
            if key in seen: continue
            seen.add(key)
            widths=sorted(len(x) for x in cy)
            N=reduce(lcm, widths, 1)
            cong=is_congruence(s,t,N)
            sig="(0;"+",".join(["2"]*e2+["3"]*e3)+";"+str(c)+")"
            found.append((n,e2,e3,c,widths,sig,cong,N,key))
            print(f"{n:3d} {e2:3d}{e3:3d}{c:6d}  {str(widths):16s} {sig:22s} {str(cong):10s} {N}")
    print()
    nc=[f for f in found if not f[6]]
    print("total genus-0 four-special-point subgroups of PSL_2(Z), index<=%d : %d"%(nmax,len(found)))
    print("NON-CONGRUENCE among them: %d"%len(nc))
    for f in nc: print("   index",f[0],"sig",f[5],"widths",f[4],"Wohlfahrt",f[7])

if __name__=="__main__":
    run(int(sys.argv[1]) if len(sys.argv)>1 else 12)
