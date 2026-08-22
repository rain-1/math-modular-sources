#!/usr/bin/env python3
"""Per-class summary: hit counts, family structure, and the deep analysis of the
small (non-family) classes."""
import os, glob, sys, json, math
from fractions import Fraction as F
from math import gcd
from collections import defaultdict

def dnl(N):
    D=[1]*(N+2)
    for n in range(1,N+2): D[n]=D[n-1]*n//gcd(D[n-1],n)
    return D

def seq(A,be,B,M,j1,j2,C,N,u0,u1):
    de=C*M*M; ep=-C*M*(j1+j2); ze=C*j1*j2
    u=[F(u0),F(u1)]
    for n in range(1,N):
        P=A*n*n+be*n+B; Qn=de*n*n+ep*n+ze
        u.append((P*u[n]-Qn*u[n-1])/F((n+1)**2))
    return u

def main():
    data=defaultdict(list)
    for p in sorted(glob.glob('out/c_*.txt'))+['out/a0.txt']:
        if not os.path.exists(p): continue
        for line in open(p):
            f=line.split()
            if len(f)!=6: continue
            M,j1,j2,A,B,C=map(int,f)
            data[(M,j1,j2)].append((A,B,C))
    print("%-12s %-14s %8s %8s   %s"%("class","(rho;delta)","hits","distinct(A,B)","structure"))
    summary={}
    for cls in sorted(data):
        M,j1,j2=cls
        rho=F(j1+j2,2*M); dl=F(abs(j2-j1),M)
        hits=data[cls]
        ab=defaultdict(int)
        for (A,B,C) in hits: ab[(A,B)]+=1
        struct = "FAMILY (C free)" if max(ab.values())>20 else "sporadic"
        print("(%d;%d,%d)".ljust(12)%cls + " (%s;%s)".ljust(15)%(rho,dl)
              + " %8d %8d   %s"%(len(hits),len(ab),struct))
        summary[str(cls)]=dict(rho=str(rho),delta=str(dl),hits=len(hits),ab=len(ab),struct=struct)
        if struct=="FAMILY":
            tops=sorted(ab.items(), key=lambda kv:-kv[1])[:5]
            print("      top (A,B) by count:", tops)
    json.dump(summary,open('out/classsummary.json','w'),indent=1)
main()
