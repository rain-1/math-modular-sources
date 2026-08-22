#!/usr/bin/env python3
"""Herfurtner's list of rational elliptic surfaces over P^1 with four singular
fibres (Math. Ann. 291 (1991) 319-342, Table 3), and the cross-ratio matching
test for second-order Apery-like rows.

For a row (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1} the four singular points of
   L = th^2 - t P(th) + t^2 Q(th+1)
are 0 (MUM), t_1, t_2 = the roots of R(t) = 1 - A t + D t^2, and infinity.
The Moebius invariant of the labelled quadruple (0, oo; t_1, t_2) is
   z = t_1/t_2 = lam_2/lam_1,      I := (1+z)^2/z = A^2 / D,
symmetric under z <-> 1/z, i.e. under swapping 0 <-> oo and t_1 <-> t_2.
For a Herfurtner configuration with base points p = (p_a,p_b,p_c,p_d) and an
assignment (p_0, p_oo | p_1, p_2) the same invariant is
   z = ((p_1-p_0)(p_2-p_oo)) / ((p_1-p_oo)(p_2-p_0)),   I = (1+z)^2/z.
A row can be the Picard-Fuchs system of the configuration only if the two
invariants agree.
"""
from sympy import (Rational as Q, sqrt, I as ii, simplify, nsimplify, oo, S,
                   Symbol, exp, pi, radsimp, expand, together, cancel, N as num)
from itertools import permutations

eta = (-S(1) + sqrt(-3))/2          # primitive cube root of unity
INF = oo

# Kodaira Euler numbers
def euler(f):
    if f[0]=='I' and f.endswith('*'):  return int(f[1:-1])+6
    if f[0]=='I' and f[1].isdigit():   return int(f[1:])
    return {'II':2,'III':3,'IV':4,'IV*':8,'III*':9,'II*':10}[f]

# PSL_2(Z) local type: 'c' = parabolic (cusp), 2 or 3 = elliptic order, 1 = trivial
def psl(f):
    if f=='I0*': return 1
    if f[0]=='I' and (f[1].isdigit()):
        return 'c'
    if f.endswith('*') and f[0]=='I': return 'c'
    return {'II':3,'III':2,'IV':3,'IV*':3,'III*':2,'II*':3}[f]

w5p=(1+sqrt(5))/2; w5m=(1-sqrt(5))/2
# ---- Herfurtner Table 3, rows 13-50: the 38 rigid configurations -------------
RIGID = [
 (13,['I1','I1','I1','I9'],  [1,eta,eta**2,INF], 12),
 (14,['I1','I1','I2','I8'],  [-1,1,0,INF], 12),
 (15,['I1','I2','I3','I6'],  [4,Q(-1,2),0,INF], 12),
 (16,['I1','I1','I5','I5'],  [w5p**5, w5m**5, 0, INF], 12),
 (17,['I2','I2','I4','I4'],  [-1,1,0,INF], 12),
 (18,['I3','I3','I3','I3'],  [1,eta,eta**2,INF], 12),
 (19,['I1','I1','I8','II'],  [-(1+ii*sqrt(2))**4/3, -(1-ii*sqrt(2))**4/3, INF, 0], 10),
 (20,['I1','I2','I7','II'],  [Q(-9,4),Q(-8,9),INF,0], 10),
 (21,['I1','I4','I5','II'],  [-10,0,INF,Q(1,8)], 10),
 (22,['I2','I3','I5','II'],  [Q(-5,9),0,INF,3], 10),
 (23,['I1','I1','I7','III'], [((1+ii*sqrt(7))/2)**7/4, ((1-ii*sqrt(7))/2)**7/4, INF, 0], 9),
 (24,['I1','I2','I6','III'], [4,1,INF,0], 9),
 (25,['I1','I3','I5','III'], [Q(-25,3),0,INF,Q(1,5)], 9),
 (26,['I2','I3','I4','III'], [Q(-1,3),0,INF,1], 9),
 (27,['I1','I1','I6','IV'],  [1,-1,INF,0], 8),
 (28,['I1','I2','I5','IV'],  [Q(-27,4),Q(-1,2),INF,0], 8),
 (29,['I3','I3','I2','IV'],  [INF,0,-1,1], 8),
 (30,['I1','I7','II','II'],  [0,INF,-(-1+3*ii*sqrt(3))**2/4, -(-1-3*ii*sqrt(3))**2/4], 8),
 (31,['I2','I6','II','II'],  [0,INF,1,-1], 8),
 (32,['I4','I4','II','II'],  [-1+sqrt(3),-1-sqrt(3),Q(1,2),-4], 8),
 (33,['I1','I6','II','III'], [-(3*ii*sqrt(3)+1)*Q(2,7), INF, 1, 0], 7),
 (34,['I2','I5','II','III'], [Q(125,14),INF,0,Q(27,2)], 7),
 (35,['I3','I4','II','III'], [0,INF,-27,1], 7),
 (36,['I1','I5','II','IV'],  [Q(-16,3),INF,3,0], 6),
 (37,['I2','I4','II','IV'],  [Q(1,9),INF,1,0], 6),
 (38,['I1','I5','III','III'],[Q(-11,2),INF,ii,-ii], 6),
 (39,['I2','I4','III','III'],[0,INF,1,-1], 6),
 (40,['I3','I3','III','III'],[3+2*sqrt(3),3-2*sqrt(3),0,INF], 6),
 (41,['I1','I4','III','IV'], [Q(-27,5),INF,1,0], 5),
 (42,['I2','I3','III','IV'], [Q(1,5),INF,1,0], 5),
 (43,['I2','I2','IV','IV'],  [0,INF,1,-1], 4),
 (44,['I2','IV','III','III'],[INF,0,1,-1], 2),
 (45,['I3','III','III','III'],[INF,eta,eta**2,1], 3),
 (46,['I3','II','III','IV'], [INF,-3,1,0], 3),
 (47,['I4','IV','II','II'],  [INF,0,1,-1], 4),
 (48,['I4','II','III','III'],[INF,-5,ii*sqrt(2),-ii*sqrt(2)], 4),
 (49,['I5','III','II','II'], [INF,0,(1+ii*sqrt(15))**3/8,(1-ii*sqrt(15))**3/8], 5),
 (50,['I6','II','II','II'],  [INF,eta,eta**2,1], 6),
]

# the 18 configurations with exactly one T^- fibre (1-parameter families)
FAMILIES = [
 ['I4','I1','I1','I0*'], ['I2','I2','I2','I0*'], ['I3','I1','II','I0*'],
 ['I2','I1','III','I0*'],['I1','I1','IV','I0*'], ['I1','II','III','I0*'],
 ['I2','II','II','I0*'],
 ['I1','I1','I1','I3*'], ['I1*','I1','I1','I3'],
 ['I1','I1','I2','I2*'], ['I1*','I1','I2','I2'],
 ['I1','I1','I1','III*'],['I1*','I1','I1','III'],
 ['I1','I1','I2','IV*'], ['I1*','I1','I2','II'], ['I1','I1','I2*','II'],
 ['I1','I1','II','IV*'], ['I1*','I1','II','II'],
]

def crossinv(p0, pinf, p1, p2):
    """I = (1+z)^2/z with z = ((p1-p0)(p2-pinf))/((p1-pinf)(p2-p0)); handles oo."""
    def d(a,b):
        if a is INF or b is INF: return None
        return a-b
    num_, den_ = [], []
    # z = (p1-p0)/(p1-pinf) * (p2-pinf)/(p2-p0)
    def ratio(a, u, v):           # (a-u)/(a-v)
        if a is INF: return S(1)
        if u is INF: return S(1)/0 if False else None
        if v is INF: return None
        return (a-u)/(a-v)
    # handle infinities by cases
    pts = [p0,pinf,p1,p2]
    if INF in pts:
        k = pts.index(INF)
        # send oo -> a finite generic value by w = 1/(p - c) with c not a base point
        c = S(7)/13
        while any((not (p is INF)) and simplify(p-c)==0 for p in pts):
            c = c + 1
        pts = [ (S(0) if p is INF else 1/(p-c)) for p in pts ]
    P0,Pi,P1,P2 = pts
    z = cancel(((P1-P0)*(P2-Pi))/((P1-Pi)*(P2-P0)))
    return simplify(cancel((1+z)**2/z))

def assignments(fib, pts):
    """all (i0, iinf, i1, i2) with fib[i0] an I_n (n>=1) cusp and
       psl(fib[i1]) == psl(fib[i2]) and i1<i2."""
    out=[]
    n=len(fib)
    for i0 in range(n):
        f0=fib[i0]
        if not (f0[0]=='I' and f0[1].isdigit() and int(f0[1:].rstrip('*'))>=1 and not f0.endswith('*')):
            continue
        rest=[i for i in range(n) if i!=i0]
        for iinf in rest:
            r2=[i for i in rest if i!=iinf]
            i1,i2=r2
            if psl(fib[i1])!=psl(fib[i2]): continue
            out.append((i0,iinf,i1,i2))
    return out

def table():
    rows=[]
    for num_,fib,pts,dj in RIGID:
        assert sum(euler(f) for f in fib)==12, (num_,fib)
        for (i0,iinf,i1,i2) in assignments(fib,pts):
            Iv = crossinv(pts[i0],pts[iinf],pts[i1],pts[i2])
            rows.append(dict(row=num_,fib=fib,degJ=dj,
                             at0=fib[i0],atinf=fib[iinf],att=(fib[i1],fib[i2]),
                             rho=psl(fib[i1]),dinf=psl(fib[iinf]),I=Iv))
    return rows

if __name__=='__main__':
    T=table()
    print("assignments with rho_1=rho_2 (t_1,t_2 of equal PSL_2 type):",len(T))
    print()
    hdr="%-4s %-24s %-6s %-7s %-7s %-6s %-6s %s"
    print(hdr%("row","configuration","degJ","at 0","at oo","t1,t2","rho","invariant I = A^2/D"))
    for r in sorted(T,key=lambda r:(r['row'],)):
        print(hdr%(r['row']," ".join(r['fib']),r['degJ'],r['at0'],r['atinf'],
                   r['att'][0]+","+r['att'][1],str(r['rho']),
                   str(r['I'])+("  ~ %.6f"%complex(num(r['I'])).real if complex(num(r['I'])).imag==0 else "  ~ %s"%complex(num(r['I'])))))

# ---------------------------------------------------------------- matching ---
def row_types(M,j1,j2):
    """PSL_2 local types of a row of class (M;j1,j2):
       at 0: cusp; at t_1,t_2: from rho = (j1+j2)/(2M); at oo: from delta=(j2-j1)/M."""
    rho = Q(j1+j2,2*M); dinf = Q(abs(j2-j1),M)
    def typ(x):
        f = x - x.__floor__()
        if f == 0: return 'c'
        if f == Q(1,2): return 2
        if f in (Q(1,3),Q(2,3)): return 3
        return None
    return typ(rho), typ(dinf), rho, dinf

def match_row(M,j1,j2,A,C, T=None):
    if T is None: T = table()
    D = C*M*M
    Iv = Q(A*A, D)
    tt, ti, rho, dinf = row_types(M,j1,j2)
    hits=[]
    seen=set()
    for r in T:
        if r['rho']!=tt or r['dinf']!=ti: continue
        try:
            if simplify(r['I'] - Iv) != 0: continue
        except Exception:
            continue
        key=(r['row'],r['at0'],r['atinf'],r['att'])
        if key in seen: continue
        seen.add(key)
        hits.append(r)
    return Iv, (tt,ti), hits
