from math import comb
from collections import defaultdict
from itertools import product

def apery(n): return sum(comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1))

# N = (x+y)(z+1)(x+y+z)(y+z+1) as dict on (i,j,k) exponents of x,y,z
def mulnn(a,b,cap):
    r=defaultdict(int)
    for (i1,j1,k1),c1 in a.items():
        for (i2,j2,k2),c2 in b.items():
            i,j,k=i1+i2,j1+j2,k1+k2
            if i<=cap and j<=cap and k<=cap: r[(i,j,k)]+=c1*c2
    return {e:c for e,c in r.items() if c}

fac=[{(1,0,0):1,(0,1,0):1},{(0,0,1):1,(0,0,0):1},{(1,0,0):1,(0,1,0):1,(0,0,1):1},{(0,1,0):1,(0,0,1):1,(0,0,0):1}]
def Npoly(cap):
    p={(0,0,0):1}
    for f in fac: p=mulnn(p,f,cap)
    return p

print("Straub CT check, n<=30:")
allok=True
for n in range(0,31):
    Np=Npoly(n)
    acc={(0,0,0):1}
    for _ in range(n): acc=mulnn(acc,Np,n)
    c=acc.get((n,n,n),0)
    if c!=apery(n): allok=False; print("  MISMATCH n=",n)
print("  all match n<=30:",allok)

# Newton polytope of Lambda = N/(xyz): vertices = exponents of N shifted by (-1,-1,-1)
Nfull={}
p={(0,0,0):1}
for f in fac:
    r=defaultdict(int)
    for e1,c1 in p.items():
        for e2,c2 in f.items():
            e=tuple(e1[i]+e2[i] for i in range(3)); r[e]+=c1*c2
    p={e:c for e,c in r.items() if c}
exps=[tuple(e[i]-1 for i in range(3)) for e in p]
print("Lambda has",len(exps),"monomials; support:",sorted(exps))
# interior lattice points of conv(support): brute force over small box using LP-free test
import itertools
from fractions import Fraction
pts=exps
lo=[min(e[i] for e in pts) for i in range(3)]; hi=[max(e[i] for e in pts) for i in range(3)]
print("bbox",lo,hi)
# convex hull via scipy if available
try:
    from scipy.spatial import ConvexHull
    import numpy as np
    hull=ConvexHull(np.array(pts,dtype=float))
    eqs=hull.equations  # A x + b <= 0
    interior=[]; boundary=[]
    for c in itertools.product(*[range(lo[i],hi[i]+1) for i in range(3)]):
        v=np.array(c,dtype=float)
        d=eqs[:,:3]@v+eqs[:,3]
        if (d< -1e-9).all(): interior.append(c)
        elif (d<1e-9).all(): boundary.append(c)
    print("interior lattice points:",interior)
    print("num boundary lattice points:",len(boundary))
except Exception as ex:
    print("scipy unavailable:",ex)
