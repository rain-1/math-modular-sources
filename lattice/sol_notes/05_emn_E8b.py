"""Bonus for E8: minimax over polynomials phi(h) of degree d in h=g^2 of  256*sup_g W(g)|phi(g^2)|,
   normalised phi(0)=1.  Compares (1-4h)^2 [E7], 1-ch [E8a], q=1-10h+23h^2 [E8b], and the true optima."""
import numpy as np
from numpy.polynomial import polynomial as Pn
# 256 W(g) = 16 g^2 on (0,1/2];  16 (1-g)^4/g^2 on [1/2,1)
def sup_abs(coefs):
    """coefs = [c0,c1,...] of phi(h); returns sup over (0,1) of 256 W(g)|phi(g^2)|"""
    # phi(g^2) as poly in g
    ph = np.zeros(2*len(coefs)-1); ph[0::2]=coefs
    best=0.0; barg=0.0
    # piece A: 16 g^2 phi(g^2)  -> poly
    A = Pn.polymul([0,0,16.0], ph)
    dA = Pn.polyder(A)
    rts = [r.real for r in np.roots(dA[::-1]) if abs(r.imag)<1e-9] if len(dA)>1 else []
    for x in [0.0,0.5]+[r for r in rts if 0<r<0.5]:
        v=abs(Pn.polyval(x,A))
        if v>best: best,barg=v,x
    # piece B: 16 (1-g)^4 phi(g^2)/g^2 ; use N(g)=16(1-g)^4 phi(g^2), F=N/g^2, F'=(N' g -2N)/g^3
    N = Pn.polymul(16*np.array(Pn.polypow([1.0,-1.0],4)), ph)
    num = Pn.polysub(Pn.polymul(Pn.polyder(N),[0,1.0]), 2*N)
    rts = [r.real for r in np.roots(num[::-1]) if abs(r.imag)<1e-9]
    for x in [0.5, 1.0-1e-12]+[r for r in rts if 0.5<r<1.0]:
        v=abs(Pn.polyval(x,N)/x**2)
        if v>best: best,barg=v,x
    return best,barg

def minimise(d, start):
    cur=np.array(start,dtype=float); curv=sup_abs(np.r_[1.0,cur])[0]; step=1.0
    while step>1e-12:
        improved=False
        for i in range(d):
            for sgn in (1,-1):
                t=cur.copy(); t[i]+=sgn*step
                v=sup_abs(np.r_[1.0,t])[0]
                if v<curv-1e-15: cur,curv=t,v; improved=True
        if not improved: step/=2
    return cur,curv

print("256*sup_g W(g)|phi(g^2)| for various phi(h), h=g^2 :")
for nm,cf in [("1                     ",[1.0]),
              ("(1-4h)^2 = 1-8h+16h^2 ",[1.0,-8,16]),
              ("1 - 4.6240231347 h    ",[1.0,-4.6240231346590769]),
              ("1-10h+23h^2  (q)      ",[1.0,-10,23])]:
    s,a=sup_abs(cf); print(f"   phi = {nm}: sup = {s:.12f}   at g = {a:.10f}")
print()
for d,start in [(1,[-4.0]),(2,[-10.,23.]),(3,[-14.,52.,-60.]),(4,[-18.,100.,-220.,160.])]:
    cf,v=minimise(d,start)
    print(f"   optimal degree-{d} phi(h): 1 + " + " + ".join(f"({c:.8f})h^{i+1}" for i,c in enumerate(cf)) + f"   -> {v:.10f}")
