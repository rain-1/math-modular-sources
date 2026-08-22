"""3-adic alignment test for the chi_{-3} well-poised hypergeometric rows
against Zagier row C (engine).  Prediction (Theorem F + Conj D):
   xi_3(ours) = zeta_3(2) = 2 * xi_3(C),   since xi_inf(ours) = L, xi_inf(C) = L/2,
   and the 3-Euler factor E_3(2) = 1 - chi_{-3}(3)/9 = 1.
"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_family import row

def rowC(N):
    a=[F(1),F(3)]; b=[F(0),F(1)]
    for n in range(1,N):
        a.append(((10*n*n+10*n+3)*a[n]-9*n*n*a[n-1])/F((n+1)**2))
        b.append(((10*n*n+10*n+3)*b[n]-9*n*n*b[n-1])/F((n+1)**2))
    return a,b

def v3(x):
    if x==0: return 10**9
    n,d=x.numerator,x.denominator; v=0
    while n%3==0: n//=3; v+=1
    while d%3==0: d//=3; v-=1
    return v

if __name__=="__main__":
    N=int(sys.argv[1]) if len(sys.argv)>1 else 160
    a,b=rowC(N)
    # 3-adic Cauchy check for row C
    print("row C 3-adic increments v3(b_n/a_n - b_{n-1}/a_{n-1}):",
          [v3(b[n]/a[n]-b[n-1]/a[n-1]) for n in range(N-4,N)])
    xiC = b[N-1]/a[N-1]           # ~ zeta_3(2)/2 to 3-adic precision ~2N
    prec = 2*(N-1)-10
    print("row C: xi_3 known to about 3^%d" % prec)
    print()
    for cname, af in [("a=2m", lambda m:2*m), ("a=3m/2", lambda m:(3*m)//2),
                      ("a=m", lambda m:m), ("a=7m/4", lambda m:(7*m)//4)]:
        print("--- family %s ---" % cname)
        print("   m   v3(Q)  v3(P)  v3(P/Q - 2*xiC)   [prediction: -> infinity]")
        prev=None
        for m in range(2,25,2):
            Q,P=row(m,af(m))
            d = P/Q - 2*xiC
            print("  %3d  %5d  %5d   %s" % (m, v3(Q), v3(P), v3(d)))
