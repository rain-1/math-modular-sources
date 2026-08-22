"""The arithmetic bridge: cross-determinant divisibility between the chi_{-3}
hypergeometric row (decayer) and Zagier row C (engine), at balanced sampling.

Delta(m,n) = 2*Q_m*b_n - a_n*P_m   (rows: Q_m L - P_m ,  a_n (L/2) - b_n)
v_3(Delta) = v_3(Q_m)+v_3(a_n) + min( sigma_C * n , sigma_dec * m ) + O(log).
"""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from chi3_padic import rowC, v3

if __name__=="__main__":
    p,q = int(sys.argv[1]), int(sys.argv[2])     # alpha = p/q
    rat = sys.argv[3] if len(sys.argv)>3 else "1:1"   # n : k  sampling
    A,B = map(int, rat.split(":"))
    N = int(sys.argv[4]) if len(sys.argv)>4 else 12
    a,b = rowC(A*N+2)
    print("alpha=%d/%d  sampling  n=%d*s (row C), k=%d*s (hypergeometric)"%(p,q,A,B))
    print("   s   v3(Q)  v3(a_n)  v3(Delta)  v3(Delta)-v3(Q)-v3(a_n)")
    for s in range(1,N+1):
        k = B*s; n = A*s
        Q,P = row(q*k, p*k)
        D = 2*Q*b[n] - a[n]*P
        print("  %3d  %6d %7d  %9d  %s"%(s, v3(Q), v3(a[n]), v3(D), v3(D)-v3(Q)-v3(a[n])))
