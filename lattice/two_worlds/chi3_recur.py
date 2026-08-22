"""Fit a linear recurrence with polynomial coefficients for the chi_{-3} rows."""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row

def nullspace(M):
    """exact nullspace over Q; M list of rows (lists of Fractions). Returns list of basis vectors."""
    M=[r[:] for r in M]; rows=len(M); cols=len(M[0])
    piv=[]; r=0
    for c in range(cols):
        pr=None
        for i in range(r,rows):
            if M[i][c]!=0: pr=i; break
        if pr is None: continue
        M[r],M[pr]=M[pr],M[r]
        pv=M[r][c]
        M[r]=[x/pv for x in M[r]]
        for i in range(rows):
            if i!=r and M[i][c]!=0:
                f=M[i][c]; M[i]=[a-f*b for a,b in zip(M[i],M[r])]
        piv.append(c); r+=1
        if r==rows: break
    free=[c for c in range(cols) if c not in piv]
    basis=[]
    for fc in free:
        v=[F(0)]*cols; v[fc]=F(1)
        for i,c in enumerate(piv): v[c]=-M[i][fc]
        basis.append(v)
    return basis

def fit(seq, order, deg, m0=0):
    """find c_{i,j} with sum_{i=0..order} (sum_j c_{i,j} m^j) seq[m+i] = 0"""
    nunk=(order+1)*(deg+1)
    rows=[]
    m=m0
    while len(rows) < nunk+4 and m+order < len(seq):
        row_=[]
        for i in range(order+1):
            for j in range(deg+1):
                row_.append(F(m)**j*seq[m+i])
        rows.append(row_); m+=1
    if len(rows)<nunk+2: return None
    ns=nullspace(rows)
    return ns

if __name__=="__main__":
    fam = sys.argv[1] if len(sys.argv)>1 else '2m'
    af = {'2m':lambda m:2*m, '3m2':lambda m:(3*m)//2, 'm':lambda m:m}[fam]
    M = 46
    Qs=[]; Ps=[]
    for m in range(M):
        Q,P = row(m, af(m))
        Qs.append(Q); Ps.append(P)
    print("Q_0..Q_6 =", [str(x) for x in Qs[:7]])
    for order in range(2,9):
        for deg in range(1, 16):
            if (order+1)*(deg+1) + 6 > M-order: continue
            ns = fit(Qs, order, deg)
            if ns:
                print("FOUND recurrence: order=%d, poly degree=%d, nullity=%d" % (order,deg,len(ns)))
                v=ns[0]
                for i in range(order+1):
                    print("   c_%d(m) = %s" % (i, [str(v[i*(deg+1)+j]) for j in range(deg+1)]))
                sys.exit(0)
    print("no recurrence found in the searched box")
