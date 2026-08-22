"""Fit the joint (Q,P) recurrence: the true Apery recurrence of the row."""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from chi3_recur import nullspace

def fitjoint(Qs, Ps, order, deg):
    nunk=(order+1)*(deg+1); rows=[]
    for seq in (Qs,Ps):
        m=0
        while m+order < len(seq):
            rows.append([F(m)**j*seq[m+i] for i in range(order+1) for j in range(deg+1)])
            m+=1
    if len(rows)<nunk+4: return None
    return nullspace(rows)

if __name__=="__main__":
    fam=sys.argv[1] if len(sys.argv)>1 else '2m'
    af={'2m':lambda m:2*m,'3m2':lambda m:(3*m)//2,'m':lambda m:m,'7m4':lambda m:(7*m)//4,'5m4':lambda m:(5*m)//4}[fam]
    M=int(sys.argv[2]) if len(sys.argv)>2 else 40
    Qs=[];Ps=[]
    for m in range(M):
        Q,P=row(m,af(m)); Qs.append(Q); Ps.append(P)
    for order in range(2,7):
        for deg in range(1,14):
            if (order+1)*(deg+1)+8 > 2*(M-order): continue
            ns=fitjoint(Qs,Ps,order,deg)
            if ns:
                print("family %s: JOINT recurrence order=%d deg=%d nullity=%d"%(fam,order,deg,len(ns)))
                v=ns[0]
                import sympy as sp
                m=sp.Symbol('m')
                polys=[]
                for i in range(order+1):
                    pl=sum(sp.Rational(v[i*(deg+1)+j].numerator, v[i*(deg+1)+j].denominator)*m**j for j in range(deg+1))
                    polys.append(sp.factor(sp.simplify(pl)))
                    print("   c_%d(m) = %s"%(i,polys[-1]))
                # characteristic roots: leading coefficients
                lead=[sp.limit(p/m**sp.degree(sp.Poly(sp.expand(p),m)), m, sp.oo) if p!=0 else 0 for p in polys]
                degs=[sp.degree(sp.Poly(sp.expand(p),m)) if p!=0 else -1 for p in polys]
                D=max(degs)
                ch=[sp.LC(sp.Poly(sp.expand(p),m)) if (p!=0 and sp.degree(sp.Poly(sp.expand(p),m))==D) else 0 for p in polys]
                print("   degrees:",degs," char poly coeffs (top degree %d):"%D, ch)
                x=sp.Symbol('x')
                cp=sum(ch[i]*x**i for i in range(order+1))
                print("   characteristic polynomial:", sp.factor(cp), " roots:", sp.nsolve if False else sp.solve(cp,x))
                sys.exit(0)
    print("none found")
