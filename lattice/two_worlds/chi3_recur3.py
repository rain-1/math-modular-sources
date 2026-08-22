"""Joint (Q,P) recurrence for the family a = p*k, b = q*k (index k, alpha=p/q)."""
import sys
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
from chi3_recur import nullspace
import sympy as sp

def fitjoint(Qs, Ps, order, deg):
    nunk=(order+1)*(deg+1); rows=[]
    for seq in (Qs,Ps):
        k=0
        while k+order < len(seq):
            rows.append([F(k)**j*seq[k+i] for i in range(order+1) for j in range(deg+1)])
            k+=1
    if len(rows)<nunk+4: return None
    return nullspace(rows)

def analyse(p,q,K=34,maxorder=5,maxdeg=16):
    Qs=[];Ps=[]
    for k in range(K):
        Q,P=row(q*k, p*k)      # b=q*k, a=p*k
        Qs.append(Q); Ps.append(P)
    for order in range(2,maxorder+1):
        for deg in range(1,maxdeg+1):
            if (order+1)*(deg+1)+8 > 2*(K-order): continue
            ns=fitjoint(Qs,Ps,order,deg)
            if ns:
                v=ns[0]; m=sp.Symbol('m'); x=sp.Symbol('x')
                polys=[sp.expand(sum(sp.Rational(v[i*(deg+1)+j].numerator,v[i*(deg+1)+j].denominator)*m**j
                                     for j in range(deg+1))) for i in range(order+1)]
                D=max(sp.degree(sp.Poly(pp,m)) for pp in polys if pp!=0)
                ch=[sp.LC(sp.Poly(pp,m)) if (pp!=0 and sp.degree(sp.Poly(pp,m))==D) else 0 for pp in polys]
                cp=sp.Poly(sum(ch[i]*x**i for i in range(order+1)),x)
                roots=sp.nroots(cp) if cp.degree()>0 else []
                print("alpha=%d/%d : order=%d deg=%d nullity=%d  char=%s"%(p,q,order,deg,len(ns),sp.factor(cp.as_expr())))
                rr=sorted([complex(r) for r in roots], key=lambda z:-abs(z))
                print("           roots(|.| desc): %s"%["%.6f%+.6fi (|.|=%.6f)"%(r.real,r.imag,abs(r)) for r in rr])
                prod=1
                for r in rr: prod*=r
                print("           product = %.6f   log|lam1|=%.6f  log|lam2|=%.6f"%(
                    prod.real, sp.log(abs(rr[0])), sp.log(abs(rr[-1]))))
                return
    print("alpha=%d/%d : none found"%(p,q))

if __name__=="__main__":
    for (p,q) in [(2,1),(7,4),(5,3),(3,2),(4,3),(1,1),(3,4),(1,2)]:
        try: analyse(p,q)
        except Exception as e: print("alpha=%d/%d error %s"%(p,q,e))
