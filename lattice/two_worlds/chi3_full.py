"""Full arithmetic/archimedean profile of the chi_{-3} well-poised family
   a = p*k, b = q*k (alpha = p/q).  Rates from successive ratios (fast convergence)."""
import sys, math
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
import mpmath as mp

def vp(x,p):
    if x==0: return None
    n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
def den_prime_to(x,p):
    d=x.denominator
    while d%p==0: d//=p
    return d

def profile(p,q,K=26,dps=None):
    if dps is None: dps = 60 + int(8*K*max(1,p/q))
    mp.mp.dps = dps
    L = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])
    Qs=[];Ps=[];Ss=[]
    for k in range(K):
        Q,P=row(q*k,p*k); Qs.append(Q); Ps.append(P)
        Ss.append(mp.mpf(Q.numerator)/Q.denominator*L - mp.mpf(P.numerator)/P.denominator)
    out={}
    out['lam1']=[mp.nstr(Qs[k+1]/Qs[k],8) for k in range(K-4,K-1)]
    out['lam2']=[mp.nstr(Ss[k+1]/Ss[k],8) for k in range(K-4,K-1)]
    v3Q=[vp(Q,3) for Q in Qs]
    V=[Qs[k+1]*Ps[k]-Qs[k]*Ps[k+1] for k in range(K-1)]
    v3V=[vp(x,3) for x in V]
    out['kappa']=[-v3Q[k]/k for k in range(K-3,K)]
    out['w']=[v3V[k]/k for k in range(K-4,K-1)]
    out['nu']=[math.log(den_prime_to(Qs[k],3))/k if den_prime_to(Qs[k],3)>1 else 0.0 for k in range(K-3,K)]
    out['kd']=[math.log(den_prime_to(Ps[k],3))/k for k in range(K-3,K)]
    return out

if __name__=="__main__":
    fams=[(2,1),(7,4),(5,3),(3,2),(4,3),(1,1),(3,4),(1,2)]
    print("%-6s %-24s %-24s %-16s %-16s %-8s %-8s"%("alpha","lambda_1 (last ratios)","lambda_2 (last ratios)","kappa_3","w","nu","k_d"))
    for (p,q) in fams:
        try:
            o=profile(p,q)
            print("%-6s %-24s %-24s %-16s %-16s %-8s %-8s"%("%d/%d"%(p,q),
                ",".join(o['lam1']), ",".join(o['lam2']),
                ",".join("%.3f"%x for x in o['kappa']),
                ",".join("%.3f"%x for x in o['w']),
                ",".join("%.2f"%x for x in o['nu']),
                ",".join("%.2f"%x for x in o['kd'])))
        except Exception as e:
            print("%-6s ERROR %s"%("%d/%d"%(p,q),e))
