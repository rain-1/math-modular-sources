"""Profile the chi_{-3} well-poised family and evaluate the two-row master formula
   against Zagier row C (the engine).  All rates are per unit of the row's own index k.

Engine row C: (a,b,c)=(10,3,9); a_n in Z, d_n^2 b_n in Z (k_eng=2);
  Lambda_eng=9, |lambda_2^eng| = 1  =>  log rho_2^eng = 0;  sigma_3^C = v_3(9) = 2.
Decayer: our row at index k.
  F = 1/2 [ S + eta + alpha_s*log rho_2^eng + gamma*log lambda - Gdiv ]
  S = max(2*alpha_s, kd*gamma), eta = gamma*(kappa*log3 + nu),
  Gdiv = min(2*alpha_s, sigma*gamma)*log3,  sigma = w + 2*kappa
  H = F + gamma*(logLambda - log lambda),  delta = 1 - F/H
"""
import sys, math
from fractions import Fraction as F
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
import mpmath as mp
LOG3 = math.log(3)

def vp(x,p):
    if x==0: return None
    n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
def dpt(x,p):
    d=x.denominator
    while d%p==0: d//=p
    return d

def profile(p,q,K):
    mp.mp.dps = 60 + int(12*K*p/q)
    L = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])
    Qs=[];Ps=[];Ss=[]
    for k in range(K):
        Q,P=row(q*k,p*k); Qs.append(Q); Ps.append(P)
        Ss.append(mp.mpf(Q.numerator)/Q.denominator*L - mp.mpf(P.numerator)/P.denominator)
    K1=K-1
    logL1 = float(mp.log(abs(Qs[K1]/Qs[K1-1])))
    logL2 = float(mp.log(abs(Ss[K1]/Ss[K1-1])))
    kappa  = -vp(Qs[K1],3)/K1
    V=[Qs[k+1]*Ps[k]-Qs[k]*Ps[k+1] for k in range(K-1)]
    w = vp(V[K1-1],3)/K1
    nu = math.log(dpt(Qs[K1],3))/K1 if dpt(Qs[K1],3)>1 else 0.0
    kd = math.log(dpt(Ps[K1],3))/K1
    return dict(logLam=logL1, loglam=logL2, kappa=kappa, w=w, nu=nu, kd=kd, sigma=w+2*kappa, K=K1)

def design(pr, rgrid=None):
    """optimise over the sampling ratio r = alpha_s/gamma (gamma = 1)."""
    best=None
    if rgrid is None: rgrid=[0.05*i for i in range(1,120)]
    for r in rgrid:
        S = max(2*r, pr['kd'])
        eta = pr['kappa']*LOG3 + pr['nu']
        Gd  = min(2*r, pr['sigma'])*LOG3
        Fv = 0.5*(S + eta + 0.0 + pr['loglam'] - Gd)
        Hv = Fv + (pr['logLam'] - pr['loglam'])
        if Hv<=0: continue
        d = 1 - Fv/Hv
        if best is None or d>best[0]: best=(d,r,Fv,Hv)
    return best

if __name__=="__main__":
    fams=[(2,1,30),(7,4,15),(5,3,18),(3,2,26),(4,3,18),(5,4,15),(1,1,30),(3,4,20),(1,2,26)]
    print("%-6s %-4s %9s %9s %8s %8s %7s %7s %8s | %7s %7s %8s %8s"%(
        "alpha","k","logLam","loglam","kappa3","w","nu","k_d","sigma3","delta","r_opt","F","H"))
    for (p,q,K) in fams:
        try:
            pr=profile(p,q,K)
            d,r,Fv,Hv=design(pr)
            print("%-6s %-4d %9.4f %9.4f %8.3f %8.3f %7.3f %7.3f %8.3f | %7.4f %7.2f %8.3f %8.3f"%(
                "%d/%d"%(p,q),pr['K'],pr['logLam'],pr['loglam'],pr['kappa'],pr['w'],pr['nu'],pr['kd'],pr['sigma'],
                d,r,Fv,Hv))
        except Exception as e:
            print("%-6s ERROR %s"%("%d/%d"%(p,q),e))
