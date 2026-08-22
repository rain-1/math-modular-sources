"""Final design-rule table: chi_{-3} hypergeometric decayers x modular engines at p=3.
Rates averaged over the last few indices to damp the s_3-digit noise."""
import sys, math
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/two_worlds')
from chi3_row import row
import mpmath as mp
LOG3=math.log(3)

def vp(x,p):
    n,d=x.numerator,x.denominator; v=0
    while n%p==0: n//=p; v+=1
    while d%p==0: d//=p; v-=1
    return v
def dpt(x,p):
    d=x.denominator
    while d%p==0: d//=p
    return d

ENGINES = {   # name: (k_eng, sigma_3, log rho_2^eng, log Lambda_eng)
 'C  (10,3,9)'  : (2, 2, 0.0,          math.log(9)),
 'B  (9,3,27)'  : (2, 3, 0.5*math.log(27), 0.5*math.log(27)),
 'F  (17,6,72)' : (2, 2, math.log(8),  math.log(9)),
}

def profile(p,q,K,nav=4):
    mp.mp.dps = 80 + int(12*K*p/q)
    L = mp.nsum(lambda n: (mp.mpf(1)/(3*n+1)**2 - mp.mpf(1)/(3*n+2)**2), [0, mp.inf])
    Qs=[];Ps=[];Ss=[]
    for k in range(K):
        Q,P=row(q*k,p*k); Qs.append(Q); Ps.append(P)
        Ss.append(mp.mpf(Q.numerator)/Q.denominator*L - mp.mpf(P.numerator)/P.denominator)
    V=[Qs[k+1]*Ps[k]-Qs[k]*Ps[k+1] for k in range(K-1)]
    av=lambda f: sum(f(k) for k in range(K-nav,K))/nav
    return dict(
      logLam = av(lambda k: float(mp.log(abs(Qs[k]/Qs[k-1])))),
      loglam = av(lambda k: float(mp.log(abs(Ss[k]/Ss[k-1])))),
      kappa  = av(lambda k: -vp(Qs[k],3)/k),
      w      = av(lambda k: vp(V[k-1],3)/k),
      nu     = av(lambda k: math.log(dpt(Qs[k],3))/k if dpt(Qs[k],3)>1 else 0.0),
      kd     = av(lambda k: math.log(dpt(Ps[k],3))/k),
      K=K-1)

def design(pr, eng):
    keng, sig_e, lrho, lLe = eng
    sig_d = pr['w'] + 2*pr['kappa']
    best=None
    for i in range(1,4000):
        r = 0.005*i
        S   = max(keng*r, pr['kd'])
        eta = pr['kappa']*LOG3 + pr['nu']
        Gd  = min(sig_e*r, sig_d)*LOG3
        Fv  = 0.5*(S + eta + r*lrho + pr['loglam'] - Gd)
        Hv  = Fv + (pr['logLam'] - pr['loglam'])
        if Hv<=0: continue
        d = 1 - Fv/Hv
        if best is None or d>best[0]: best=(d,r,Fv,Hv)
    return best, sig_d

if __name__=="__main__":
    fams=[(2,1,30),(9,5,13),(7,4,16),(5,3,19),(3,2,26),(4,3,19),(5,4,16),(1,1,30)]
    print("%-6s %-3s %8s %9s %7s %7s %7s %7s %7s"%("alpha","k","logLam","loglam","kap3","w","nu","k_d","sig3"))
    profs={}
    for (p,q,K) in fams:
        pr=profile(p,q,K); profs[(p,q)]=pr
        print("%-6s %-3d %8.4f %9.4f %7.3f %7.3f %7.3f %7.3f %7.3f"%(
            "%d/%d"%(p,q),pr['K'],pr['logLam'],pr['loglam'],pr['kappa'],pr['w'],pr['nu'],pr['kd'],
            pr['w']+2*pr['kappa']))
    print()
    print("%-6s | %-28s | %-28s | %-28s"%("alpha","engine C","engine B","engine F"))
    for (p,q,K) in fams:
        cells=[]
        for name in ['C  (10,3,9)','B  (9,3,27)','F  (17,6,72)']:
            (d,r,Fv,Hv),sd = design(profs[(p,q)], ENGINES[name])
            cells.append("delta=%.4f r=%.2f F=%.2f H=%.1f"%(d,r,Fv,Hv))
        print("%-6s | %-28s | %-28s | %-28s"%("%d/%d"%(p,q),*cells))
