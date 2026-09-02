"""Exact decomposition of c^{(a)}[log(1-t)/(1-t)] into Li2 terms + elementary part."""
from kummer import *
from bloch import *
from mpmath import conj, arg, nint

def decompose(h, a):
    """returns (points, lambdas, E) with c_D^{(a)} = sum_z lambda_z Li2(z) + E."""
    k, D, s, N = h.k, h.D, h.sigma, h.N
    al = h.alphas()
    coef = {}   # key -> [z, lam]
    def add(z, lam):
        for key,(zz,ll) in coef.items():
            if abs(zz-z) < mpf(10)**(-mp.dps+20):
                coef[key][1] += lam; return
        coef[len(coef)] = [z, lam]
    E = mpc(0)
    for j in range(k):
        for l in range(k):
            c = (-mpf(s)/D)*al[j]**a
            if j == l:
                E += c*(log(1-al[j])**2 - log(-al[j])**2)/2
            else:
                d = al[j]-al[l]
                P1 = (1-al[l])/d; P0 = -al[l]/d
                add(P1,  c); add(P0, -c)
                E += c*( log(1-al[l])*log(1-P1) - log(-al[l])*log(1-P0) )
    E -= (log(N) + mpc(0,1)*pi*h.w())*h.cB(a)
    pts = [coef[i][0] for i in range(len(coef))]
    lams= [coef[i][1] for i in range(len(coef))]
    return pts, lams, E

def parts(h, a):
    pts, lams, E = decompose(h,a)
    total = sum(l*polylog(2,z) for z,l in zip(pts,lams)) + E
    LiRe   = sum(re(l)*re(polylog(2,z)) for z,l in zip(pts,lams))
    Bloch  = -sum(im(l)*BW(z) for z,l in zip(pts,lams))
    ArgP   = sum(im(l)*arg(1-z)*log(abs(z)) for z,l in zip(pts,lams))
    return dict(pts=pts,lams=lams,E=E,total=total,LiRe=LiRe,Bloch=Bloch,ArgP=ArgP,ReE=re(E))
