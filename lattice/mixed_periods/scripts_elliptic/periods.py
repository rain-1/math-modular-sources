from mpmath import mp, mpf, sqrt, log, cos, sin, pi, quad, mpmathify
import json, itertools

mp.dps = 90

KERNELS = {
 'B': lambda t: 1/(1-t),
 'D': lambda t: log(1-t)/(1-t),
 'L': lambda t: log(1-t),
}

def host_roots(alist):
    rs = sorted([mpf(1)/a for a in alist])
    return rs

def int_0_to_r1(alist, kern):
    """int_0^{r1} k(t)/sqrt(P(t)) dt , P=prod(1-a t), r1 = smallest root."""
    a = sorted(alist)            # largest a  <-> smallest root
    a1 = mpf(a[-1]); r1 = 1/a1
    rest = [mpf(z) for z in a[:-1]]
    S = sqrt(r1)
    def g(s):
        t = r1 - s*s
        Q = mpf(1)
        for z in rest: Q *= (1 - z*t)
        return 2*kern(t)/(sqrt(a1)*sqrt(Q))
    return quad(g, [0, S], maxdegree=12)

def int_between(alist, i, j, kern):
    """int_{r_i}^{r_j} k/sqrt(|P|) over consecutive roots r_i<r_j (indices into sorted roots)."""
    rs = host_roots(alist)
    A = sorted([mpf(z) for z in alist], reverse=True)   # a sorted so that root order matches
    # root rs[k] corresponds to a-value A[k]
    ri, rj = rs[i], rs[j]
    ai, aj = A[i], A[j]
    others = [A[k] for k in range(len(A)) if k not in (i,j)]
    def g(phi):
        u = (1-cos(phi))/2
        t = ri + (rj-ri)*u
        Q = mpf(1)
        for z in others: Q *= abs(1 - z*t)
        return kern(t)/(sqrt(ai*aj)*sqrt(Q))
    return quad(g, [0, pi], maxdegree=12)

HOSTS = {}
for m in [1,2,3,4,5,6,12]:
    HOSTS['E_%d'%m] = ([1,9,4*m], {'d1': mpf(1)/9, 'd2': mpf(1)/(4*m)})
for m in [1,2,3]:
    HOSTS['F_%d'%m] = ([1,25,4*m], {'d1': mpf(1)/25, 'd2': mpf(1)/(4*m)})
HOSTS['G'] = ([4,8,12], {})

out = {}
for name,(alist,_) in HOSTS.items():
    rs = host_roots(alist)
    rec = {'alist': alist, 'roots': [mp.nstr(r,30) for r in rs]}
    for kn,kf in KERNELS.items():
        A = int_0_to_r1(alist, kf)
        B = int_between(alist, 0, 1, kf)
        rec['A_'+kn] = A
        rec['B_'+kn] = B
        if len(rs)==3 and rs[2] < 1:   # G: third fold also in (0,1)
            C = int_between(alist, 1, 2, kf)
            rec['C_'+kn] = C
    out[name]=rec

import pickle
with open('periods.pkl','wb') as f:
    pickle.dump({k:{kk:(str(vv) if not isinstance(vv,list) else vv) for kk,vv in v.items()} for k,v in out.items()}, f)

def S(x): return mp.nstr(x, 62)
for name,(alist,folds) in HOSTS.items():
    rec = out[name]
    print('='*78)
    print(name, ' P(x) = ', ' * '.join('(1-%dx)'%a for a in alist), '  roots:', rec['roots'])
    for kn in 'BDL':
        print('  k=%s  A = int_0^{r1}        = %s' % (kn, S(rec['A_'+kn])))
        print('        B = int_{r1}^{r2}     = %s' % S(rec['B_'+kn]))
        if 'C_'+kn in rec:
            print('        C = int_{r2}^{r3}     = %s' % S(rec['C_'+kn]))
    print()
