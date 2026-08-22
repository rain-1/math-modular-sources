"""WARNING (see consolidation/ONE_CLASS_TWO_WORLDS.md Sec.7): the log-Lambda branch below
is WRONG.  Q_m = 9*sum_j A_j has a sign change at j ~ b/2 and the residues cancel, so
max_y H(y) overestimates log|Q_m| by a constant.  The log-lambda branch (max_x G) is the
correct exponential rate of the linear form, up to the polynomial factor m^{2a-4b}.
Kept as a record; the rates actually used are the measured ones in chi3_design2.py.

Saddle-point rates for the chi_{-3} well-poised family.
R_m(t) = (2t+b+1) prod_{j=1..a}(t-j+1)(t+b+j) / prod_{j=0..b}(t+j+1/3)^2 (t+j+2/3)^2
a = alpha*m, b = beta*m.  Linear form  S_m = Q_m L(2,chi-3) - P_m = sum_{t>=0} R_m(t).
  log lambda = max_{x>alpha} G(x),  log Lambda = max_{0<y<beta} H(y)
"""
import mpmath as mp
mp.mp.dps = 30

def xlx(u):
    u = mp.mpf(u)
    return 0 if u == 0 else u*mp.log(u)

def G(x, al, be):
    return 5*xlx(x) - xlx(x-al) + xlx(x+al+be) - 5*xlx(x+be) - 2*al + 4*be

def H(y, al, be):
    return xlx(y+al) + xlx(al+be-y) - 5*xlx(y) - 5*xlx(be-y) - 2*al + 4*be

def maxG(al, be):
    # G'(x) = log( x^5 (x+al+be) / ((x-al)(x+be)^5) ) = 0
    f = lambda x: 5*mp.log(x) - mp.log(x-al) + mp.log(x+al+be) - 5*mp.log(x+be)
    lo, hi = al*mp.mpf('1.0000001')+mp.mpf('1e-12'), mp.mpf(1000)
    try:
        x0 = mp.findroot(f, [lo, mp.mpf(10)+al], solver='anderson')
    except Exception:
        # scan
        best = None
        xs = [al + mp.mpf(k)/200 for k in range(1, 4000)]
        for x in xs:
            v = G(x, al, be)
            if best is None or v > best[1]: best = (x, v)
        return best
    return x0, G(x0, al, be)

def maxH(al, be):
    best = None
    N = 4000
    for k in range(1, N):
        y = be*mp.mpf(k)/N
        v = H(y, al, be)
        if best is None or v > best[1]: best = (y, v)
    return best

if __name__ == "__main__":
    print(" alpha  beta   logLambda   loglambda   Lambda*lambda   budget(k=2,kap3=?)")
    for be in [1]:
        for al10 in range(2, 21):
            al = mp.mpf(al10)/10
            if al > 2*be: continue
            xg, lg = maxG(al, be)
            yh, lh = maxH(al, be)
            print("  %4.2f  %4.2f  %10s  %10s  %10s" % (
                al, be, mp.nstr(lh,8), mp.nstr(lg,8), mp.nstr(lh+lg,8)))
