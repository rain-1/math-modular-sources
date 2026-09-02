from mpmath import mp, mpf, sqrt, log, quad, mpc, pi, exp, mpmathify
mp.dps=70
# independent check for E_1: P=(1-x)(1-9x)(1-4x); r1=1/9
def P(t): return (1-t)*(1-9*t)*(1-4*t)
# direct tanh-sinh with singular endpoint
for kn,kf in [('B',lambda t:1/(1-t)),('D',lambda t:log(1-t)/(1-t)),('L',lambda t:log(1-t))]:
    v=quad(lambda t: kf(t)/sqrt(P(t)), [0, mpf(1)/9], maxdegree=15)
    print('E_1 A_%s direct  '%kn, mp.nstr(v,45))
    w=quad(lambda t: kf(t)/sqrt(abs(P(t))), [mpf(1)/9, mpf(1)/4], maxdegree=15)
    print('E_1 B_%s direct  '%kn, mp.nstr(w,45))

# branch check: continue H along a semicircle above r1=1/9 and confirm H -> i*|P|^{-1/2}
def Hpath(t, eps):
    z = mpc(t, eps)
    return 1/(sqrt(1-z)*sqrt(1-9*z)*sqrt(1-4*z))
for t in [mpf('0.15'), mpf('0.2')]:
    v = Hpath(t, mpf('1e-25'))
    print('t=',t,' H(t+i0+) =', mp.nstr(v,20), '   i/sqrt|P| =', mp.nstr(mpc(0,1)/sqrt(abs(P(t))),20))
# G: past two branch points
def PG(t): return (1-4*t)*(1-8*t)*(1-12*t)
def HG(t,eps):
    z=mpc(t,eps); return 1/(sqrt(1-4*z)*sqrt(1-8*z)*sqrt(1-12*z))
for t in [mpf('0.1'), mpf('0.2')]:
    print('G t=',t,' H=',mp.nstr(HG(t,mpf('1e-25')),20),'  pred=', mp.nstr((mpc(0,1)**(1 if t<0.125 else 2))/sqrt(abs(PG(t))),20))
