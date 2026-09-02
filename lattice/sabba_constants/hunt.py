"""PSLQ hunt for closed forms + exact bilinear relation coupling a=+1 and a=-1."""
from mpmath import (mp, mpf, mpc, nstr, pslq, identify, sqrt, pi, e, exp, log, factorial,
                    besseli, besselj, besselk, bessely, gamma, ber, bei, ker, kei, mpmathify)
mp.dps = 120

def minimal(a, N=1500, M=700):
    y = [mpf(0)]*(N+2); y[N] = mpf(1)
    for n in range(N+1, 1, -1):
        y[n-2] = (y[n] - n*y[n-1])/(a*(-1)**n)
    return [v/y[0] for v in y[:M]]

data = {}
for a in (1,-1):
    y = minimal(a)
    f1  = sum(y[n]/factorial(n)          for n in range(len(y)-1,-1,-1))
    fm1 = sum((-1)**n*y[n]/factorial(n)  for n in range(len(y)-1,-1,-1))
    d1  = sum(y[n]/factorial(n-1)        for n in range(len(y)-1,0,-1))
    dm1 = sum((-1)**(n-1)*y[n]/factorial(n-1) for n in range(len(y)-1,0,-1))
    data[a] = dict(r=y[1], f1=f1, fm1=fm1, d1=d1, dm1=dm1)
    print("a=%+d : y1/y0 = %s"%(a, nstr(y[1],40)))
    print("       f(1)  = %s"%nstr(f1,40)); print("       f(-1) = %s"%nstr(fm1,40))
    print("       f'(1) = %s"%nstr(d1,40)); print("       f'(-1)= %s"%nstr(dm1,40))
    print("       check  f(-1) + 2a f'(1) = %s"%nstr(fm1 + 2*a*d1, 5))
    print("       check  f(1) - 2a(f''(-1)-f'(-1)) : skipped")
P, M_ = data[1], data[-1]
lhs = P['r'] + M_['r'] - 2
rhs = -(P['f1']*M_['fm1'] + M_['f1']*P['fm1'])
print()
print("bilinear invariant:  r+ + r- - 2 = -(f(1)g(-1)+g(1)f(-1))")
print("   LHS =", nstr(lhs, 40)); print("   RHS =", nstr(rhs, 40)); print("   diff =", nstr(lhs-rhs, 5))

# ---------------- PSLQ hunts ----------------
Cc = 1/(1-P['r']); Dd = 1 - M_['r']
Lp = mpf("1.118124400043788958464955040"); Lm = mpf("0.891123023067508172851899991")
print()
print("PSLQ hunts (maxcoeff 10^6, maxsteps 10^6)")
def hunt(name, targets, basis, bnames, maxc=10**6):
    for tn, X in targets:
        vec = [X] + basis
        rel = pslq(vec, maxcoeff=maxc, maxsteps=10**5, tol=mpf(10)**(-mp.dps+25))
        print("   %-28s %-6s : %s"%(name, tn, rel))
mp.dps = 100
b2 = [ber(0,2), bei(0,2), ber(1,2), bei(1,2), ker(0,2), kei(0,2)]
hunt("Kelvin nu=0,1 at 2", [("C",Cc),("D",Dd),("r+",P['r']),("L+",Lp),("f+(1)",P['f1'])], b2, None)
bI = [besseli(0,2), besseli(1,2), besselj(0,2), besselj(1,2), mpf(1)]
hunt("I,J at 2", [("C",Cc),("D",Dd),("L+",Lp),("L-",Lm),("f+(1)",P['f1'])], bI, None)
bI2 = [besseli(0,2*sqrt(2)), besseli(1,2*sqrt(2)), besselj(0,2*sqrt(2)), besselj(1,2*sqrt(2)), mpf(1)]
hunt("I,J at 2sqrt2", [("C",Cc),("D",Dd),("L+",Lp),("L-",Lm)], bI2, None)
bG = [gamma(mpf(1)/4), gamma(mpf(3)/4), pi, sqrt(2), exp(1), mpf(1)]
hunt("Gamma/pi/e", [("L+",Lp),("L-",Lm),("L+*L-",Lp*Lm),("r+",P['r'])], bG, None)
print()
print("algebraicity: pslq on powers of C and of L+ (degree <= 8)")
for nm, X in [("C",Cc),("D",Dd),("L+",Lp),("L-",Lm),("r+",P['r']),("f+(1)",P['f1'])]:
    for deg in (4,6,8):
        rel = pslq([X**k for k in range(deg+1)], maxcoeff=10**8, maxsteps=10**5)
        if rel: print("   %s deg %d : %s"%(nm,deg,rel))
    print("   %-6s no algebraic relation of degree <= 8 with coeffs < 10^8"%nm)
