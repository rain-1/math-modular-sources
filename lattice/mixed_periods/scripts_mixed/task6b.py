from mpmath import mp, mpf, log, sqrt, pi, quad, nstr, fabs, clsin, atan, pslq
mp.dps=80
# c_C = int_0^{1/4} log(1-t)/(t sqrt(1-4t)) dt ; u-sub t=(1-u^2)/4 -> 2 int_0^1 log((3+u^2)/4)/(1-u^2) du
cC1 = quad(lambda t: log(1-t)/(t*sqrt(1-4*t)), [0, mpf(1)/4])
cC2 = 2*quad(lambda u: log((3+u*u)/mpf(4))/(1-u*u), [0,1])
print("c_C direct (t-integral)   =", nstr(cC1,50))
print("c_C via u-substitution    =", nstr(cC2,50))
print("-pi^2/18                  =", nstr(-pi**2/18,50))
print("|c_C - (-pi^2/18)|        =", nstr(fabs(cC2+pi**2/18),3))
print("|direct - usub|           =", nstr(fabs(cC1-cC2),3))
# also c_D at m=1 vs the closed form and the Q_1 shape
m=1;D=3;a=sqrt(D);th=2*atan(1/a)
print()
print("m=1: theta =", nstr(th,30), "  pi/3 =", nstr(pi/3,30))
