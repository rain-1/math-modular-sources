from mpmath import mp, mpf, nstr, identify, pslq, sqrt, pi, e, exp, log, besseli, besselj, besselk, bessely, gamma, mpc, ber, bei, ker, kei, sinh, cosh, sin, cos
mp.dps = 60
rp = mpf("-0.590491135253101731210234442266594724662621239783171178656390")
rm = mpf("0.423366610269815607602108250217086074203850564724289160158896")
Lp = mpf("1.11812440004378895846495504008")
Lm = mpf("0.891123023067508172851899991326")
print("identify tries (mpmath identify, small constant set):")
for nm,v in [("r+",rp),("r-",rm),("L+",Lp),("L-",Lm),("L+*L-",Lp*Lm),("L+/L-",Lp/Lm),
             ("1/L+",1/Lp),("L+^2",Lp**2),("L-^2",Lm**2),("L+^4",Lp**4),("L+^2*L-^2",(Lp*Lm)**2)]:
    print("  %-8s = %-32s  identify: %s"%(nm, nstr(v,25), identify(v, ['pi','exp(1)','sqrt(2)','log(2)'])))
print()
print("mp.dps=30 identify of L with Bessel-ish constants")
mp.dps = 28
for nm,v in [("L+",Lp),("L-",Lm)]:
    print("  ",nm, identify(+v))
