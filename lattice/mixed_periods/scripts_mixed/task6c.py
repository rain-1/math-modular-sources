from mpmath import mp, mpf, sqrt, log, pi, quad, clsin, nstr, fabs, pslq, atan, zeta
mp.dps=120
def Lval2():
    return 3**(-mpf(2))*(zeta(2,mpf(1)/3)-zeta(2,mpf(2)/3))
L1=pi/(3*sqrt(3)); L2=Lval2()
m=1;D=3;a=sqrt(D);th=2*atan(1/a)
cB=th/a; cD=(th*log(mpf(D)/m)-2*clsin(2,pi-th))/a
print("== cross-check of our c_B,c_D against CDT's monodromy coefficients (m=1) ==")
print("c_B                            =",nstr(cB,45))
print("L(1,chi_-3) = pi/(3 sqrt3)     =",nstr(L1,45),"  diff",nstr(fabs(cB-L1),3))
print("c_D                            =",nstr(cD,45))
print("L(1,chi_-3)log3 - L(2,chi_-3)  =",nstr(L1*log(3)-L2,45),"  diff",nstr(fabs(cD-(L1*log(3)-L2)),3))
print("L(2,chi_-3)                    =",nstr(L2,45))
print()
def period(m,i,j,l):
    D=4*m-1
    return quad(lambda u: ((1-u*u)/(4*m))**i*((D+u*u)/mpf(4*m))**(-j)*log((D+u*u)/mpf(4*m))**l,[0,1])/(2*m)
E=period(1,0,1,2)
print("== bonus: c[log^2(1-t)/(1-t)] at m=1 =",nstr(E,50))
Cl2=clsin(2,pi/3); z3=zeta(3); l3=log(3)
cand=[("E",E),("1",mpf(1)),("pi/sqrt3",pi/sqrt(3)),("Cl2(pi/3)/sqrt3",Cl2/sqrt(3)),
      ("pi log3/sqrt3",pi*l3/sqrt(3)),("z3",z3),("z3/sqrt3",z3/sqrt(3)),
      ("pi^2 log3",pi**2*l3),("pi^2 log3/sqrt3",pi**2*l3/sqrt(3)),
      ("log^3 3",l3**3),("pi Cl2(pi/3)",pi*Cl2),("pi Cl2(pi/3)/sqrt3",pi*Cl2/sqrt(3)),
      ("Cl2(pi/3) log3",Cl2*l3),("Cl2(pi/3)log3/sqrt3",Cl2*l3/sqrt(3)),
      ("pi^3",pi**3),("pi^3/sqrt3",pi**3/sqrt(3)),("pi log^2 3",pi*l3**2),("pi log^2 3/sqrt3",pi*l3**2/sqrt(3))]
for k in range(4,len(cand)+1):
    r=pslq([c[1] for c in cand[:k]],maxcoeff=10**5,maxsteps=10**6,tol=mpf(10)**-80)
    if r and r[0]!=0:
        print("  relation with first %d basis elts:"%k, {cand[i][0]:r[i] for i in range(k) if r[i]!=0}); break
else:
    print("  no relation found over the 18-element weight<=3 level-6 basis (maxcoeff 1e5, tol 1e-80)")
