from mpmath import mp, mpf, quad, atanh, catalan, log, sqrt, binomial
mp.dps = 40
print("== A3'. Phi(1)=sum a_n : partial sums with 1/N tail extrapolation")
def af(n):
    return mpf(1)/(mpf(2)**(n+1)*(n+1))*sum(binomial(n,k)/mpf(2*k+1) for k in range(n+1))
ps=mpf(0); marks={}
for n in range(400):
    ps+=af(n)
    if n+1 in (100,200,400): marks[n+1]=ps
p1,p2,p4=marks[100],marks[200],marks[400]
rich = p4 + (p4-p2)   # 1/N tail: S = S_N + c/N ;  S ~ 2*S_2N - S_N
print("   S_100=",p1,"\n   S_200=",p2,"\n   S_400=",p4)
print("   Richardson 2*S_400-S_200 =", 2*p4-p2, "\n   G                        =", catalan,
      "\n   diff =", 2*p4-p2-catalan)
print("\n== A7'. branch of H at z=1, H(z)=2*int_0^q artanh u/(1+u^2)du, q=sqrt(z/(2-z))")
def Hclosed(zz):
    q = sqrt(zz/(2-zz))
    return 2*quad(lambda t: atanh(t)/(1+t*t), [0,q])
print("   H(1)=",Hclosed(mpf(1)),"  G=",catalan)
for e in ['1e-3','1e-5','1e-7','1e-9','1e-11']:
    eps=mpf(e); zz=1-eps; dif=Hclosed(zz)-catalan
    print("   1-z=%-6s H-G=%-28s (H-G)/((1-z)log(1-z)) = %s"%(e,mp.nstr(dif,12),mp.nstr(dif/(eps*log(eps)),12)))
