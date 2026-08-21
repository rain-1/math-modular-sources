from fractions import Fraction
import math, mpmath as mp, json
ROWS = [
 ("Apery-Sym1",      [136,68,10,16,-16,4]),
 ("T-Sym1",          [24,12,2,16,-16,4]),
 ("Domb-Sym1",       [20,10,2,64,-64,16]),
 ("AZ(9,3,-27)-Sym1",[72,36,6,-432,432,-108]),
 ("Cooper-s7-Sym1",  [26,13,2,-27,27,-6]),
 ("Cooper-s10-Sym1", [24,12,2,-256,256,-60]),
 ("Cooper-s18-Sym1", [56,28,6,768,-768,180]),
]
N=900
res={}
for nm,rec in ROWS:
    A,Ap,B,C,Cp,Cpp=[Fraction(x) for x in rec]
    a=[Fraction(1),B]; b=[Fraction(0),Fraction(1)]
    for n in range(1,N-1):
        a.append(((A*n*n+Ap*n+B)*a[n]-(C*n*n+Cp*n+Cpp)*a[n-1])/Fraction((n+1)**2))
        b.append(((A*n*n+Ap*n+B)*b[n]-(C*n*n+Cp*n+Cpp)*b[n-1])/Fraction((n+1)**2))
    d=float(A*A-4*C)
    if d<0:
        print(nm,"complex roots -> no limit"); continue
    l1=(float(A)+math.sqrt(d))/2; l2=(float(A)-math.sqrt(d))/2
    if abs(l2)>abs(l1): l1,l2=l2,l1
    dig=int((N-3)*math.log10(abs(l1/l2)))
    mp.mp.dps=min(400,max(50,dig+10))
    n0=N-3
    lim=mp.mpf(b[n0].numerator)/mp.mpf(b[n0].denominator)/(mp.mpf(a[n0].numerator)/mp.mpf(a[n0].denominator))
    s=mp.nstr(lim, min(350, max(40, dig-5)))
    res[nm]=s
    print(nm, "digits~", min(350,dig), s[:50])
json.dump(res, open("sym1_limits.json","w"), indent=1)
with open("sym1_limits.txt","w") as f:
    for nm,s in res.items(): f.write(nm+" "+s+"\n")
