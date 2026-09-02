"""Final summary of all numerical data (regenerates every number quoted in REPORT.md)."""
from mpmath import mp, mpf, nstr, pslq, factorial, identify
mp.dps = 130
def cf(b, aa, N=600):
    Am2, Am1 = mpf(1), b(0); Bm2, Bm1 = mpf(0), mpf(1)
    for n in range(1, N+1):
        An=b(n)*Am1+aa(n)*Am2; Bn=b(n)*Bm1+aa(n)*Bm2
        Am2,Am1,Bm2,Bm1 = Am1,An,Bm1,Bn
    return Am1/Bm1
C = cf(lambda n: mpf(n), lambda n: mpf(1) if n==1 else mpf((-1)**n))
D = cf(lambda n: mpf(n+1), lambda n: mpf((-1)**n))
print("C = %s"%nstr(C,100)); print("D = %s"%nstr(D,100))
print("OEIS A244279/A244280 comment value 0.628736607098954801603428 -- differs from C at the 16th decimal.")
def minimal(a, N=1500, M=700):
    y=[mpf(0)]*(N+2); y[N]=mpf(1)
    for n in range(N+1,1,-1): y[n-2]=(y[n]-n*y[n-1])/(a*(-1)**n)
    return [v/y[0] for v in y[:M]]
def Lc(a, K=12, par=0):
    # NB: (n+1)!|y_n| has a PARITY-dependent 1/n correction, so extrapolate within
    # one parity class only (see Lconst2.py for the 44-digit computation).
    y=minimal(a); ns=[n for n in range(300,700) if n%2==par][:K]
    us=[factorial(n+1)*abs(y[n]) for n in ns]; h=[mpf(1)/n for n in ns]
    T=[[mpf(0)]*K for _ in range(K)]
    for i in range(K): T[i][0]=us[i]
    for j in range(1,K):
        for i in range(K-j): T[i][j]=(T[i+1][j-1]*(-h[i])-T[i][j-1]*(-h[i+j]))/(h[i+j]-h[i])
    return T[0][K-1]
res={}
for a in (1,-1):
    y=minimal(a)
    f1=sum(y[n]/factorial(n) for n in range(len(y)-1,-1,-1))
    fm1=sum((-1)**n*y[n]/factorial(n) for n in range(len(y)-1,-1,-1))
    res[a]=(y[1], f1, fm1, Lc(a))
    print("\na = %+d"%a)
    print("  y_n/y_0, n=0..7:", [nstr(y[n],22) for n in range(8)])
    print("  r_a = y1/y0 = %s"%nstr(y[1],60))
    print("  L_a (lim (n+1)!|y_n|) = %s"%nstr(res[a][3],28))
    print("  f_a(1)  = %s"%nstr(f1,40)); print("  f_a(-1) = %s"%nstr(fm1,40))
print("\nC - 1/(1-r_+) = %s"%nstr(C - 1/(1-res[1][0]),5))
print("D - (1 - r_-) = %s"%nstr(D - (1-res[-1][0]),5))
print("\nPSLQ among the problem's own constants (looking for internal relations):")
mp.dps=100
basis=[mpf(1), res[1][0], res[-1][0], res[1][1], res[1][2], res[-1][1], res[-1][2]]
names=["1","r+","r-","f+(1)","f+(-1)","f-(1)","f-(-1)"]
for nm,v in [("L+",res[1][3]),("L-",res[-1][3]),("L+*L-",res[1][3]*res[-1][3]),("1/L+",1/res[1][3])]:
    print("   %-6s : %s"%(nm, pslq([v]+basis, maxcoeff=10**5, maxsteps=10**5)))
print("   f+(1)*f-(1)  vs 1,L+,L-:", pslq([res[1][1]*res[-1][1], mpf(1), res[1][3], res[-1][3]], maxcoeff=10**6, maxsteps=10**5))
