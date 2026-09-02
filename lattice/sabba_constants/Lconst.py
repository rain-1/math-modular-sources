"""L_a = lim (n+1)! |y_n^min| -- there is a parity-dependent 1/n correction, so
   extrapolate along each parity class separately."""
from mpmath import mp, mpf, nstr, factorial, pslq
mp.dps = 200
def minimal(a, N=2200, M=900):
    y=[mpf(0)]*(N+2); y[N]=mpf(1)
    for n in range(N+1,1,-1): y[n-2]=(y[n]-n*y[n-1])/(a*(-1)**n)
    return [v/y[0] for v in y[:M]]
def neville(ns, us):
    K=len(ns); h=[mpf(1)/n for n in ns]; T=[[mpf(0)]*K for _ in range(K)]
    for i in range(K): T[i][0]=us[i]
    for j in range(1,K):
        for i in range(K-j): T[i][j]=(T[i+1][j-1]*(-h[i])-T[i][j-1]*(-h[i+j]))/(h[i+j]-h[i])
    return T[0][K-1]
Ls={}
for a in (1,-1):
    y=minimal(a); out=[]
    for par in (0,1):
        for K in (12,16,20):
            ns=[n for n in range(300,900) if n%2==par][:K]
            out.append(neville(ns,[factorial(n+1)*abs(y[n]) for n in ns]))
    print("a=%+d  L estimates:"%a, [nstr(v,26) for v in out])
    Ls[a]=out[1]
print()
print("L+ =", nstr(Ls[1],26)); print("L- =", nstr(Ls[-1],26))
mp.dps=24
Lp, Lm = +Ls[1], +Ls[-1]
r=[mpf("-0.59049113525310173121023444226659472466"), mpf("0.42336661026981560760210825021708607420")]
f=[mpf("0.32731153258738437894971585523227487616"), mpf("1.49247183746941474541746407273851419424"),
   mpf("1.34098472187255069339515926870705705476"), mpf("0.50637565930402707711186898154497553630")]
print("\ninternal PSLQ (dps=24, maxcoeff 10^4):")
for nm,v in [("L+",Lp),("L-",Lm),("L+L-",Lp*Lm),("L+/L-",Lp/Lm)]:
    print("  %-6s vs [1,r+,r-,f+(1),f+(-1),f-(1),f-(-1)] : %s"%(nm,
        pslq([v,mpf(1)]+r+f, maxcoeff=10**4, maxsteps=10**5)))
