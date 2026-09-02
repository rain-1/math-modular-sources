from mpmath import mp, mpf, nstr, factorial, pslq
mp.dps = 400
def minimal(a, N=6000, M=3000):
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
    y=minimal(a); est=[]
    for par in (0,1):
        for K in (30,40):
            ns=[n for n in range(1400,3000) if n%2==par][:K]
            est.append(neville(ns,[factorial(n+1)*abs(y[n]) for n in ns]))
    print("a=%+d :"%a, [nstr(v,52) for v in est])
    Ls[a]=est[1]
mp.dps=45
Lp,Lm = +Ls[1], +Ls[-1]
print("\nL+ =", nstr(Lp,44)); print("L- =", nstr(Lm,44))
r=[mpf("-0.5904911352531017312102344422665947246626212398"), mpf("0.4233666102698156076021082502170860742038505647")]
f=[mpf("0.3273115325873843789497158552322748761624"), mpf("1.4924718374694147454174640727385141942380"),
   mpf("1.3409847218725506933951592687070570547610"), mpf("0.5063756593040270771118689815449755363031")]
print("\ninternal PSLQ at 40+ digits (maxcoeff 10^6):")
for nm,v in [("L+",Lp),("L-",Lm),("L+L-",Lp*Lm),("L+/L-",Lp/Lm),("L+^2",Lp**2)]:
    print("  %-6s : %s"%(nm, pslq([v,mpf(1)]+r+f, maxcoeff=10**6, maxsteps=10**5)))
