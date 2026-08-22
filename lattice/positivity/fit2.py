import math,statistics
P="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad/cone_all.csv"
rows=[]
for line in open(P):
    f=line.strip().split(',')
    if len(f)<13: continue
    rows.append(dict(k=float(f[0]),n=int(f[1]),kap=float(f[2]),hc=float(f[3]),
                     l1=float(f[4]),l2=float(f[5]),lc=float(f[6]),rat=float(f[7]),
                     lq=float(f[8]),incone=int(f[9])))
import math
sig=lambda k:12+k*math.log(2)
F=lambda k:0.5*(13.0995887908+14.3931452672-12-k*math.log(2))
for k in (22.4,23.0,23.9):
    R=sorted([r for r in rows if r['k']==k],key=lambda r:r['n'])
    print(f"=== k={k}  sigma={sig(k):.5f}  F={F(k):+.5f} ===")
    # kappa deficit fit:  sigma - kappa_n  vs 1/n and log n / n
    xs=[1.0/r['n'] for r in R]; ys=[sig(k)-r['kap'] for r in R]
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx; b=my-a*mx
    print(f"  sigma-kappa_n = {a:.3f}/n + {b:+.4f}   (deficit at n=80: {sig(k)-R[-1]['kap']:.4f})")
    xs2=[math.log(r['n'])/r['n'] for r in R]
    mx=sum(xs2)/n; sxx=sum((x-mx)**2 for x in xs2); sxy=sum((x-mx)*(y-my) for x,y in zip(xs2,ys))
    a2=sxy/sxx; b2=my-a2*mx
    print(f"  sigma-kappa_n = {a2:.3f}*log n/n + {b2:+.4f}")
    # lambda2/lambda1
    r21=[r['n']*(r['l2']-r['l1']) for r in R]
    print(f"  log(lam2/lam1): max={max(r21):.3f} (n={R[[x for x in r21].index(max(r21))]['n']}), median={statistics.median(r21):.3f}")
    lr=[math.log(r['rat']) for r in R]
    print(f"  ratio<=1.05 at n = {[r['n'] for r in R if r['rat']<=1.05]}")
    print(f"  ratio<=1.5  count = {sum(1 for r in R if r['rat']<=1.5)}/{len(R)};  ratio>10 at n = {[r['n'] for r in R if r['rat']>10]}")
    # liminf statistic:  min of log(ratio) over sliding windows of 10
    ns=[r['n'] for r in R]
    wins=[(ns[i],ns[i+9],min(lr[i:i+10])) for i in range(0,len(R)-9,10)]
    print("  min log(ratio) per block of 10:", " ".join(f"[{a}-{b}]{c:.3f}" for a,b,c in wins))
    print(f"  n^-1 log(cone-min): n=40 {[r['lc'] for r in R if r['n']==40][0]:+.4f}  n=60 {[r['lc'] for r in R if r['n']==60][0]:+.4f}  n=80 {[r['lc'] for r in R if r['n']==80][0]:+.4f}")
    print(f"  n^-1 log(covol)/2 : n=40 {[r['hc'] for r in R if r['n']==40][0]:+.4f}  n=60 {[r['hc'] for r in R if r['n']==60][0]:+.4f}  n=80 {[r['hc'] for r in R if r['n']==80][0]:+.4f}")
    print(f"  delta = -log|form|/log q at n=80: {-R[-1]['lc']/R[-1]['lq']:.5f}")
    print()
