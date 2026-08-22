import csv, math, statistics
rows=[]
for line in open(SP:="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad/cone_all.csv"):
    f=line.strip().split(',')
    if len(f)<13: continue
    rows.append(dict(k=float(f[0]),n=int(f[1]),kap=float(f[2]),hc=float(f[3]),
                     l1=float(f[4]),l2=float(f[5]),lc=float(f[6]),rat=float(f[7]),
                     lq=float(f[8]),incone=int(f[9])))
def ols(xs,ys):
    n=len(xs); mx=sum(xs)/n; my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs); sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx; b=my-a*mx
    res=[y-(a*x+b) for x,y in zip(xs,ys)]
    s2=sum(r*r for r in res)/(n-2); se=math.sqrt(s2/sxx)
    return a,b,se,math.sqrt(s2)
for k in sorted(set(r['k'] for r in rows)):
    R=[r for r in rows if r['k']==k]; R.sort(key=lambda r:r['n'])
    n=[r['n'] for r in R]; lr=[math.log(r['rat']) for r in R]
    a,b,se,sd=ols(n,lr)
    # also log-ratio vs log n
    a2,b2,se2,sd2=ols([math.log(x) for x in n],lr)
    # cone-min relative to sqrt(covol):  n*(lc - hc)
    exc=[r['n']*(r['lc']-r['hc']) for r in R]
    a3,b3,se3,sd3=ols(n,exc)
    print(f"k={k}")
    print(f"  N={len(R)}  ratio: min={min(r['rat'] for r in R):.3f} med={statistics.median([r['rat'] for r in R]):.3f} max={max(r['rat'] for r in R):.3f}")
    print(f"  log(ratio) = ({a:+.6f} +- {se:.6f})*n + {b:.4f}   resid sd={sd:.3f}   [slope/se = {a/se:+.2f}]")
    print(f"  log(ratio) = ({a2:+.4f} +- {se2:.4f})*log n + {b2:.4f}  resid sd={sd2:.3f}")
    print(f"  log(cone/sqrt(covol)) = ({a3:+.6f} +- {se3:.6f})*n + {b3:.4f}  resid sd={sd3:.3f}")
    # halves comparison
    h=len(R)//2
    print(f"  mean log(ratio): first half n<={n[h-1]}: {statistics.mean(lr[:h]):.3f}   second half: {statistics.mean(lr[h:]):.3f}")
    print(f"  max log(ratio) over n<=40: {max(x for x,nn in zip(lr,n) if nn<=40):.3f}   over n>40: {max(x for x,nn in zip(lr,n) if nn>40):.3f}")
    print(f"  shortest vector in cone (up to +-): {sum(r['incone'] for r in R)}/{len(R)} = {100*sum(r['incone'] for r in R)/len(R):.0f}%")
    print(f"  n^-1 log(cone-min) at n=80: {[r['lc'] for r in R if r['n']==80][0]:+.4f}   kappa_n={[r['kap'] for r in R if r['n']==80][0]:.4f}")
    print()
