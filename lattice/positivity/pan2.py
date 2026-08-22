import math, statistics
SP="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad"
R={}
for l in open(SP+"/pairs.csv"):
    f=l.strip().split(',')
    if len(f)!=10 or f[7]!='1': continue
    R[(int(f[0]),int(f[1]),int(f[2]))]=dict(lc=float(f[3]),l1=float(f[4]),rat=float(f[5]),
        lq=float(f[6]),v2=float(f[8]),lidx=float(f[9]))
def ols(xs,ys):
    n=len(xs);mx=sum(xs)/n;my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs);sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx;b=my-a*mx;res=[y-(a*x+b) for x,y in zip(xs,ys)]
    s2=sum(r*r for r in res)/(n-2)
    return a,b,math.sqrt(s2/sxx),math.sqrt(s2)
print("FIXED RULE  j0=round(0.30*m), j1=j0+1  (no selection bias)")
print("  m | rate/m  | lam1/m  | ratio | logq/m | delta | v2(h)/m bits")
xs=[];ys=[];ds=[];rs=[]
for m in range(10,45):
    j0=round(0.30*m); k=(m,j0,j0+1)
    if k not in R: continue
    r=R[k]; xs.append(m); ys.append(r['lc']); ds.append(-r['lc']/r['lq']); rs.append(math.log(r['rat']))
    if m%2==0 or m>36: print(f" {m:3d} | {r['lc']:+7.4f} | {r['l1']:+7.4f} | {r['rat']:5.2f} | {r['lq']:6.4f} | {-r['lc']/r['lq']:.4f} | {r['v2']/math.log(2):.3f}")
sel=[(x,y) for x,y in zip(xs,ys) if x>=15]
a,b,se,sd=ols([1/x for x,_ in sel],[y for _,y in sel]); print(f"\n  rate = {b:+.4f} {a:+.2f}/m   (1/m fit, sd {sd:.4f})")
a2,b2,se2,sd2=ols([x for x,_ in sel],[y for _,y in sel]); print(f"  rate = ({a2:+.5f} +- {se2:.5f})*m {b2:+.4f}   (linear fit, sd {sd2:.4f})")
a3,b3,se3,sd3=ols([x for x,_ in zip(xs,rs) if x>=15],[y for x,y in zip(xs,rs) if x>=15])
print(f"  log(ratio) = ({a3:+.5f} +- {se3:.5f})*m {b3:+.4f}   (sd {sd3:.4f})")
print("\nZudilin x Nesterenko pair (0, m/3):")
print("  m | rate/m  | ratio | logq/m | delta | v2(h)/m bits  |  best-pair advantage (nats/m)")
for m in range(6,45,3):
    k=(m,0,m//3)
    if k not in R: continue
    r=R[k]
    S=[v for kk,v in R.items() if kk[0]==m]; b=min(S,key=lambda v:v['lc'])
    print(f" {m:3d} | {r['lc']:+7.4f} | {r['rat']:5.2f} | {r['lq']:6.4f} | {-r['lc']/r['lq']:.4f} | {r['v2']/math.log(2):.3f} | {r['lc']-b['lc']:+.4f}")
print("\nv2(h)/m in bits, over ALL pairs, by m:")
for m in [12,20,28,36,40,44]:
    S=[v['v2']/math.log(2) for kk,v in R.items() if kk[0]==m]
    if S: print(f"  m={m:3d}: min {min(S):.2f} median {statistics.median(S):.2f} max {max(S):.2f}   (per unit n=m/3: median {3*statistics.median(S):.2f})")
