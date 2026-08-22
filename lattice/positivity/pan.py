import math, statistics
SP="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad"
R=[]
for l in open(SP+"/pairs.csv"):
    f=l.strip().split(',')
    if len(f)!=10: continue
    if f[7]!='1': continue
    R.append(dict(m=int(f[0]),j1=int(f[1]),j2=int(f[2]),lc=float(f[3]),l1=float(f[4]),
                  rat=float(f[5]),lq=float(f[6]),v2=float(f[8]),lidx=float(f[9])))
print("pairs:",len(R))
print("\n--- best pair at each m (min of m^-1 log(cone-min)) ---")
print("  m | best (j1,j2) | m^-1 log|qG-p| | m^-1 log q | delta | v2(h)log2/m | ZudxNest (0,m/3) | Zud-only best")
xs=[];ys=[];ds=[]
for m in range(4,41):
    S=[r for r in R if r['m']==m]
    if not S: continue
    b=min(S,key=lambda r:r['lc'])
    bd=max(S,key=lambda r:-r['lc']/r['lq'])
    zn=[r for r in S if r['j1']==0 and r['j2']*3==m]
    xs.append(m); ys.append(b['lc']); ds.append(-bd['lc']/bd['lq'])
    print(f" {m:3d} | ({b['j1']:3d},{b['j2']:3d}) | {b['lc']:+9.4f} | {b['lq']:8.4f} | {-b['lc']/b['lq']:.4f} | {b['v2']:7.4f} | "
          + (f"{zn[0]['lc']:+8.4f} (d={-zn[0]['lc']/zn[0]['lq']:.4f})" if zn else "    --    ")
          + f" | best delta {ds[-1]:.4f} at ({bd['j1']},{bd['j2']})")
def ols(xs,ys):
    n=len(xs);mx=sum(xs)/n;my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs);sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx;b=my-a*mx;res=[y-(a*x+b) for x,y in zip(xs,ys)]
    s2=sum(r*r for r in res)/(n-2)
    return a,b,math.sqrt(s2/sxx),math.sqrt(s2)
sel=[(x,y) for x,y in zip(xs,ys) if x>=15]
a,b,se,sd=ols([1/x for x,_ in sel],[y for _,y in sel])
print(f"\nbest m^-1 log(cone-min) = {b:+.4f} {a:+.3f}/m  (m>=15, se {se:.3f}, sd {sd:.4f})")
a2,b2,se2,sd2=ols([x for x,_ in sel],[y for _,y in sel])
print(f"                       = {a2:+.5f}*m {b2:+.4f}  (linear, se {se2:.5f}, sd {sd2:.4f})")
seld=[(x,y) for x,y in zip(xs,ds) if x>=15]
a3,b3,se3,sd3=ols([1/x for x,_ in seld],[y for _,y in seld])
print(f"best delta = {b3:+.4f} {a3:+.3f}/m  (m>=15, sd {sd3:.4f})")
print("\n--- ratio (cone-min / lambda_1) over all 3913 pair lattices ---")
rr=[r['rat'] for r in R]; lrr=[math.log(x) for x in rr]
print(f" min {min(rr):.4f}  median {statistics.median(rr):.4f}  mean log {statistics.mean(lrr):.4f}  90th pct {sorted(rr)[int(.9*len(rr))]:.3f}  max {max(rr):.1f}")
a,b,se,sd=ols([r['m'] for r in R], lrr)
print(f" log(ratio) = ({a:+.5f} +- {se:.5f})*m {b:+.4f}   [slope/se={a/se:+.2f}]  resid sd {sd:.3f}")
for m in [4,8,12,16,20,24,28,32,36,40]:
    S=[r['rat'] for r in R if r['m']==m]
    if S: print(f"   m={m:3d}: N={len(S):4d} min={min(S):.4f} median={statistics.median(S):.3f} mean log={statistics.mean([math.log(x) for x in S]):.3f} frac<=2: {sum(1 for x in S if x<=2)/len(S):.2f}")
