import math
SP="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad"
R=[]
for l in open(SP+"/grid.csv"):
    f=l.strip().split(',')
    if len(f)!=9: continue
    m,j,s=int(f[0]),int(f[1]),int(f[2])
    R.append(dict(m=m,j=j,s=s,ld=float(f[3]),lv=float(f[4]),obj=float(f[5]),objn=float(f[6]),lq=float(f[7]),dl=float(f[8])))
print("grid points:",len(R),"  m range",min(r['m'] for r in R),max(r['m'] for r in R))
print("\n--- best (smallest) obj = m^-1 log(den*|form|) at each m ---")
print(" m | min_j obj  at j | j=0 (Zudilin) | j=m/3 (Nesterenko, if int) | best delta")
for m in list(range(1,13))+[15,20,24,30,36,40,45,48,50,54,57,60]:
    S=[r for r in R if r['m']==m]
    if not S: continue
    b=min(S,key=lambda r:r['obj'])
    z=[r for r in S if r['j']==0][0]
    ne=[r for r in S if r['j']*3==m]
    bd=max(S,key=lambda r:r['dl'])
    print(f"{m:3d} | {b['obj']:+8.4f} at j={b['j']:3d} | {z['obj']:+8.4f} | "
          + (f"{ne[0]['obj']:+8.4f}" if ne else "    --   ")
          + f" | {bd['dl']:+.5f} (m={m},j={bd['j']})")
print("\n--- global best over the whole grid ---")
b=min(R,key=lambda r:r['obj']); print("min obj per unit m:",b)
b=max(R,key=lambda r:r['dl']); print("max delta        :",b)
print("\n--- trend of min_j obj (per unit m) ---")
xs=[];ys=[]
for m in range(1,61):
    S=[r for r in R if r['m']==m]
    b=min(S,key=lambda r:r['obj']); xs.append(m); ys.append(b['obj'])
for m in range(1,61,4): print(f"  m={m:3d} min_j obj={ys[m-1]:+.4f}  argmin j={min([r for r in R if r['m']==m],key=lambda r:r['obj'])['j']}")
# fit min obj -> asymptote  a + b/m
def ols(xs,ys):
    n=len(xs);mx=sum(xs)/n;my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs);sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx;b=my-a*mx; res=[y-(a*x+b) for x,y in zip(xs,ys)]
    return a,b,math.sqrt(sum(r*r for r in res)/(n-2))
sel=[(x,y) for x,y in zip(xs,ys) if x>=20]
a,b,sd=ols([1/x for x,_ in sel],[y for _,y in sel])
print(f"\nmin_j obj = {b:+.4f} + {a:.3f}/m   (m>=20, resid sd {sd:.4f})  -> per unit n (=m/3): {3*b:+.4f}")
