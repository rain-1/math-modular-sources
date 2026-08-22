import math
SP="/tmp/claude-1000/-home-ubuntu-code-math-modular-sources/6e1f664d-43c6-4d6c-a4b9-d6a20e36eaa5/scratchpad"
mn={}
for l in open(SP+"/mn.csv"):
    if "," not in l: continue
    a,b=l.split(','); mn[int(a)]=float(b)
rows=[]
for line in open(SP+"/cone_all.csv"):
    f=line.strip().split(',')
    if len(f)<13 or f[0]!="22.4000": continue
    rows.append((int(f[1]),float(f[2])))
sig=12+22.4*math.log(2)
def ols(xs,ys):
    n=len(xs);mx=sum(xs)/n;my=sum(ys)/n
    sxx=sum((x-mx)**2 for x in xs);sxy=sum((x-mx)*(y-my) for x,y in zip(xs,ys))
    a=sxy/sxx;b=my-a*mx;res=[y-(a*x+b) for x,y in zip(xs,ys)]
    s2=sum(r*r for r in res)/(n-2)
    return a,b,math.sqrt(s2/sxx),math.sqrt(s2)
print(" n  sigma(k)  logM/n   kappa_n   PNTgap   gcdloss")
gl=[];pnt=[];ns=[]
for n,kap in rows:
    p=sig-mn[n]; g=mn[n]-kap
    ns.append(n);pnt.append(p);gl.append(g)
    if n%10==0 or n<8: print(f"{n:3d} {sig:.5f} {mn[n]:.5f} {kap:.5f}  {p:+.5f}  {g:+.5f}")
a,b,se,sd=ols([1/x for x in ns],gl); print(f"\ngcd-loss = {a:.3f}/n {b:+.5f} (se {se:.4f}, sd {sd:.4f})")
a,b,se,sd=ols([math.log(x)/x for x in ns],gl); print(f"gcd-loss = {a:.3f}*log n/n {b:+.5f} (sd {sd:.4f})")
a,b,se,sd=ols([1/x for x in ns],pnt); print(f"PNT gap  = {a:.3f}/n {b:+.5f} (sd {sd:.4f})")
print("n=80: gcd-loss",gl[-1]," PNT gap",pnt[-1])
