#!/usr/bin/env python3
"""ap_analyse.py -- post-process ap_out_<N>.txt (5-adic valuations of a_n,b_n,e_n,h_n).
   The valuations are EXACT (from Z/5^MM arithmetic in PARI); the least-squares
   fits below are ordinary floating point."""
import sys, math
p = sys.argv[1]
NN=None; MM=None
rows={}
for line in open(p):
    line=line.strip()
    if line.startswith('#'):
        if line.startswith('# NN='):
            parts=line.split(); NN=int(parts[1][3:]); MM=int(parts[2][3:])
        continue
    f=line.split()
    if len(f)==5:
        rows[int(f[0])]=tuple(int(z) for z in f[1:])
print(f"N={NN}  MM={MM}   (working precision 5^{MM})")

def lsq(ns, ds):
    # d = A + B log n
    n=len(ns); sx=sum(math.log(k) for k in ns); sy=sum(ds)
    sxx=sum(math.log(k)**2 for k in ns); sxy=sum(math.log(k)*d for k,d in zip(ns,ds))
    B=(n*sxy-sx*sy)/(n*sxx-sx*sx); A=(sy-B*sx)/n
    return A,B

# ---- (1b) h_n
top=int(0.8*NN)
grid=[50,100,150,200,250,300,350,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600]
grid=[g for g in grid if g<=top]
print("\n(1b)  h_n = a_n + eta_N b_n")
print(" n     v5(h_n)   v5(h_n)/n    d_n=v5(h_n)-3n")
for g in grid:
    va,vb,ve,vh=rows[g]
    print(f"{g:5d}  {vh:7d}   {vh/g:.6f}   {vh-3*g:+d}")
ns=[k for k in range(2,top+1)]
ds=[rows[k][3]-3*k for k in ns]
A,B=lsq(ns,ds)
print(f"\n  least squares over 2<=n<={top}:  d_n = {A:.4f} + {B:.4f} log n")
print(f"  max d_n = {max(ds)} (at n={ns[ds.index(max(ds))]}),  min d_n = {min(ds)} (at n={ns[ds.index(min(ds))]})")
# restrict to n>=50 too
ns2=[k for k in range(50,top+1)]; ds2=[rows[k][3]-3*k for k in ns2]
A2,B2=lsq(ns2,ds2)
print(f"  least squares over 50<=n<={top}: d_n = {A2:.4f} + {B2:.4f} log n ; max={max(ds2)} min={min(ds2)}")
print(f"  v5(h_n)/n at n={top}: {rows[top][3]/top:.6f}")

# ---- (1c) b_n and e_n
for idx,name in ((1,'b_n = [x^n] E_2^*'),(2,'e_n = [x^n] E\'_{-2}'),(0,'a_n = [x^n] E_2^* E\'_{-2}')):
    vs=[rows[k][idx] for k in range(1,NN+1)]
    print(f"\n(1c)  {name}")
    print("   n :", "  ".join(f"{g}:{rows[g][idx]}" for g in [10,50,100,200,400,600,800,1000,1200,1400,1600,1800,2000] if g in rows))
    print(f"   max v5 over 1<=n<={NN}: {max(vs)}   mean: {sum(vs)/len(vs):.4f}   #(v=0): {vs.count(0)}/{len(vs)}")
    print(f"   v5/n at n={NN}: {rows[NN][idx]/NN:.6f}   at n={top}: {rows[top][idx]/top:.6f}")
    ns3=[k for k in range(50,NN+1)]; ds3=[rows[k][idx] for k in ns3]
    A3,B3=lsq(ns3,ds3)
    print(f"   fit v5 = {A3:.4f} + {B3:.4f} log n  over 50<=n<={NN}")
