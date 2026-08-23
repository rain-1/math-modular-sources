"""Best achievable CDT margin on the EMN z-line host C\\{1} (base 0)."""
import sys, math, cmath, numpy as np
sys.path.insert(0,'/home/ubuntu/code/math-modular-sources/lattice/cdt_finder')
from cdt_bound import tau

def BC(R, N=768):
    zs = np.array([cmath.exp(2j*math.pi*j/N) for j in range(N)])
    fs = 1-np.exp(-R*zs); fp = R*np.exp(-R*zs)
    tot=0.0
    for j in range(N):
        d=fs-fs[j]; dz=zs-zs[j]
        with np.errstate(divide='ignore', invalid='ignore'):
            q=np.abs(d/dz)
        q[j]=abs(fp[j]); tot+=np.sum(np.log(q))
    return tot/N**2

BCcache={}
def bc(R):
    R=round(R,4)
    if R not in BCcache: BCcache[R]=BC(R)
    return BCcache[R]

Rs=[2.2,2.6,3,3.5,4,5,6,7,8,9.5,11,13,16,20,25,30,40,60,90,140,220,350]
print("  the inventory: {1} + Li_1..Li_w (lcm-free, e=k) + H, theta H (type n[1..2n])")
print(f"{'w':>3} {'m':>3} {'tau':>8} {'best entry>0':>13} {'BC':>9} {'bound m<=':>11} {'best margin':>12}")
for w in range(0, 13):
    m = w+3
    cols=[(w+1, 2)]
    e = [0]+list(range(1,w+1))+[1,1]
    T = tau(m, cols, e)['tau']
    bestmar=(-1e9,None); bestbd=(1e9,None)
    for R in Rs:
        ent = math.log(R)-T
        if ent<=0: continue
        B = bc(R)
        mar = m*ent - B
        if mar>bestmar[0]: bestmar=(mar,R)
        if B/ent<bestbd[0]: bestbd=(B/ent,R)
    if bestmar[1] is None:
        print(f"{w:>3} {m:>3} {T:>8.4f}   entry never > 0 on the tested grid"); continue
    Rb=bestmar[1]
    print(f"{w:>3} {m:>3} {T:>8.4f} {math.log(Rb)-T:>13.4f} {bc(Rb):>9.4f} {bestbd[0]:>11.3f} {bestmar[0]:>12.4f}")
print()
print("  margin > 0 would be a contradiction (i.e. an irrationality proof).")
print("  For comparison (CATALAN_MU4 sec.6): mu_4 host A margin -14.04, level 8 entry -0.077.")
