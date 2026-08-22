import sys, math
sys.path.insert(0, 'lattice/adelic_holonomy')
from adelic_bound import adelic
L2 = math.log(2)
# level-8 Catalan host (Zagier E): s = 1/lambda_2 = 1/4
CEIL = math.log(256*0.25)                      # 4.158883  uniformisation ceiling
REAL = math.log(256*0.25*0.6292232680)         # 3.695614  CDT's realised contour loss
BCC  = 11.845 + math.log(0.25)                 # 10.458706
def ev(m, cols, e, slopes, lp, BC, lab):
    r = adelic(m, cols, e, ({2: slopes} if slopes else {}), lp, BC)
    print(f"  {lab:<52s} m={m:3d} tau={r['tau']:.4f} gamma={r['gamma']:+.4f} "
          f"entry={r['entry']:+.4f} margin={r['margin']:+9.3f}")
    return r
print("="*112)
print("LEVEL-8 CATALAN HOST  (Zagier E (12,4,32); lambda_1=8, lambda_2=4, s=1/4, k=2,"
      " b=(2,2), CDT inventory u=(1,3))")
print(f"  ceiling log|phi'(0)| = log(256 s) = {CEIL:.6f}   realised = {REAL:.6f}   BC = {BCC:.6f}")
eCDT = [0,0,1,0,0,0,0,0,0,1,1,1,1,1]
print("\n[A] CDT architecture, 7 pure + 7 conditional")
ev(14,[(1,2),(3,2)],eCDT,None,CEIL,BCC,"archimedean only, at the ceiling")
ev(14,[(1,2),(3,2)],eCDT,None,REAL,BCC,"archimedean only, CDT contour")
ev(14,[(1,2),(3,2)],eCDT,[2]*7+[0]*7,CEIL,BCC,"ADELIC (pure orbit 2-adic slope 2), ceiling")
ev(14,[(1,2),(3,2)],eCDT,[2]*7+[0]*7,REAL,BCC,"ADELIC (pure orbit 2-adic slope 2), CDT contour")
print("\n[B] hypothetical: add d doubly-small functions of 2-adic slope sigma_y")
for sy in (4,8):
    for d in (2,4,7):
        m=14+d; e=[0]*m
        for j in range(6): e[-(j+1)]=1
        ev(m,[(1,2),(3,2)],e,[sy]*d+[2]*7+[0]*7,CEIL,BCC,f"  sigma_y={sy}, d={d}, ceiling")
print("\n[C] best-conceivable pure inventory (u_1=u_2=m/2), pure orbit = m/2")
for m in (14,20,30,50):
    e=[0]*m; h=m//2
    for j in range(6): e[-(j+1)]=1
    ev(m,[(h,2),(h,2)],e,[2]*h+[0]*(m-h),CEIL,BCC,f"  u=(m/2,m/2), adelic, ceiling")
    ev(m,[(h,2),(h,2)],e,None,CEIL,BCC,f"  u=(m/2,m/2), archimedean only, ceiling")
print("\n"+"="*112)
print("X_1(5) Sym^2 HOST (Beukers Thm 4; lambda_1=phi^5, lambda_2=phi^{-5} unit, c=-1, k=3)")
print("  c = -1  =>  sigma_p = v_p(c) = 0 at EVERY prime (Cor 3.2 of CRYSTAL_THEOREM_F);")
print("  lambda_2 is a unit  =>  the pure module Li_j(x/s) has v_p-slope v_p(lambda_2) = 0 too.")
print("  Hence gamma_p = 0 at every prime and the adelic bound coincides with CDT's:")
CE5, BC5 = math.log(256), 11.845
for m,cols,lab in ((14,[(1,2),(3,2)],"CDT inventory, k=3 (b=(2,2,2))"),):
    cols3=[(1,2),(3,2),(0,2)]
    e=[0]*m
    for j in range(6): e[-(j+1)]=1
    ev(m,cols3,e,None,CE5,BC5,lab+", ceiling (normalised)")
    ev(m,cols3,e,None,math.log(256*0.6292232680),BC5,lab+", CDT contour")
