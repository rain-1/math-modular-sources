from kummer import *
from bloch import *
import sys
setup(110)

HOSTS = [ (3,9,'M'), (3,63,'P'), (3,126,'M'), (4,80,'P'), (4,624,'P'),
          (5,1025,'P'), (5,7775,'M'), (6,15624,'P') ]
for (k,N,fam) in HOSTS:
    h = Host(k,N,fam); D=h.D; n = mp.nint(mpf(D)**(mpf(1)/k))
    print(f"\n===== k={k} N={N} fam={fam} D={D} = {int(n)}^{k}  (sigma={h.sigma}, w={mp.nstr(h.w(),3)}) =====")
    al = h.alphas()
    print("  alphas / n :", [mp.nstr(a/n,8) for a in al])
    P0=[];P1=[]
    for j in range(k):
        for l in range(k):
            if j==l: continue
            d=al[j]-al[l]
            P0.append((l,j,-al[l]/d)); P1.append((l,j,(1-al[l])/d))
    # distinct points
    def distinct(L):
        out=[]
        for (l,j,z) in L:
            if not any(abs(z-q)<mpf(10)**-80 for q in out): out.append(z)
        return out
    d0=distinct(P0); d1=distinct(P1)
    print(f"  #distinct P0 points = {len(d0)}   #distinct P1 points = {len(d1)}")
    print("  P0 points (= zeta_l/(zeta_l-zeta_j), independent of n):")
    for z in d0: print(f"     {mp.nstr(z,12)}     |z|={mp.nstr(abs(z),8)}   D(z)={mp.nstr(BW(z),20)}")
    print("  P1 points (= (n*zeta_l-1)/(n*(zeta_l-zeta_j))):")
    for z in d1: print(f"     {mp.nstr(z,12)}     |z|={mp.nstr(abs(z),8)}   D(z)={mp.nstr(BW(z),20)}")
