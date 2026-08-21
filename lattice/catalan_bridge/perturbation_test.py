import csv
from math import gcd, log2
def v2(x):
    x=abs(x); return None if x==0 else (x&-x).bit_length()-1
Z={int(r['n']):(int(r['X_n']),int(r['Y_n'])) for r in csv.DictReader(open('zudilin_rows.csv'))}
N={int(r['n']):(int(r['V_n']),int(r['U_n'])) for r in csv.DictReader(open('nesterenko_rows.csv'))}
D={}; cur=1
for i in range(1,600): cur=cur*i//gcd(cur,i); D[i]=cur
print("Is v2(h_n) robust to perturbing U_n by +-1?  (if robust -> v2 is NOT arithmetic)")
print("  n |  v2(h) |  v2(h) with U+1 | v2(h) with U-1 | v2(h) with V+1 | 24n")
for n in (20,40,60,80,98):
    X,Y=Z[n]; V,U=N[n]; S=D[6*n]**2
    a1,a2=X//S,V//S
    base=v2(a1*U-a2*Y)
    p1=v2(a1*(U+1)-a2*Y)
    m1=v2(a1*(U-1)-a2*Y)
    # perturb V -> alpha2 changes
    a2b=(V+S)//S
    pv=v2(a1*U-a2b*Y)
    print("  %3d | %6d | %15d | %14d | %14d | %5d"%(n,base,p1,m1,pv,24*n))
