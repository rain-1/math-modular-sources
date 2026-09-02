from ident2 import *
import sys, itertools
setup(150)
CL=im(polylog(2,exp(mpc(0,1)*pi/3)))
def CM(a,b,c):
    """Cassaigne-Maillot: pi*m(a+b x+c y)"""
    a,b,c=mpf(a),mpf(b),mpf(c)
    if a>b+c or b>a+c or c>a+b: return pi*log(max(a,b,c))
    cl=lambda x: mpf(-1) if x<-1 else (mpf(1) if x>1 else x)
    al=mp.acos(cl((b*b+c*c-a*a)/(2*b*c))); be=mp.acos(cl((a*a+c*c-b*b)/(2*a*c))); ga=mp.acos(cl((a*a+b*b-c*c)/(2*a*b)))
    return re(BW((b/a)*exp(mpc(0,1)*ga))) + al*log(a)+be*log(b)+ga*log(c)
print("TASK 6: Mahler-measure hunt, k=3\n")
print("control: pi*m(1+x+y) =",mp.nstr(CM(1,1,1),25),"   Cl2(pi/3) =",mp.nstr(CL,25),
      "   ratio",mp.nstr(CM(1,1,1)/CL,25))
for (N,fam) in [(9,'M'),(63,'P'),(126,'M')]:
    h=Host(3,N,fam); n=mpf(h.D)**(mpf(1)/3)
    print(f"\n=== host (1{'-' if fam=='M' else '+'}{N}x)^(-1/3), D={h.D}, D^(1/3)={mp.nstr(n,12)} ===")
    for a in [1,2]:
        P=parts(h,a); B=P['Bloch']; cd=re(h.cD(a))
        print(f"  a={a}: c_D={mp.nstr(cd,25)}  Bloch={mp.nstr(B,25)}")
        # triangles with sides from this set
        S=[mpf(1),mpf(2),mpf(3),mpf(4),n,n*n,2*n,3*n,n+1,n-1,n/2,sqrt(n)]
        Snm=["1","2","3","4","r","r^2","2r","3r","r+1","r-1","r/2","sqrt(r)"]
        hits=[]
        for i,j,l in itertools.combinations_with_replacement(range(len(S)),3):
            aa,bb,cc=S[i],S[j],S[l]
            if aa>bb+cc or bb>aa+cc or cc>aa+bb: continue
            val=CM(aa,bb,cc)
            for tgt,tn in [(B,"Bloch"),(cd,"c_D")]:
                r=pslq([tgt,val,pi*log(2),pi*log(3),pi*log(n),pi**2],tol=mpf(10)**-120,maxcoeff=10**6,maxsteps=200000)
                if r and r[0]!=0 and r[1]!=0:
                    hits.append((tn,Snm[i],Snm[j],Snm[l],r))
        if hits:
            for hh in hits[:12]: print("      HIT:",hh)
        else: print("      no hit among the triangle set (PSLQ 120 digits, maxcoeff 1e6)")
        sys.stdout.flush()
