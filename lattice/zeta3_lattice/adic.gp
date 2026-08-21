read("/home/ubuntu/code/math-modular-sources/lattice/zeta3_lattice/rows.gp");
D=rowsI(10,4,64); T=rowsI(12,4,16);
N=700; s2(n)=hammingweight(n);
vD=vector(N+1,i,valuation(D[1][i],2)); vT=vector(N+1,i,valuation(T[1][i],2));
print("== v2(a_n) ==");
print("max v2(a^D_n) n<=",N,": ",vecmax(vD),"   max v2(a^T_n): ",vecmax(vT));
{my(bad=0);for(n=1,N, if(vT[n+1]!=3*s2(n)-(n%2),bad++));print("v2(a^T_n)=3s2(n)-(n mod 2): exceptions n<=",N,": ",bad);}
{my(mx=-99,mn=99);for(n=1,N, my(d=vD[n+1]-3*s2(n)); mx=max(mx,d);mn=min(mn,d));
 print("v2(a^D_n)-3s2(n) in [",mn,",",mx,"]");}
print("v2(a^D_n) <= 4*log2(n) for 2<=n<=",N,"? ", if(vecmin(vector(N-1,j,4*log(j+1)/log(2)-vD[j+2]))>=0,"yes","NO"));
incD(n)=6*(n-1)-3*valuation(n,2)-vD[n+1]-vD[n];
incT(n)=4*(n-1)-3*valuation(n,2)-vT[n+1]-vT[n];
tailD(m)=vecmin(vector(40,j,incD(m+j)));
tailT(m)=vecmin(vector(60,j,incT(m+j)));
print("\n== tails ==");
{for(j=1,6,my(m=25*j);print("m=",m," tailD=",tailD(m)," def=",6*m-tailD(m),"   tailT=",tailT(m)," def=",4*m-tailT(m)));}
{my(dD=vector(400,m,6*m-tailD(m)),dT=vector(400,m,4*m-tailT(m)));
 print("deficit 6m-tailD in [",vecmin(dD),",",vecmax(dD),"];  4m-tailT in [",vecmin(dT),",",vecmax(dT),"]");}
\\ ---- cross determinant, computed 2-adically mod 2^K ----
K=4000; M=2^K;
fc3=vector(N+1); fc3[1]=Mod(1,M); for(n=1,N, fc3[n+1]=fc3[n]*Mod(n,M)^3);
DA=vector(N+1,i,Mod(D[1][i],M)); DC=vector(N+1,i,Mod(D[2][i],M));
TA=vector(N+1,i,Mod(T[1][i],M)); TC=vector(N+1,i,Mod(T[2][i],M));
del2(m,n)={my(x=lift(3*TA[n+1]*DC[m+1]*fc3[n+1] - 4*DA[m+1]*TC[n+1]*fc3[m+1]));
  if(x==0,return(oo)); valuation(x,2)-3*(m-s2(m))-3*(n-s2(n))};
pred(m,n)=vT[n+1]+vD[m+1]+min(tailD(m),2+tailT(n));
print("\n== cross determinant Delta(m,n)=3 aT_n bD_m - 4 aD_m bT_n ==");
{my(bad=0,tot=0);for(m=1,50,for(n=1,50, tot++; if(del2(m,n)!=pred(m,n), bad++; if(bad<8,print("  mismatch (",m,",",n,") act=",del2(m,n)," pred=",pred(m,n))))));
 print("exhaustive m,n<=50 mismatches: ",bad,"/",tot);}
print("\n== 2:3 sampling: v2(Delta(2n,3n)) - 12n ==");
{my(w=vector(100,n,del2(2*n,3*n)-12*n));print("range over n<=100: [",vecmin(w),",",vecmax(w),"]");
 for(j=1,8,my(n=12*j);print("n=",n,"  v2=",del2(2*n,3*n),"  -12n=",del2(2*n,3*n)-12*n));}
print("\n== v2(Delta(a n, c n))/n for various a:c at n=40 ==");
{for(a=1,4,for(c=1,6, print1("  ",a,":",c,"->",del2(a*40,c*40)/40.,"")));print();}
\q
