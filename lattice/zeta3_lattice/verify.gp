read("/home/ubuntu/code/math-modular-sources/lattice/zeta3_lattice/rows.gp");
D=rowsI(10,4,64); T=rowsI(12,4,16); N=700; s2(n)=hammingweight(n);
vD=vector(N+1,i,valuation(D[1][i],2)); vT=vector(N+1,i,valuation(T[1][i],2));
incD(n)=6*(n-1)-3*valuation(n,2)-vD[n+1]-vD[n];
incT(n)=4*(n-1)-3*valuation(n,2)-vT[n+1]-vT[n];
tailD(m)=vecmin(vector(40,j,incD(m+j)));
tailT(m)=vecmin(vector(60,j,incT(m+j)));
pred(m,n)=vT[n+1]+vD[m+1]+min(tailD(m),2+tailT(n));
K=4000;M=2^K;
fc3=vector(N+1);fc3[1]=Mod(1,M);for(n=1,N,fc3[n+1]=fc3[n]*Mod(n,M)^3);
DA=vector(N+1,i,Mod(D[1][i],M));DC=vector(N+1,i,Mod(D[2][i],M));
TA=vector(N+1,i,Mod(T[1][i],M));TC=vector(N+1,i,Mod(T[2][i],M));
del2(m,n)={my(x=lift(3*TA[n+1]*DC[m+1]*fc3[n+1]-4*DA[m+1]*TC[n+1]*fc3[m+1]));
 if(x==0,return(oo));valuation(x,2)-3*(m-s2(m))-3*(n-s2(n))};
print("A) is pred always a LOWER bound for v2(Delta)? m,n<=60");
{my(bad=0,eqc=0,tot=0);for(m=1,60,for(n=1,60,tot++;my(a=del2(m,n),p=pred(m,n));
  if(a<p,bad++;print("  VIOLATION (",m,",",n,")"));if(a==p,eqc++)));
 print("  violations: ",bad,"/",tot,"   exact equalities: ",eqc,"/",tot);}
print("\nB) v2(a_n) bounds");
print("  v2(aT_n)=3s2(n)-(n odd) exact for n<=",N,": ",
  if(sum(n=1,N,vT[n+1]!=3*s2(n)-(n%2))==0,"YES","no"));
{my(mx=0);for(n=1,N,mx=max(mx,vD[n+1]-3*s2(n)));print("  v2(a^D_n) <= 3*s2(n)+",mx," for n<=",N);}
print("\nC) 2:3 sampling deficit  12n - v2(Delta(2n,3n))");
{my(w=vector(110,n,12*n-del2(2*n,3*n)));
 print("  range n<=110: [",vecmin(w),",",vecmax(w),"]  ; max deficit/log2(n) = ",
   vecmax(vector(110,n,w[n]/(log(n+1)/log(2))))*1.);}
print("\nD) diagonal v2(Delta(n,n)) - (4n-3)");
{print("  ",vector(10,j,del2(20*j,20*j)-(4*20*j-3)));}
print("\nE) 2-adic evidence for eps := 3*xi^D_2 - 4*xi^T_2 = 0");
{my(n=330);print("  v2(Delta(n,n)) at n=",n," is ",del2(n,n),
  " ; since v2(a^T)+v2(a^D)<=",vT[n+1]+vD[n+1],", v2(eps) >= ",del2(n,n)-vT[n+1]-vD[n+1]);}
print("\nF) closed form incT(n) = 4n-6-6*s2(n-1) ?");
{my(bad=0);for(n=2,N-1,if(incT(n)!=4*n-6-6*s2(n-1),bad++));print("  exceptions 2<=n<=",N-1,": ",bad);}
\q
