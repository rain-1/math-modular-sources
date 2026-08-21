default(parisizemax,2000000000); default(realprecision,200);
N=600; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
print("a: ",vector(6,i,A[i]));
lim=B[N+1]*1./A[N+1]; lim2=B[N]*1./A[N]; print("limit=",lim); print("conv=",abs(lim-lim2));
print("log|b-lim a|/n at n=300: ",log(abs(B[301]-lim*A[301]))/300,"  log a/n: ",log(A[301])/300, "  log lambda2=",log(0.11774900609143765751894));
dn=1;k1=1;k2=1;for(n=1,300,dn=lcm(dn,n);if(denominator(dn*B[n+1])!=1,k1=0);if(denominator(dn^2*B[n+1])!=1,k2=0));
print("d_n b_n in Z: ",k1,"  d_n^2 b_n in Z: ",k2);
\\ what is the minimal denominator of b_n? compare with d_n
print("den(b_n) for n=10..15: ",vector(6,i,denominator(B[10+i])), "  d_n: ",vector(6,i,lcm(vector(9+i,j,j))));
print("score with k=1: ",-log(0.11774900609143765751894)-1, "  with k=2: ",-log(0.11774900609143765751894)-2);
G=sumalt(k=0,(-1)^k/(2*k+1)^2); L8=lfun(lfuncreate(-8),2); L3=lfun(lfuncreate(-3),2); L24=lfun(lfuncreate(-24),2);
print(lindep([lim,1,zeta(2),zeta(3),G,L8,L3,L24,Pi*log(2),Pi*log(1+sqrt(2)),Pi*sqrt(2),Pi*sqrt(2)*log(1+sqrt(2)),log(1+sqrt(2))^2,Pi^2*sqrt(2)])~);
\q
