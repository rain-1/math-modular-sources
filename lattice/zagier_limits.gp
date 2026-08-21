default(realprecision,80);
par=[[7,2,-8],[9,3,27],[10,3,9],[11,3,-1],[12,4,32],[17,6,72]]; fam=["A","B","C","D","E","F"];
N=600;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((a*n^2+a*n+b)*A[n+1]-c*n^2*A[n])/(n+1)^2;B[n+2]=((a*n^2+a*n+b)*B[n+1]-c*n^2*B[n])/(n+1)^2);[A,B]};
G=sumalt(k=0,(-1)^k/(2*k+1)^2); L3=lfun(lfuncreate(-3),2); z2=zeta(2);
basis=[1,z2,G,L3,Pi*log(2),Pi*sqrt(3)*log(3), log(2)^2,Pi^2*log(2)];
{for(i=1,6,my(R=rows(par[i][1],par[i][2],par[i][3]),A=R[1],B=R[2]);
 my(l1=B[N+1]*1./A[N+1],l2=B[N]*1./A[N]);
 print(fam[i],"  lim=",l1,"  |diff consecutive|=",abs(l1-l2));
 print("   lindep [1,z2,G,L(2,chi-3),pi log2,pi sqrt3 log3,log^2 2, pi^2 log2]: ",lindep(concat([l1],basis))~));}
\q
