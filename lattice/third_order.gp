default(realprecision,60);
\\ (n+1)^3 u_{n+1} = (2n+1)(a n^2 + a n + b) u_n - c n^3 u_{n-1}
nm=["AZ(7,3,81)","(9,3,-27)","Domb(10,4,64)","(11,5,125)","(12,4,16)","Apery(17,5,1)"];
par=[[7,3,81],[9,3,-27],[10,4,64],[11,5,125],[12,4,16],[17,5,1]];
N=400;
rows(a,b,c)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=b;B[1]=0;B[2]=1;
 for(n=1,N-1,A[n+2]=((2*n+1)*(a*n^2+a*n+b)*A[n+1]-c*n^3*A[n])/(n+1)^3;
            B[n+2]=((2*n+1)*(a*n^2+a*n+b)*B[n+1]-c*n^3*B[n])/(n+1)^3);[A,B]};
R=vector(6,i,rows(par[i][1],par[i][2],par[i][3]));
z3=zeta(3); L3m3=lfun(lfuncreate(-3),3); L3m4=Pi^3/32; 
basis=[1,z3,L3m3,L3m4,Pi^2,Pi^2*log(2),Pi^2*log(3),Pi*sqrt(3)*zeta(2),log(2)^3];
{for(i=1,6,my(A=R[i][1],B=R[i][2],c=par[i][3]);
 my(l1=B[N+1]*1./A[N+1],l2=B[N]*1./A[N]);
 my(sl=vector(4,j,my(p=[2,3,5,7][j]);valuation(B[N+1]/A[N+1]-B[N]/A[N],p)*1./N));
 my(mv=vector(4,j,my(p=[2,3,5,7][j]);vecmax(vector(N,n,valuation(A[n+1],p))))); 
 print(nm[i],"  c=",c," roots=",polroots(x^2-2*par[i][1]*x+c)~);
 print("   lim=",l1,"  conv=",abs(l1-l2));
 print("   lindep[1,z3,L(3,chi-3),pi^3/32,pi^2,pi^2log2,pi^2log3,pi sqrt3 z2,log^3 2]=",lindep(concat([l1],basis))~);
 print("   slopes p=2,3,5,7: ",sl,"   max v_p(a_n): ",mv));}
\q
