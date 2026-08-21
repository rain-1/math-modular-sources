default(realprecision,120);
N=300;
rowsLf(N)={my(A=vector(N+1),B=vector(N+1));A[1]=1;A[2]=2;B[1]=0;B[2]=1;
 for(n=1,N-1,
   my(P=20*n^2+10*n+2, Q=16*(2*n-1)^2);
   A[n+2]=(P*A[n+1]-Q*A[n])/(n+1)^2;
   B[n+2]=(P*B[n+1]-Q*B[n])/(n+1)^2);
 [A,B]};
R=rowsLf(N); A=R[1]; B=R[2];
print("v2(A_n) at n=50,100,150,200,250,300: ", vector(6,i,valuation(A[50*i+1],2)));
print("A_n integral? denom check n=1..20: ", vector(20,i,denominator(A[i])));
for(k=1,6, my(n=50*k); print("n=",n," v2(diff)=",valuation(B[n+1]/A[n+1]-B[n]/A[n],2)));
\q
