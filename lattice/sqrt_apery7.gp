default(parisizemax,2000000000); default(realprecision,60);
N=400; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
\\ exact: r_n = (b_n a_N - b_N a_n)/a_N
for(k=1,5, n=40*k; r=(B[n+1]*A[N+1]-B[N+1]*A[n+1])/A[N+1]; print("n=",n,"  log|r_n|/n = ",log(abs(r))/n));
print("expected log lambda2 = ",log(4*(17-12*sqrt(2))));
print("Casoratian check n=10: ",A[11]*B[10]-A[10]*B[11]," vs binom(20,10)^2*W1: ",binomial(20,10)^2*(A[2]*B[1]-A[1]*B[2]));
\q
