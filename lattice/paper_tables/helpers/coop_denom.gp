default(realprecision,60);
f10(n,u1,u0)=2*(2*n+1)*(3*n^2+3*n+1)*u1 + 4*n*(16*n^2-1)*u0;
rowsC(fn,b,N)={my(A=vector(N+1));A[1]=1;A[2]=b;
 for(n=1,N-1,A[n+2]=fn(n,A[n+1],A[n])/(n+1)^3);A};
A=rowsC(f10,2,300);
print("denom A_n s10, n=0..20: ", vector(21,i,denominator(A[i])));
print("denom A_300: ", denominator(A[301]));
print("v2(denom A_300)=",valuation(denominator(A[301]),2)," v3=",valuation(denominator(A[301]),3));
\q
