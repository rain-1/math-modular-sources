M=120; default(seriesprecision,M+5);
et(k)=eta(q^k+O(q^(M+5)));
u = 3*q*et(2)^2*et(12)^4/(et(4)^4*et(6)^2);
ps1 = et(2)^6*et(3)/(et(1)^3*et(6)^2);
F12 = ps1*(3-4*u-3*u^2)/(3*(1+u)^2);
t = u/3;
tt = serreverse(t);
Ft = subst(F12, q, tt);
A = vector(M, n, polcoeff(Ft, n-1));
print("A_n: ",vector(14,i,A[i]));
print("integral? ",vector(M,i,denominator(A[i]))==vector(M,i,1));
d=2; rowsM=[];
{for(n=1,M-2, my(r=[]); for(j=0,2, for(e=0,d, r=concat(r,[n^e*A[n+2-j]]))); rowsM=concat(rowsM,[r]));}
Mx=matrix(#rowsM,3*(d+1),i,j,rowsM[i][j]);
K=matker(Mx); print("kernel dim: ",#K);
if(#K>0, v=K[,1]; v=v/content(v); print("p0,p1,p2 coeffs (1,n,n^2): ",vector(3,j,vector(d+1,e,v[(j-1)*(d+1)+e]))));
\q
