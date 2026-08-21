default(parisizemax,2000000000); default(realprecision,150);
M=220; Ap=vector(M);Ap[1]=1;Ap[2]=5;for(n=1,M-2,Ap[n+2]=((2*n+1)*(17*n^2+17*n+5)*Ap[n+1]-n^3*Ap[n])/(n+1)^3);
s=Ser(Ap,x)+O(x^M); r=sqrt(s); a=vector(M,n,4^(n-1)*polcoeff(r,n-1));
print("integral? ",vector(M,i,denominator(a[i]))==vector(M,i,1));
\\ fit (n+1)^2 u_{n+1} = (A n^2 + A n + Bc) u_n - C n^2 u_{n-1}  (Zagier form) or general deg-2 3-term
rows=[];for(n=1,M-2,rows=concat(rows,[[n^2*a[n+2],n*a[n+2],a[n+2],n^2*a[n+1],n*a[n+1],a[n+1],n^2*a[n],n*a[n],a[n]]]));
Mx=matrix(#rows,9,i,j,rows[i][j]); K=matker(Mx); print("kerdim ",#K); v=K[,1]; v=v/content(v); print("p0,p1,p2 (n^2,n,1): ",vector(3,j,vector(3,e,v[(j-1)*3+e])));
\\ now companion with this recurrence, normalised b_0=0,b_1=1
p0(n)=v[1]*n^2+v[2]*n+v[3]; p1(n)=v[4]*n^2+v[5]*n+v[6]; p2(n)=v[7]*n^2+v[8]*n+v[9];
N=400; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=a[2];B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=-(p1(n)*A[n+1]+p2(n)*A[n])/p0(n);B[n+2]=-(p1(n)*B[n+1]+p2(n)*B[n])/p0(n));
print("A matches a? ",vector(M,i,A[i])==a);
lim=B[N+1]*1./A[N+1]; print("limit=",lim);
print("char roots: ",polroots(v[1]*x^2+v[4]*x+v[7])~);
print("log|b-lim a|/n: ",log(abs(B[N+1]-lim*A[N+1]))/N,"   log a/n: ",log(A[N+1])/N);
dn=1;ok2=1;ok1=1;for(n=1,200,dn=lcm(dn,n);if(denominator(dn^2*B[n+1])!=1,ok2=0);if(denominator(dn*B[n+1])!=1,ok1=0));
print("d_n^2 b_n in Z (n<=200): ",ok2,"   d_n b_n in Z: ",ok1);
\q
