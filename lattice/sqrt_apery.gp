default(realprecision,120);
\\ Sym^1 of Apery: (n+1)^2 u_{n+1} = (136 n^2+136 n+10) u_n - 4 n^2 u_{n-1}
N=400; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+136*n+10)*A[n+1]-4*n^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+136*n+10)*B[n+1]-4*n^2*B[n])/(n+1)^2);
print("a_n: ",vector(8,i,A[i]));
print("integral to N? ",vector(N+1,i,denominator(A[i]))==vector(N+1,i,1));
\\ check sqrt relation with Apery zeta(3) generating function: A_n apery
Ap=vector(60);Ap[1]=1;Ap[2]=5;for(n=1,58,Ap[n+2]=((2*n+1)*(17*n^2+17*n+5)*Ap[n+1]-n^3*Ap[n])/(n+1)^3);
s=Ser(Ap,x); r=sqrt(s); print("4^n [x^n] sqrt(F_Apery): ",vector(8,n,4^(n-1)*polcoeff(r,n-1)));
lim=B[N+1]*1./A[N+1]; print("limit=",lim);
print("rate of |b-lim a|: ",log(abs(B[N+1]-lim*A[N+1]))/N, "  log a_n/n=",log(A[N+1])/N);
dn(n)=lcm(vector(n,i,i)); print("d_n^2 b_n integral for n<=200? ",vector(200,n,denominator(dn(n)^2*B[n+1]))==vector(200,n,1));
print("d_n b_n integral? ",vector(60,n,denominator(dn(n)*B[n+1]))==vector(60,n,1));
\\ identification attempts
G=sumalt(k=0,(-1)^k/(2*k+1)^2);
basis=[1,zeta(2),zeta(3),G,Pi^2,Pi*log(2),log(2)^2,Pi*sqrt(2)*log(1+sqrt(2)),Pi*sqrt(3),sqrt(2)*zeta(2),sqrt(3)*zeta(2),Pi*log(1+sqrt(2)),log(1+sqrt(2))^2,Pi*sqrt(2),Pi*sqrt(6)];
print(lindep(concat([lim],basis))~);
\q
