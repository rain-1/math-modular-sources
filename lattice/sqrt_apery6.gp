default(parisizemax,2000000000); default(realprecision,400);
N=700; A=vector(N+1);B=vector(N+1);A[1]=1;A[2]=10;B[1]=0;B[2]=1;
for(n=1,N-1,A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);
xi=B[N+1]*1./A[N+1];
\\ decay check with 400-digit xi
print("log|b_n - xi a_n|/n at n=60,100: ",log(abs(B[61]-xi*A[61]))/60," ",log(abs(B[101]-xi*A[101]))/100,"  expected ",log(4*(17-12*sqrt(2))));
dn=1;ok=1;for(n=1,N,dn=lcm(dn,n);if(denominator(dn^2*B[n+1])!=1,ok=0;print("fail at ",n);break));print("d_n^2 b_n in Z to ",N,": ",ok);
default(realprecision,150);
L(D,s)=lfun(lfuncreate(D),s);
G=L(-4,2); s2=sqrt(2); l12=log(1+s2);
b1=[xi,1,zeta(2),zeta(3),G,L(8,2),L(-8,2),L(12,2),L(-3,2),L(24,2),L(-24,2),Pi*log(2),Pi*log(3),Pi*l12,Pi*s2*l12,l12^2,zeta(2)*log(2),zeta(2)*l12];
print(lindep(b1)~);
b2=[xi*s2,1,zeta(2),zeta(3),G,L(8,2),L(-8,2),L(-3,2),L(24,2),L(-24,2),Pi*log(2),Pi*l12,l12^2,zeta(2)*l12];
print(lindep(b2)~);
print("xi^2, xi*pi etc: ",lindep([xi^2,1,zeta(2),zeta(3),G,L(-8,2),L(8,2),Pi*l12])~);
\q
