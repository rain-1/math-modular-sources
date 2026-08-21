default(realprecision, 200);
N=400;
a=vector(N); b=vector(N); a[1]=1; a[2]=2; b[1]=0; b[2]=1;
{for(n=1,N-2,
  a[n+2]=((20*n^2+10*n+2)*a[n+1]-16*(2*n-1)^2*a[n])/(n+1)^2;
  b[n+2]=((20*n^2+10*n+2)*b[n+1]-16*(2*n-1)^2*b[n])/(n+1)^2);}
L=b[N]/a[N]*1.0;
print("L=",L);
print("alg2 ", lindep([1,L,L^2]));
print("alg3 ", lindep([1,L,L^2,L^3]));
print("z2 ", lindep([1,zeta(2),L]));
print("z3 ", lindep([1,zeta(3),L]));
print("G  ", lindep([1,Catalan,L]));
print("l2 ", lindep([1,log(2),L]));
print("pi ", lindep([1,Pi,L]));
print("pil2", lindep([1,Pi*log(2),L]));
print("sq2 ", lindep([1,sqrt(2),L]));
print("big ", lindep([1,zeta(2),zeta(3),Catalan,Pi^3,log(2),Pi*log(2),L]));
