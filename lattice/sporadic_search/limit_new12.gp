default(realprecision, 80);
N=520;
a=vector(N); b=vector(N);
a[1]=1; a[2]=2; b[1]=0; b[2]=1;
{for(n=1,N-2,
  a[n+2]=((20*n^2+10*n+2)*a[n+1]-16*(2*n-1)^2*a[n])/(n+1)^2;
  b[n+2]=((20*n^2+10*n+2)*b[n+1]-16*(2*n-1)^2*b[n])/(n+1)^2);}
print("a integral? ", vecmax(vector(N,i,denominator(a[i]))));
print("b denominators max ", vecmax(vector(N,i,denominator(b[i]))));
L = b[N]/a[N]*1.0;
print("limit approx = ", L);
print("prev = ", b[N-1]/a[N-1]*1.0);
