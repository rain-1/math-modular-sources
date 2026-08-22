default(parisizemax,8000000000);
N=1200;
A=vector(N+2); A[1]=1; A[2]=5;
{for(n=1,N, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
Fs = sum(n=0,N, A[n+1]*t^n) + O(t^(N+1));
sig = sqrt(1-34*t+t^2+O(t^(N+1)));
H = Fs*sig;
bad=0; {for(i=0,N, if(denominator(polcoeff(H,i))!=1, bad++));}
print("F*sigma integral to t^",N,"? bad=",bad);
R = 1/H;
r = vector(N+1, i, polcoeff(R,i-1));
bad=0; {for(i=1,N+1, if(denominator(r[i])!=1, bad++));}
print("R=1/(F sigma) integral? bad=",bad);
print("r_0..r_8 = ", vector(9,i,r[i]));
\\ Gauss congruences r_{m p^a} = r_{m p^{a-1}} mod p^a
fails=0; tested=0;
{forprime(p=2,60, my(pa,a); a=1; pa=p;
  while(pa<=N,
    for(m=1,N\pa, my(n=m*pa, nn=m*pa\p);
      tested++;
      if((r[n+1]-r[nn+1])%(p^a)!=0, fails++; if(fails<6,print("FAIL p=",p," a=",a," m=",m))));
    a++; pa*=p));}
print("Gauss congruences: tested ",tested," instances (p<=60, n<=",N,"), failures = ",fails);
quit;
