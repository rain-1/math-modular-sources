default(parisizemax,10000000000);
N=4000;
A=vector(N+2); A[1]=1; A[2]=5;
{for(n=1,N, A[n+2] = ((2*n+1)*(17*n^2+17*n+5)*A[n+1] - n^3*A[n])/(n+1)^3);}
print("A computed");
Fs = sum(n=0,N, A[n+1]*t^n) + O(t^(N+1));
sig = sqrt(1-34*t+t^2+O(t^(N+1)));
b=0; {for(i=0,N, if(denominator(polcoeff(sig,i))!=1, b++));}
print("sigma in Z[[t]] to t^",N,": violations=",b);
H = Fs*sig;
b=0; {for(i=0,N, if(denominator(polcoeff(H,i))!=1, b++));}
print("F*sigma in Z[[t]] to t^",N,": violations=",b);
R = 1/H; r = vector(N+1,i,polcoeff(R,i-1));
b=0; {for(i=1,N+1, if(denominator(r[i])!=1, b++));}
print("1/(F sigma) in Z[[t]]: violations=",b);
fails=0; tested=0;
{forprime(p=2,200, my(pa=p,a=1);
  while(pa<=N,
    for(m=1,N\pa, my(n=m*pa); tested++;
      if((r[n+1]-r[n\p+1])%(p^a)!=0, fails++; if(fails<6,print("GAUSS FAIL p=",p," a=",a," m=",m))));
    a++; pa*=p));}
print("Gauss congruences r_{mp^a}=r_{mp^{a-1}} mod p^a: tested=",tested," fails=",fails," (p<=200, n<=",N,")");
\\ q integrality directly
u = sum(n=1,N, r[n+1]/n*t^n) + O(t^(N+1));
q = t*exp(u);
b=0; {for(i=0,N, if(denominator(polcoeff(q,i))!=1, b++));}
print("q(t) in t Z[[t]] to t^",N,": violations=",b);
tq = serreverse(q);
b=0; {for(i=0,N, if(denominator(polcoeff(tq,i))!=1, b++));}
print("t(q) in q Z[[q]] to q^",N,": violations=",b);
Fq = subst(Fs,t,tq);
b=0; {for(i=0,N, if(denominator(polcoeff(Fq,i))!=1, b++));}
print("F(t(q)) in 1+q Z[[q]] to q^",N,": violations=",b);
print("t(q) coeffs 1..12: ", vector(12,i,polcoeff(tq,i)));
print("F(t(q)) coeffs 0..12: ", vector(13,i,polcoeff(Fq,i-1)));
quit;
