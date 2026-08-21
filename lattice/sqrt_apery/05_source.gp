default(parisizemax,6000000000);
M=120; default(seriesprecision,M+6);
et(k)=eta(q^k+O(q^(M+6)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5);
t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
f=sqrt(F); Dq(g)=q*deriv(g,q);
\\ (0) identity  Dq t = t*F*sqrt(1-34t+t^2)
P3 = 1-34*t+t^2;
print("[chk] Dq(t) - t*F*sqrt(P3) == 0 ? ", Dq(t)-t*F*sqrt(P3)==0);
\\ (1) the source Psi = f^3 * t/4
Psi = f^3*t/4;
print("psi(m), m=1..8: ",vector(8,i,polcoeff(Psi,i)));
print("v_2(4^m psi(m)) m=1..40 (>=0 required): ",vector(40,m,valuation(4^m*polcoeff(Psi,m),2)));
print("4^m psi(m) integral for m<=",M-2,"? ",vector(M-2,m,denominator(4^m*polcoeff(Psi,m)))==vector(M-2,m,1));
print("4^m [q^m] f integral? ",vector(M,m,denominator(4^m*polcoeff(f,m-1)))==vector(M,m,1));
Theta = sum(m=1,M, polcoeff(Psi,m)/m^2*q^m) + O(q^(M+1));
qt = serreverse(t);
print("[chk] q as series in t integral? ",vector(M,i,denominator(polcoeff(qt,i)))==vector(M,i,1));
G = subst(f*Theta, q, qt);
bb = vector(M-4, n, 4^n*polcoeff(G,n));
N=M; A=vector(N+2); B=vector(N+2); A[1]=1;A[2]=10;B[1]=0;B[2]=1;
{for(n=1,N, A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;
            B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);}
print("Eichler b_n n=1..8: ",vector(8,i,bb[i]));
print("recur  b_n n=1..8: ",vector(8,i,B[i+1]));
print("MATCH to n=",M-4,"? ", vector(M-4,n,bb[n])==vector(M-4,n,B[n+1]));
\\ (2) e_{n,m} = 4^n [t^n](f q^m) integral and divisible by 4^m
ft = subst(f,q,qt);
ok=1; okd=1;
{for(n=1,40, my(S=ft*qt^0); for(m=1,n, my(e=4^n*polcoeff(ft*qt^m,n));
   if(denominator(e)!=1, ok=0; print("nonint n=",n," m=",m));
   if(e%4^m!=0, okd=0; print("not div by 4^m: n=",n," m=",m," v2=",valuation(e,2)," need ",2*m))));}
print("e_{n,m} integral (n<=40): ",ok,"   divisible by 4^m: ",okd);
\q
