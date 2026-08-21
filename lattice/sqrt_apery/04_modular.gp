default(parisizemax,6000000000);
M=90; default(seriesprecision,M+6);
et(k)=eta(q^k+O(q^(M+6)));
F=et(2)^7*et(3)^7/(et(1)^5*et(6)^5);
t=q*(et(1)*et(6)/(et(2)*et(3)))^12;
print("F: ",vector(8,i,polcoeff(F,i-1)));
print("t: ",vector(8,i,polcoeff(t,i)));
f=sqrt(F);
print("f=sqrt(F) coeffs q^0..q^12: ",vector(13,i,polcoeff(f,i-1)));
print("v_2 of den(f_m), m=0..40: ",vector(41,i,valuation(denominator(polcoeff(f,i-1)),2)));
Dq(g)=q*deriv(g,q);
Phi = f*Dq(t)/4;
print("Phi=f*Dq(t)/4, c(m) m=1..12: ",vector(12,i,polcoeff(Phi,i)));
print("v_2 den c(m), m=1..40: ",vector(40,i,valuation(denominator(polcoeff(Phi,i)),2)));
\\ Theta = Dq^{-2} Phi
Theta = sum(m=1,M, polcoeff(Phi,m)/m^2*q^m) + O(q^(M+1));
\\ b_n =? [u^n] (f*Theta),  u=t/4  => [u^n]X = 4^n [t^n]X
G = f*Theta;
\\ express in t : reversion
qt = serreverse(t);      \\ q as series in t
Gt = subst(G, q, qt);
bb = vector(M-2, n, 4^n*polcoeff(Gt,n));
print("candidate b_n n=1..8: ",vector(8,i,bb[i]));
N=M; A=vector(N+2); B=vector(N+2); A[1]=1;A[2]=10;B[1]=0;B[2]=1;
{for(n=1,N, A[n+2]=((136*n^2+68*n+10)*A[n+1]-4*(2*n-1)^2*A[n])/(n+1)^2;
            B[n+2]=((136*n^2+68*n+10)*B[n+1]-4*(2*n-1)^2*B[n])/(n+1)^2);}
print("true    b_n n=1..8: ",vector(8,i,B[i+1]));
print("ratio: ",vector(8,i,bb[i]/B[i+1]));
\\ also check A_n from f
ft=subst(f,q,qt); print("4^n[t^n]f  n=0..8: ",vector(9,i,4^(i-1)*polcoeff(ft,i-1)));
\q
