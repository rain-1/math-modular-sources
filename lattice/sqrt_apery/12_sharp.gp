default(parisizemax,6000000000);
M=60; Ap=vector(M);Ap[1]=1;Ap[2]=5;
{for(n=1,M-2,Ap[n+2]=((2*n+1)*(17*n^2+17*n+5)*Ap[n+1]-n^3*Ap[n])/(n+1)^3);}
s=Ser(Ap,x)+O(x^M); r=sqrt(s);
print("2^n [t^n] sqrt(F), n=0..10: ",vector(11,i,2^(i-1)*polcoeff(r,i-1)));
print("v_2 of [t^n]sqrt(F) den, n=0..30: ",vector(31,i,valuation(denominator(polcoeff(r,i-1)),2)));
print("a_n mod 4, n=0..12: ",vector(13,i,lift(Mod(4^(i-1)*polcoeff(r,i-1),4))));
\q
