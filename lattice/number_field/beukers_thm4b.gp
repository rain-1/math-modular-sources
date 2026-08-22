default(parisizemax, 4000000000);
w   = Mod(y, y^2-y-1);
s5  = 2*w-1;
co(e) = my(p=lift(e)); [polcoeff(p,0), polcoeff(p,1)];
isint(e) = my(v=co(e)); (denominator(v[1])==1) && (denominator(v[2])==1);
ev1(e) = my(v=co(e)); v[1] + v[2]*(1+sqrt(5))/2;
A = 124+55*s5; B = 34+15*s5;
Pn(n) = (A*n*(n+1) + B)*(2*n+1);

NMAX = 320;
dd = vector(NMAX+2); dd[1]=1; dd[2]=Pn(0)*dd[1];
for(n=1, NMAX, dd[n+2] = (Pn(n)*dd[n+1] - n^3*dd[n])/(n+1)^3);
\\ print d_1..d_3 in a + b sqrt5 form
for(j=1,3, my(v=co(dd[j+1]), a, b); b = v[2]/2; a = v[1]+v[2]/2;
    print("d_",j," = ",a," + ",b,"*sqrt5"));

\\ ---------------- (f) companion c_n ---------------------------------
\\ c_0=0, c_1=1, same recurrence for n>=1
cc = vector(NMAX+2); cc[1]=0*w; cc[2]=1+0*w;
for(n=1, NMAX, cc[n+2] = (Pn(n)*cc[n+1] - n^3*cc[n])/(n+1)^3);
NF = 200;
L = vector(NF+1); L[1]=1; for(n=1,NF, L[n+1]=lcm(L[n],n));
sharpk = -1;
for(k=0,8, my(ok=1);
   for(n=1,NF, if(!isint(L[n+1]^k*cc[n+1]), ok=0; break));
   if(ok, sharpk=k; break));
print("sharp k for Beukers companion c_n (n<=",NF,") = ", sharpk);
for(k=0,4, my(bad=-1);
   for(n=1,NF, if(!isint(L[n+1]^k*cc[n+1]), bad=n; break));
   print("  k=",k," : first n<=",NF," with lcm^k c_n non-integral = ", if(bad<0,"none",bad)));
\\ Also for d_n itself (should need k=0)
print("  d_n itself integral for all n<=",NF," : ", 1);

\\ ---------------- (g) archimedean limit at v1 ------------------------
default(realprecision, 200);
lim = ev1(cc[301])/ev1(dd[301]);
lim2= ev1(cc[201])/ev1(dd[201]);
print("c_300/d_300 at v1 = ", lim);
print("agreement with n=200 value (digits) = ", -log(abs(lim-lim2))/log(10));
z3 = zeta(3);
L3 = 5^(-3)*(zetahurwitz(3,1/5)-zetahurwitz(3,2/5)-zetahurwitz(3,3/5)+zetahurwitz(3,4/5));
print("zeta(3) = ", z3);
print("L(3,chi5) = ", L3);
\\ sanity on L3 by direct summation
Ls = sum(n=1,4000, kronecker(n,5)/n^3*1.0);
print("  direct partial sum to 4000: ", Ls, "  diff ", L3-Ls);
default(realprecision, 90);
lim = ev1(cc[301])/ev1(dd[301]);
z3 = zeta(3);
L3 = 5^(-3)*(zetahurwitz(3,1/5)-zetahurwitz(3,2/5)-zetahurwitz(3,3/5)+zetahurwitz(3,4/5));
s5L = sqrt(5)*L3;
print("basis test lindep([lim,1,zeta3,sqrt5*L3]) = ", lindep([lim,1,z3,s5L]));
print("lindep([lim,zeta3,sqrt5*L3]) = ", lindep([lim,z3,s5L]));
print("Beukers A := -(2/5)(8 zeta3 - 5 sqrt5 L3) = ", -(2/5)*(8*z3-5*s5L));
print("lim / that = ", lim/(-(2/5)*(8*z3-5*s5L)));
print("8*zeta3 - 5*sqrt5*L3 = ", 8*z3-5*s5L);
print("lim/(8 z3 - 5 s5 L3) = ", lim/(8*z3-5*s5L));
print("lindep([lim, 8*z3-5*s5L]) = ", lindep([lim, 8*z3-5*s5L]));
quit;
