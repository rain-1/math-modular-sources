default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 300);
NQ = 402; NA = 380;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
NUM = 16416*x^7 + 14148*x^6 - 15914*x^5 - 5945*x^4 - 591*x^3 - 7*x^2 + x;
DEN = (72*x^2+17*x+1)^4;
print("### the cusp^4 magnetic source on Gamma_0(6)");
print("numerator   = ", factor(NUM));
print("gcd(N,Q)    = ", gcd(NUM,DEN));
rho = NUM/DEN;
print("rho = ", NUM/gcd(NUM,DEN), " / ", DEN/gcd(NUM,DEN));
\\ Fricke: rho(1/(Cu)) =? -rho(u)
{
my(t, RI);
RI = subst(NUM, x, 1/(C*x))/subst(DEN, x, 1/(C*x));
print("rho(1/(72u)) + rho(u) = ", simplify(RI + rho));
print("rho(1/(72u)) - rho(u) = ", simplify(RI - rho));
}
Phi = F2*subst(NUM,x,us)/subst(DEN,x,us);
Xi = Dinv(Phi, NQ-2);
print("c(1..14) = ", vector(14,j,polcoeff(Phi,j)));
print("integral? first failure m = ", intfail(Phi, NQ-2), "   magnetic? first failure m = ", magfail(Phi, NQ-2));
print("|c(m)|^(1/m) at m = 100, 200, 399 : ", abs(polcoeff(Phi,100)*1.0)^(1/100), " ", abs(polcoeff(Phi,200)*1.0)^(1/200), " ", abs(polcoeff(Phi,399)*1.0)^(1/399));
Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
Bv = peel2(Fs*Th, xs, NA, NQ);
print("B_0..B_8 = ", vector(9,j,Bv[j]));
{
my(dn, kmin);
dn = 1; kmin = 0;
for(n=1, NA, dn = lcm(dn,n); my(de, kj, tt); de = denominator(Bv[n+1]); kj = 0; tt = de; while(tt>1 && kj<12, kj++; tt = tt/gcd(tt,dn)); if(kj>kmin, kmin=kj; print("   k rises to ", kj, " at n=", n)));
print("minimal k with d_n^k B_n in Z (n <= ", NA, ") : k = ", kmin);
}
print("");
print("B_n/A_n at n = 100, 200, 300, 380:");
{
my(v);
foreach([100,200,300,380], n, v = Bv[n+1]/Av[n+1]*1.0; print("   n=", n, " : ", v));
}
xi = Bv[NA+1]/Av[NA+1]*1.0;
print("");
print("digits of agreement between n=300 and n=380: ", -log(abs(Bv[381]/Av[381]-Bv[301]/Av[301])*1.0)/log(10));
print("");
print("identification attempts for xi = ", xi);
print("  lindep([xi,1])            = ", lindep([xi,1]));
print("  lindep([xi,1,zeta(2)])    = ", lindep([xi,1,zeta(2)]));
print("  lindep([xi,1,zeta(3)])    = ", lindep([xi,1,zeta(3)]));
print("  lindep([xi,1,zeta(2),zeta(3)]) = ", lindep([xi,1,zeta(2),zeta(3)]));
print("  lindep([xi,1,log(2)])     = ", lindep([xi,1,log(2)]));
print("  lindep([xi,1,log(2),log(3)]) = ", lindep([xi,1,log(2),log(3)]));
print("  algdep(xi,4)              = ", algdep(xi,4));
quit;
