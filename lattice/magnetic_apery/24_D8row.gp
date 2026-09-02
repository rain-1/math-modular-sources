default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 500);
NQ = 262; NA = 250;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
rho = subst(-x*(72*x^2+18*x+1)*(72*x^2-1), x, us)/subst((72*x^2+16*x+1)^2*(72*x^2+17*x+1), x, us);
Phi = F2*rho;
Xi = Dinv(Phi, NQ-2);
Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
Bv = peel2(Fs*Th, xs, NA, NQ);
print("### the D8 source on Apery's host: precise convergence data");
lam1 = 17+12*sqrt(2); lam2 = 17-12*sqrt(2);
print("lam1 = ", lam1, "  lam2 = ", lam2, "  lam1*lam2 = ", lam1*lam2);
print("");
print("identification of xi = lim B_n/A_n:");
xi = Bv[NA+1]/Av[NA+1]*1.0;
print("  B_250/A_250 = ", xi);
print("  zeta(2)/8   = ", zeta(2)/8);
print("  difference  = ", xi - zeta(2)/8);
print("  lindep([xi, zeta(2)]) = ", lindep([xi, zeta(2)]));
print("");
print("n | B_n/A_n - zeta(2)/8 | ratio to previous | A_n*zeta(2)/8 - B_n");
{
my(prev, val, lf);
prev = 0;
forstep(n=20, NA, 10,
  val = Bv[n+1]/Av[n+1] - zeta(2)/8;
  lf = Av[n+1]*zeta(2)/8 - Bv[n+1];
  print("  ", n, " | ", val*1.0, " | ", if(prev==0, "-", (val/prev)^(1/10)*1.0), " | ", lf*1.0);
  prev = val;
);
}
print("");
print("the linear form L_n = A_n zeta(2)/8 - B_n :  L_n * n^(3/2) at n = 50,100,150,200,250:");
{
forstep(n=50, NA, 50, print("  n=", n, "  L_n = ", (Av[n+1]*zeta(2)/8 - Bv[n+1])*1.0, "   L_n*n^(3/2) = ", (Av[n+1]*zeta(2)/8 - Bv[n+1])*n^1.5));
}
print("");
print("d_n^2 |L_n| (the quantity that must -> 0 for irrationality):");
{
my(dn);
dn = 1;
for(n=1, NA, dn = lcm(dn,n); if(n%25==0, print("  n=", n, "   d_n^2 |L_n| = ", abs(Av[n+1]*zeta(2)/8 - Bv[n+1])*dn^2*1.0)));
}
quit;
