default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 400);
NQ = 262; NA = 250;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
print("### Apery host row: A_n = ", vector(8,j,Av[j]));
NAMES = ["FOLD", "CUSPPOLE", "D8", "D12"];
{RHO = [(-10368*x^6-4896*x^5-864*x^4-68*x^3-2*x^2)/(72*x^2-1)^4,
       (25920*x^5+6048*x^4-84*x^2-5*x)/(72*x^2+17*x+1)^3,
       (-x*(72*x^2+18*x+1)*(72*x^2-1))/((72*x^2+16*x+1)^2*(72*x^2+17*x+1)),
       (-x)/((12*x+1)^2*(9*x+1))];}
{
for(t=1, 4,
  my(rho, Phi, Xi, Th, Bv, dn, kmin, ok, ratio, lam2);
  print("");
  print("== source ", NAMES[t]);
  rho = subst(numerator(RHO[t]), x, us)/subst(denominator(RHO[t]), x, us);
  Phi = F2*rho;
  Xi = Dinv(Phi, NQ-2);
  \\ Theta = D^{-2} Xi = sum c'(m) m^{-2} q^m
  Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
  Bv = peel2(Fs*Th, xs, NA, NQ);
  print("   B_0..B_6 = ", vector(7,j,Bv[j]));
  dn = 1; kmin = 0;
  for(n=1, NA, dn = lcm(dn, n); my(de, kj, tt); de = denominator(Bv[n+1]); kj = 0; tt = de; while(tt>1 && kj<12, kj++; tt = tt/gcd(tt,dn)); if(kj>kmin, kmin=kj));
  print("   minimal k with d_n^k B_n in Z (n <= ", NA, ") : k = ", kmin);
  print("   B_n/A_n at n = 60, 120, 250:");
  print("      ", Bv[61]/Av[61]*1.0);
  print("      ", Bv[121]/Av[121]*1.0);
  print("      ", Bv[NA+1]/Av[NA+1]*1.0);
  my(df); df = abs(Bv[NA+1]/Av[NA+1] - Bv[121]/Av[121])*1.0; print("   agreement of the last two (digits): ", if(df==0, "exact", -log(df)/log(10)));
);
}
quit;
