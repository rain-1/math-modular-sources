default(parisizemax,12000000000);
read("lib.gp");
NQ = 205; MB = 200;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
print("### where magnetism dies on Apery's host (N=6, C=72, B=17)");
print("A_n = [x^n]F, x = u/(1+17u+72u^2):");
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, 40, NQ);
print("A_0..A_10 = ", vector(11,j,Av[j]));
print("v_p(A_n) - v_p(n+1) for the first failures of (n+1)|A_n:");
{
for(n=1, 30,
  if(Av[n+1] % (n+1) != 0, print("   n=", n, "  A_n=", Av[n+1], "   (n+1)=", n+1, "  gcd=", gcd(Av[n+1],n+1)); );
);
}
print("");
print("### diagnostic: Cooper-shape rho = u(1-72u^2)/(1+B'u+72u^2)^2, list of m<=60 with m ndiv c(m)");
{
for(Bp = 14, 20,
  my(rho, Phi, bad);
  rho = us*(1-C*us^2)/(1+Bp*us+C*us^2)^2;
  Phi = F2*rho;
  bad = List();
  for(m=1, 60, if(polcoeff(Phi,m) % m != 0, listput(bad, [m, gcd(polcoeff(Phi,m),m)])));
  print("  B'=", Bp, "  #bad(m<=60)=", #bad, "   first 12: ", vector(min(12,#bad), j, bad[j]));
);
}
print("");
print("### the 2-adic and 3-adic profile of Xi for B'=16 (the disc -8 pole)");
{
my(rho, Phi, Xi);
rho = us*(1-C*us^2)/(1+16*us+C*us^2)^2;
Phi = F2*rho;
print("  c(m), m=1..24: ", vector(24,j,polcoeff(Phi,j)));
print("  m/gcd(c(m),m), m=1..48: ", vector(48,j,my(cm=polcoeff(Phi,j)); j/gcd(cm,j)));
}
quit;
