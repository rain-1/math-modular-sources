default(parisizemax,12000000000);
read("lib.gp");
NQ = 1005; MC = 1000;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
rho = subst(-x*(72*x^2+18*x+1)*(72*x^2-1), x, us)/subst((72*x^2+16*x+1)^2*(72*x^2+17*x+1), x, us);
Phi = F2*rho;
cp = vector(MC, m, polcoeff(Phi,m)/m);
print("D8 source, p = 2: v_2(gcd of c'(2^i m) - eps^i c'(m)) for eps = +1, -1, and c'(3m)/c'(m):");
{
for(s=1, 2,
  my(eps, g); eps = if(s==1, 1, -1); g = 0;
  for(i=1, 3, my(pi); pi = 2^i; for(m=1, MC\pi, if(m%2==0, next); g = gcd(g, cp[pi*m] - eps^i*cp[m])));
  print("   eps=", eps, "  e = ", if(g==0,"inf",valuation(g,2)));
);
}
print("exact U_3 eigen relation?  c'(3m) - c'(m) for m=1..20 : ", vector(20,m,if(m%3==0,"-",cp[3*m]-cp[m])));
print("c'(1..16) = ", vector(16,m,cp[m]));
quit;
