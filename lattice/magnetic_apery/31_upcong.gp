default(parisizemax,12000000000);
read("lib.gp");
NQ = 1005; MC = 1000;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
rho = subst(-x*(72*x^2+18*x+1)*(72*x^2-1), x, us)/subst((72*x^2+16*x+1)^2*(72*x^2+17*x+1), x, us);
Phi = F2*rho;
cp = vector(MC, m, polcoeff(Phi,m)/m);
print("### U_p congruences for the D8 source  Phi = -F^2 u (1+18u+72u^2)(72u^2-1)/((1+16u+72u^2)^2 (1+17u+72u^2))");
print("### (= minus the canonical source of the N=6, C=81, B=14 host)");
print("test c'(p^i m) = psi(p)^i c'(m) mod p^e, p ndiv m, i=1,2,3 ; psi = 1");
print("p | e = v_p(gcd of all defects)  [inf = identically 0]");
{
forprime(p=2, 43,
  my(g, ntests);
  g = 0; ntests = 0;
  for(i=1, 3,
    my(pi); pi = p^i;
    if(pi > MC, break);
    for(m=1, MC\pi, if(m%p==0, next); g = gcd(g, cp[pi*m] - cp[m]); ntests++);
  );
  print("  p=", p, "  e = ", if(g==0, "inf", valuation(g,p)), "   (", ntests, " tests)");
);
}
print("");
print("strong p-magnetic: p^n | m  =>  p^n | c(m) ?");
{
forprime(p=2, 13,
  my(ok); ok = 1;
  for(m=1, MC, my(v); v = valuation(m,p); if(v>0 && valuation(polcoeff(Phi,m),p) < v, ok=0; print("   FAILS p=",p," m=",m); break));
  if(ok, print("   p=", p, " : holds for all m <= ", MC));
);
}
quit;
