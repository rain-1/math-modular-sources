default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 250);
NQ = 402; NA = 380; M2 = 150; M3 = 380;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
print("### the W_6 = -1 magnetic sources with poles only at the cusps 1/2, 1/3 (x = infinity)");
print("### rho = v P(w)/(72w+17)^d in u-form: rho = 72 u (72u^2-1) Pn(u) / (72^d (72u^2+17u+1)^d)");
{
for(d=3, 8,
  my(Qw, qv, Xs, T, S);
  Qw = (72*w+17)^d;
  qv = vector(poldegree(Qw)+1, t, polcoeff(Qw, t-1));
  Xs = xibasis(C, us, F2, qv, M3);
  T = xitomat(Xs, M3);
  S = magsub(T, M2, M3);
  print("");
  print("== d = ", d, "   anti-invariant dim = ", #Xs, "   magnetic vectors: ", #S);
  for(i=1, #S,
    my(al, Phi, Xi, Th, Bv, dn, kmin, xi, dg, PN);
    al = S[i];
    PN = sum(k=0, #Xs-1, al[k+1]*subst(Pnpoly(k, d, C), U, x));
    Phi = F2*subst(PN, x, us)/subst(subst(Qnpoly(qv,C),U,x), x, us);
    Xi = Dinv(Phi, NQ-2);
    Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
    Bv = peel2(Fs*Th, xs, NA, NQ);
    dn = 1; kmin = 0;
    for(n=1, NA, dn = lcm(dn,n); my(de, kj, tt); de = denominator(Bv[n+1]); kj = 0; tt = de; while(tt>1 && kj<12, kj++; tt = tt/gcd(tt,dn)); if(kj>kmin, kmin=kj));
    xi = Bv[NA+1]/Av[NA+1];
    dg = if(Bv[381]/Av[381]==Bv[301]/Av[301], 999, -log(abs(Bv[381]/Av[381]-Bv[301]/Av[301])*1.0)/log(10));
    print("   #", i, " alpha=", al);
    print("      magnetic to m=", magfail(Phi, NQ-2), " (0=all), c(1..6)=", vector(6,t,polcoeff(Phi,t)));
    print("      k = ", kmin, "   digits(n=300 vs 380) = ", dg);
    print("      xi = ", xi*1.0);
    print("      algdep(xi,2) = ", algdep(xi*1.0,2), "   lindep([xi,1,zeta(2),zeta(3)]) = ", lindep([xi*1.0,1,zeta(2),zeta(3)]));
  );
);
}
quit;
