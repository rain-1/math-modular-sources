default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 200);
NQ = 402; NA = 380; M2 = 150; M3 = 380;
HH = HOSTS[3]; C=HH[2]; B=HH[3]; dv=HH[4]; rv=HH[5];
us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
xs = us/(1+B*us+C*us^2);
Av = peel2(Fs, xs, NA, NQ);
print("### magnetic sources with polar support in the cusp pair {1/2, 1/3} (W = -17, x = infinity)");
print("### -- the ONLY polar locus on Apery's host that contributes no |x|^(-n) tail.");
print("");
print("(a) the W_6-ANTI-INVARIANT subspace: rho = v * P(w)/(72w+17)^d, deg P <= d-2, dim d-1.");
{
for(d=2, 12,
  my(Qw, qv, Xs, T, S);
  Qw = (72*w+17)^d;
  qv = vector(poldegree(Qw)+1, t, polcoeff(Qw, t-1));
  Xs = xibasis(C, us, F2, qv, M3);
  T = xitomat(Xs, M3);
  S = magsub(T, M2, M3);
  print("   d=", d, "  dim=", #Xs, "  # magnetic basis vectors surviving to m=", M3, " : ", #S, if(#S>0, Str("   ", S), ""));
);
}
print("");
print("(b) the FULL space (both Fricke parities): rho = u P(u)/(72u^2+17u+1)^d, deg P <= 2d-1, dim 2d.");
{
for(d=1, 6,
  my(Qp, qv, Xs, T, S);
  Qp = (72*x^2+17*x+1)^d;
  qv = vector(poldegree(Qp)+1, t, polcoeff(Qp, t-1));
  Xs = xibasisW(us, F2, qv, M3);
  T = xitomat(Xs, M3);
  S = magsub(T, M2, M3);
  print("   d=", d, "  dim=", #Xs, "  # magnetic : ", #S);
  for(i=1, #S,
    my(al, NUM, RI, sgn, Phi, Xi, Th, Bv, dn, kmin, xi);
    al = S[i];
    NUM = sum(k=0, #Xs-1, al[k+1]*x^(k+1));
    RI = simplify(subst(NUM,x,1/(C*x))/subst(Qp,x,1/(C*x)));
    sgn = if(RI + NUM/Qp == 0, "anti-invariant (W_6 = -1)", if(RI - NUM/Qp == 0, "invariant (W_6 = +1)", "MIXED"));
    Phi = F2*subst(NUM,x,us)/subst(Qp,x,us);
    Xi = Dinv(Phi, NQ-2);
    Th = sum(m=1, NQ-2, polcoeff(Xi,m)/m^2*q^m) + O(q^(NQ-1));
    Bv = peel2(Fs*Th, xs, NA, NQ);
    dn = 1; kmin = 0;
    for(n=1, NA, dn = lcm(dn,n); my(de, kj, tt); de = denominator(Bv[n+1]); kj = 0; tt = de; while(tt>1 && kj<12, kj++; tt = tt/gcd(tt,dn)); if(kj>kmin, kmin=kj));
    xi = Bv[NA+1]/Av[NA+1];
    print("      #", i, "  alpha=", al, "   ", sgn);
    print("          k = ", kmin, "   B_n/A_n at n=380 = ", xi*1.0);
    print("          digits of agreement n=300 vs n=380 : ", if(Bv[381]/Av[381]==Bv[301]/Av[301], "exact", -log(abs(Bv[381]/Av[381]-Bv[301]/Av[301])*1.0)/log(10)));
  );
);
}
quit;
