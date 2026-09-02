default(parisizemax,12000000000);
read("lib.gp");
default(realprecision, 120);
NQ = 405; MA=60; MB=150; MC=400;
uvalh(tt, dvv, rvv) = prod(j=1, #dvv, eta(dvv[j]*tt, 1)^rvv[j]);
prim(P) = my(d,g); d=1; for(j=0,poldegree(P), d=lcm(d,denominator(polcoeff(P,j)))); P=P*d; g=content(P); P/g;
print("### Task 4: the level-5 host and the two other level-6 hosts (and the level-10/12/18 hosts).");
print("### For each: (i) singular values of W = C u + 1/u = 1/x - B  at CM points of X_0(N),");
print("###           (ii) systematic search for magnetic sources beyond the canonical one.");
IDX = [1,2,4,5,9,10,11,12];
{
for(t=1, #IDX,
  my(HH, N, C, Bh, dv, rv, us, Fs, F2, DL, nh);
  HH = HOSTS[IDX[t]]; N=HH[1]; C=HH[2]; Bh=HH[3]; dv=HH[4]; rv=HH[5];
  print("");
  print("=================== ", HH[6], "   N=",N," C=",C," B=",Bh);
  \\ (i) CM singular values of W
  DL = [-3,-4,-7,-8,-11,-12,-15,-16,-19,-20,-23,-24,-27,-28,-32,-35,-36,-40,-43,-48,-51,-52,-64,-67,-72,-84,-96,-100,-120];
  print("  CM singular values of W = ", C, "u + 1/u   (x = 1/(W+", Bh, ") is the host coordinate)");
  print("  disc | minpoly(W) monic? | N(W+B) | max_i |W_i+B| = max 1/|x_i|");
  for(s=1, #DL,
    my(D, vals, ws, PW, PY, nrm, mx, rts);
    D = DL[s];
    vals = List();
    for(k=1, 30,
      my(a); a = N*k;
      for(b = -a+1, a,
        if((b^2-D) % (4*a) != 0, next);
        my(c, tau, u0, wv, fnd);
        c = (b^2-D)/(4*a);
        if(gcd(gcd(a,b),c) != 1, next);
        tau = (-b + sqrt(-D)*I)/(2*a);
        u0 = uvalh(tau, dv, rv);
        wv = C*u0 + 1/u0;
        fnd = 0; for(j=1,#vals, if(abs(wv-vals[j]) < 1e-60, fnd=1; break));
        if(fnd==0, listput(vals, wv));
      );
    );
    if(#vals == 0, next);
    if(#vals > 12, print("   D=", D, "  degree ", #vals, " - skipped"); next);
    PW = prim(bestappr(real(prod(j=1,#vals,(X - vals[j]))), 10^45));
    PY = subst(PW, X, Y - Bh);
    nrm = polcoeff(PY,0)*(-1)^poldegree(PY)/pollead(PY);
    rts = polroots(PY); mx = 0; for(j=1,#rts, if(abs(rts[j])>mx, mx=abs(rts[j])));
    print("   D=", D, " | ", if(pollead(PW)==1, "yes", concat("NO, lead=",pollead(PW))), " | N(W+B) = ", nrm, " | max = ", mx);
  );
  \\ (ii) scan
  us = useries(dv,rv,NQ); Fs = Fseries(dv,rv,NQ); F2 = Fs^2;
  print("  systematic scan: Q = (W+b)^e, widest space, |b| <= 300, e = 1,2,3");
  nh = 0;
  for(b=-300, 300,
    for(e=1, 3,
      my(hu, Qp, qv, Xs, T, res);
      hu = C*x^2 + b*x + 1;
      Qp = hu^e;
      qv = vector(poldegree(Qp)+1, j, polcoeff(Qp, j-1));
      Xs = xibasisW(us, F2, qv, MC);
      if(type(Xs)=="t_INT", next);
      T = xitomat(Xs, MC);
      res = scanT(T, MA, MB, MC);
      if(res[1]==1, nh=nh+1; print("    *** MAGNETIC *** Q = (", C, "u^2 + ", b, "u + 1)^", e, "  dim=", #Xs, "  alpha=", res[2]));
    );
  );
  print("  -> ", nh, " magnetic hits in the (W+b)^e family");
);
}
quit;
