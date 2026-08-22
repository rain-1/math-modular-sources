
\\ ---- T2 driver ----------------------------------------------------
\\ nm, cf, r, p, sigma, N, Nhalf : xi is taken as t_N = x_N/a_N.
t2row(nm, cf, r, p, sig, N) = {
  my(a = genseq(cf, r, [1], N));
  print("");
  print("################ ", nm, "   p=", p, "  sigma_p=", sig, "  N=", N, " ################");
  for(j = 1, r-1,
    my(x = compan(cf, r, j, N));
    my(t = vector(N+1));
    for(n = 1, N, t[n+1] = x[n+1]/a[n+1]);
    my(xi = t[N+1]);
    my(vxiprec = valuation(t[N+1]-t[N], p));
    print("");
    print("--- companion X^(", j, ") ---");
    print("  xi_p taken as t_N = x_N/a_N ; Cauchy precision v_p(t_N - t_{N-1}) = ", vxiprec);
    print("  (so xi_p is known to ", vxiprec, " p-adic digits; E is meaningful for n <~ N/2)");
    \\ (i) slope law
    my(fl = List());
    my(nlo = N\4, nhi = N\2);
    for(n = nlo, nhi, my(z = t[n+1]-xi); if(z, listput(fl, [n, valuation(z,p)])));
    my(v1 = fl[1][2], v2 = fl[#fl][2], n1 = fl[1][1], n2 = fl[#fl][1]);
    print("  (i) v_p(t_n - xi): n=", n1, " -> ", v1, " ;  n=", n2, " -> ", v2,
          " ; secant slope = ", (v2-v1)*1.0/(n2-n1));
    my(dmin = 10^9, dmax = -10^9, dsum = 0.);
    for(i=1,#fl, my(d = fl[i][2] - sig*fl[i][1]); if(d<dmin,dmin=d); if(d>dmax,dmax=d); dsum += d);
    print("      fluctuation  v_p(t_n - xi) - sigma*n  over n=", nlo, "..", nhi,
          " : min ", dmin, "  max ", dmax, "  mean ", dsum/#fl, "  range ", dmax-dmin);
    \\ also for x_n - xi a_n  (adds v_p(a_n))
    my(gmin = 10^9, gmax = -10^9);
    for(i=1,#fl, my(n=fl[i][1], d = fl[i][2] + if(a[n+1]==0,0,valuation(a[n+1],p)) - sig*n);
        if(d<gmin,gmin=d); if(d>gmax,gmax=d));
    print("      same for v_p(x_n - xi a_n) - sigma*n : min ", gmin, "  max ", gmax);
    \\ (ii) E(m,s)
    print("  (ii) E(m,s) = v_p(x_{m p^s} - xi a_{m p^s}) - sigma*m*p^s     [raw, includes v_p(a_n)]");
    print("       Er(m,s) = v_p(t_{m p^s} - xi)          - sigma*m*p^s     [ratio form]");
    my(Nsafe = N\2);
    for(m = 1, 6,
      my(rowE = List(), rowEr = List(), rowV = List());
      my(s = 0, ps = 1);
      while(m*ps <= Nsafe,
        my(n = m*ps, z = t[n+1]-xi);
        if(z == 0, listput(rowE,"INF"); listput(rowEr,"INF"); listput(rowV,"-"),
          my(vr = valuation(z,p), va = valuation(a[n+1],p));
          listput(rowEr, vr - sig*n); listput(rowE, vr+va - sig*n); listput(rowV, va));
        s++; ps *= p;
      );
      print("     m=", m, "  E : ", Vec(rowE));
      print("          Er: ", Vec(rowEr));
      print("          v_p(a_{m p^s}): ", Vec(rowV));
    );
    \\ (iii) Dwork-type congruence for the companion
    print("  (iii) v_p( t_{m p^{s+1}} - t_{m p^s} )  and  minus sigma*m*p^s");
    for(m = 1, 6,
      my(rowD = List(), rowDd = List());
      my(s = 0, ps = 1);
      while(m*ps*p <= N,
        my(z = t[m*ps*p+1] - t[m*ps+1]);
        if(z == 0, listput(rowD,"INF"); listput(rowDd,"INF"),
          my(v = valuation(z,p)); listput(rowD, v); listput(rowDd, v - sig*m*ps));
        s++; ps *= p;
      );
      print("     m=", m, "  v: ", Vec(rowD));
      print("          v - sigma*m*p^s: ", Vec(rowDd));
    );
  );
};
