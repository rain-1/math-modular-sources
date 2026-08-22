default(parisizemax, 6000000000);
read("sc_rows.gp");
\\ Is E(m,s) eventually constant in s (for fixed m)?  And is it bounded?
t2e(nm, cf, r, p, sig, N, MMAX, s0) = {
  my(av = genseq(cf, r, [1], N));
  print("");
  print("#### ", nm, "  p=", p, " sigma=", sig, " N=", N, "  s>=", s0);
  for(j = 1, r-1,
    my(xv = compan(cf, r, j, N));
    my(t = vector(N+1)); for(n=1,N, t[n+1] = xv[n+1]/av[n+1]);
    my(xi = t[N+1]);
    my(cells = 0, viol = 0, fb = 0, emin = 10^9, emax = -10^9, ms = 0);
    for(m = 1, MMAX,
      my(es = List(), s = s0, ps = p^s0);
      while(m*ps <= N\2,
        my(n = m*ps, z = t[n+1]-xi);
        if(z != 0, listput(es, [s, valuation(z,p)+valuation(av[n+1],p) - sig*n]));
        s++; ps *= p;
      );
      es = Vec(es);
      if(#es >= 2, ms++;
        for(i = 1, #es, cells++;
          if(es[i][2] < emin, emin = es[i][2]); if(es[i][2] > emax, emax = es[i][2]);
          if(es[i][2] != es[1][2], viol++; if(fb==0, fb = [m, es[i][1], es[i][2], es[1][2]])));
      );
    );
    print("   X^(",j,"): ", ms, " values of m with >=2 usable s; ", cells, " cells; ",
          viol, " violations of 's-constancy'", if(viol, Str("  first (m,s,E,E_first)=", fb), ""),
          "   E range [", emin, ", ", emax, "]");
  );
};
t2e("R2 Zagier C", R2cf, 2, 3, 2, 4000, 300, 2);
t2e("R3 AZ eta", R3cf, 2, 5, 3, 4000, 300, 2);
t2e("R5 AESZ 207", R5cf, 4, 2, 12, 4000, 300, 2);
t2e("R5 AESZ 207", R5cf, 4, 2, 12, 4000, 300, 6);
