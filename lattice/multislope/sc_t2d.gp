default(parisizemax, 6000000000);
read("sc_rows.gp");
\\ exact Casoratian law for the rank-2 rows:
\\   W_n = a_n x_{n+1} - a_{n+1} x_n = c^n / (n+1)^k     (c = 9 resp 125, k = 2 resp 3)
\\   => v_p(t_n - xi) = sigma*n - k*v_p(n+1) - v_p(a_n) - v_p(a_{n+1})   (tail dominated by first term)
\\   => E(m,s) := v_p(x_n - xi a_n) - sigma*n = -k*v_p(n+1) - v_p(a_{n+1}),  n = m p^s
chk2(nm, cf, p, sig, kk, cc, N, MMAX) = {
  my(av = genseq(cf, 2, [1], N), xv = compan(cf, 2, 1, N));
  print("");
  print("#### ", nm, " p=", p, " N=", N);
  \\ (0) verify the Casoratian formula exactly
  my(bad = 0);
  for(n = 0, N-1, if(av[n+1]*xv[n+2] - av[n+2]*xv[n+1] != cc^n/(n+1)^kk, bad++));
  print("   Casoratian  a_n x_{n+1} - a_{n+1} x_n = ", cc, "^n/(n+1)^", kk, " : ", N, " checks, ", bad, " failures");
  my(t = vector(N+1)); for(n=1,N, t[n+1] = xv[n+1]/av[n+1]);
  my(xi = t[N+1]);
  my(cnt = 0, fl = 0, fb = 0);
  for(m = 1, MMAX,
    my(s = 1, ps = p);
    while(m*ps <= N\2,
      my(n = m*ps, z = t[n+1]-xi);
      if(z != 0,
        my(e = valuation(z,p) + valuation(av[n+1],p) - sig*n);
        my(pred = -kk*valuation(n+1,p) - valuation(av[n+2],p));
        cnt++;
        if(e != pred, fl++; if(fb==0, fb=[m,s,e,pred]));
      );
      s++; ps *= p;
    );
  );
  print("   law  E(m,s) = -", kk, "*v_p(mp^s+1) - v_p(a_{mp^s+1})  (m<=", MMAX, ", s>=1, mp^s<=", N\2, "): ",
        cnt, " cells, ", fl, " failures", if(fl, Str("  first ", fb), ""));
};
chk2("R2 Zagier C", R2cf, 3, 2, 2, 9, 4000, 300);
chk2("R3 AZ eta", R3cf, 5, 3, 3, 125, 4000, 300);
