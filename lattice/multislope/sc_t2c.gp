\\ ---- T2c: test the exact stabilisation law E(m,s) = c_j -------------
\\ E(m,s) = v_p(x_{m p^s} - xi a_{m p^s}) - sigma*m*p^s, xi = t_N.
t2c(nm, cf, r, p, sig, N, MMAX, s0) = {
  my(av = genseq(cf, r, [1], N));
  print("");
  print("######## T2c ", nm, " p=", p, " sigma=", sig, " N=", N, " (E meaningful for n<=", N\2, ") ########");
  for(j = 1, r-1,
    my(xv = compan(cf, r, j, N));
    my(t = vector(N+1)); for(n=1,N, t[n+1] = xv[n+1]/av[n+1]);
    my(xi = t[N+1]);
    my(cnt = 0, fail = 0, cj = "unset", firstbad = 0, vals = List());
    for(m = 1, MMAX,
      my(s = 0, ps = 1);
      while(m*ps <= N\2,
        if(s >= s0,
          my(n = m*ps, z = t[n+1]-xi);
          if(z != 0,
            my(e = valuation(z,p) + valuation(av[n+1],p) - sig*n);
            cnt++;
            if(cj == "unset", cj = e);
            if(e != cj, fail++; if(firstbad==0, firstbad=[m,s,e]));
            listput(vals, e);
          );
        );
        s++; ps *= p;
      );
    );
    print("  X^(", j, "):  c_j = ", cj, "   cells tested (m<=", MMAX, ", s>=", s0, "): ", cnt,
          "   deviations: ", fail, if(fail, Str("  first ", firstbad), ""));
    if(fail, print("      observed E values: ", Vec(vals)));
  );
};
