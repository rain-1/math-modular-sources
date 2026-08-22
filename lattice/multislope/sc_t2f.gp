default(parisizemax, 6000000000);
read("sc_rows.gp");
\\ E(n) = v_p(x_n - xi a_n) - sigma*n  as a function of v_p(n), for ALL n (not just tower n)
t2f(nm, cf, r, p, sig, N) = {
  my(av = genseq(cf, r, [1], N));
  print("");
  print("#### ", nm, " p=", p, " sigma=", sig, "  E(n) vs v_p(n),  n = 100..", N\2);
  for(j = 1, r-1,
    my(xv = compan(cf, r, j, N));
    my(xi = xv[N+1]/av[N+1]);
    my(tab = matrix(14, 1));  \\ per v_p(n): [count, min, max]
    my(cnt = vector(14), mn = vector(14, i, 10^9), mx = vector(14, i, -10^9));
    for(n = 100, N\2,
      my(z = xv[n+1]/av[n+1] - xi);
      if(z != 0,
        my(e = valuation(z,p) + valuation(av[n+1],p) - sig*n);
        my(vn = min(valuation(n,p), 12) + 1;);
        cnt[vn]++; if(e<mn[vn], mn[vn]=e); if(e>mx[vn], mx[vn]=e);
      );
    );
    print("   X^(",j,"):");
    for(i = 1, 13, if(cnt[i], print("      v_p(n)=", i-1, " : ", cnt[i], " values of n, E in [", mn[i], ",", mx[i], "]")));
  );
};
t2f("R5 AESZ 207", R5cf, 4, 2, 12, 4000);
t2f("R2 Zagier C", R2cf, 2, 3, 2, 4000);
t2f("R3 AZ eta", R3cf, 2, 5, 3, 4000);
