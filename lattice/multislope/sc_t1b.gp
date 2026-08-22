default(parisizemax, 6000000000);
read("sc_rows.gp");
{
N = 3000;
prs = [2,3,5,7,13];
rows = [["R1 Apery zeta(3)", R1cf, 2],
        ["R2 Zagier C", R2cf, 2],
        ["R3 AZ eta", R3cf, 2],
        ["R4 Sym^3 Zagier E", R4cf, 6],
        ["R5 AESZ 207", R5cf, 4]];
for(ri = 1, #rows,
  my(nm = rows[ri][1], cf = rows[ri][2], r = rows[ri][3]);
  my(a = genseq(cf, r, [1], N));
  print("### ", nm, "  --- excess exponent  v_p(a_{m p^s} - a_{m p^{s-1}}) / s   (min over m, per s)");
  for(pi = 1, #prs,
    my(p = prs[pi]);
    my(out = List());
    my(s = 1, ps = p);
    while(ps <= N,
      my(mn = 10^9, mnm = 0);
      for(m = 1, N\ps,
        my(d = a[m*ps+1] - a[m*(ps\p)+1]);
        my(v = if(d==0, 10^6, valuation(d,p)));
        if(v < mn, mn = v; mnm = m);
      );
      listput(out, [s, mn, mn*1.0/s, mnm]);
      s++; ps *= p;
    );
    print("   p=", p, ": [s, min_m v_p(diff), ratio v/s, argmin m] = ", Vec(out));
  );
);
}
