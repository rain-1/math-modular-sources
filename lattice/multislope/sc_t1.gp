default(parisizemax, 6000000000);
read("sc_rows.gp");
{
NT1 = 3000;
prs = [2,3,5,7,13];
rows = [["R1 Apery zeta(3)", R1cf, 2],
        ["R2 Zagier C (10,3,9)", R2cf, 2],
        ["R3 AZ eta (11,5,125)", R3cf, 2],
        ["R4 Sym^3 Zagier E", R4cf, 6],
        ["R5 AESZ 207", R5cf, 4]];
for(ri = 1, #rows,
  my(nm = rows[ri][1], cf = rows[ri][2], r = rows[ri][3]);
  gettime();
  my(a = genseq(cf, r, [1], NT1));
  my(gt = gettime());
  print("### ", nm, "   (gen ", gt, " ms)");
  my(tot = 0, totf = 0);
  for(pi = 1, #prs,
    my(p = prs[pi], nt = 0, nf = 0, firstf = 0);
    my(s = 1, ps = p);
    while(ps <= NT1,
      for(m = 1, NT1\ps,
        nt++;
        if((a[m*ps+1] - a[m*(ps\p)+1]) % p^s != 0,
          nf++; if(firstf == 0, firstf = [m, s, valuation(a[m*ps+1]-a[m*(ps\p)+1], p)]);
        );
      );
      s++; ps *= p;
    );
    tot += nt; totf += nf;
    print("   p=", p, ": ", nt, " tests, ", nf, " failures",
          if(nf, concat("  first (m,s,v_p(diff)) = ", Str(firstf)), ""));
  );
  print("   TOTAL ", tot, " tests, ", totf, " failures");
);
}
