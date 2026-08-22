default(parisizemax, 6000000000);
read("sc_rows.gp");
{
N = 8192;
av = genseq(R4cf, 6, [1], N);
print("#### R4 Sym^3 Zagier E: good-prime tower, u(a,s) = eps^{-1} p^w t_{a p^{s+1}}/t_{a p^s}");
specs = [[2, 2, -1], [3, 4, -1], [7, 4, -1], [5, 4, 1], [13, 4, 1]];
for(j = 1, 5,
  xv = compan(R4cf, 6, j, N);
  print("  X^(", j, "):");
  for(si = 1, #specs,
    p = specs[si][1]; w = specs[si][2]; ep = specs[si][3];
    for(s = 0, 12,
      ps = p^s;
      if(ps*p > N, break);
      vals = List(); amax = N\(ps*p);
      for(aa = 1, amax,
        n0 = aa*ps; n1 = aa*ps*p;
        if(xv[n0+1]==0 || xv[n1+1]==0, next);
        u = ep * p^w * (xv[n1+1]/av[n1+1]) / (xv[n0+1]/av[n0+1]);
        listput(vals, if(u==1, 999, valuation(u-1, p)));
      );
      vals = Vec(vals);
      if(#vals >= 1 && s <= 6,
        print("     p=",p," eps=",ep," w=",w," s=", s, " (", #vals, " cells a=1..", amax, "): v_p(u-1) min=", vecmin(vals), " max=", vecmax(vals)));
    );
  );
);
}
