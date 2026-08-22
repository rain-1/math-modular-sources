default(parisizemax, 6000000000);
read("sc_rows.gp");
{
N = 8192; p = 13; w = 2;
av = genseq(R5cf, 4, [1], N);
print("#### R5 AESZ 207, good-prime tower at p=13:  u(a,s) = p^",w," * t_{a p^{s+1}}/t_{a p^s},  t_n = x_n/a_n");
print("     tabulating v_p(u - 1) : min / max over a coprime-free range, per level s");
for(j = 1, 3,
  xv = compan(R5cf, 4, j, N);
  print("  X^(", j, "):");
  for(s = 0, 3,
    ps = p^s;
    if(ps*p > N, break);
    vals = List(); amax = N\(ps*p);
    for(aa = 1, amax,
      n0 = aa*ps; n1 = aa*ps*p;
      if(xv[n0+1]==0 || xv[n1+1]==0, next);
      u = p^w * (xv[n1+1]/av[n1+1]) / (xv[n0+1]/av[n0+1]);
      listput(vals, if(u==1, 999, valuation(u-1, p)));
    );
    vals = Vec(vals);
    if(#vals, print("     s=", s, " (a=1..", amax, ", ", #vals, " cells): v_p(u-1) min=", vecmin(vals),
                    " max=", vecmax(vals), "  #(v>=1)=", sum(i=1,#vals,vals[i]>=1), " #(v>=2)=", sum(i=1,#vals,vals[i]>=2)));
  );
);
}
