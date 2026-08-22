default(parisizemax, 6000000000);
read("sc_rows.gp");
{
N = 6000; p = 2; sig = 12;
av = genseq(R5cf, 4, [1], N);
print("#### R5 AESZ 207, p=2, sigma=12, N=", N, "  (xi = t_N, precision v_2 ~ ", valuation(compan(R5cf,4,1,N)[N+1]/av[N+1] - compan(R5cf,4,1,N)[N]/av[N], 2), ")");
for(j = 1, 3,
  xv = compan(R5cf, 4, j, N);
  xi = xv[N+1]/av[N+1];
  cnt = vector(15); mn = vector(15,i,10^9); mx = vector(15,i,-10^9);
  for(n = 2, N\2,
    z = xv[n+1]/av[n+1] - xi;
    if(z != 0,
      e = valuation(z,2) + valuation(av[n+1],2) - sig*n;
      vn = min(valuation(n,2),13)+1;
      cnt[vn]++; if(e<mn[vn],mn[vn]=e); if(e>mx[vn],mx[vn]=e);
    );
  );
  print("  X^(",j,"):  n = 2..", N\2);
  for(i=1,14, if(cnt[i], print("     v_2(n)=",i-1," : ",cnt[i]," values, E in [",mn[i],",",mx[i],"]")));
);
}
