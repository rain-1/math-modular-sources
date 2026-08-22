default(parisizemax, 6000000000);
read("sc_rows.gp");
{
N = 800;
a = genseq(R4cf, 6, [1], N);
print("v_2(A_n) for n=1..40: ", vector(40,i,valuation(a[i+1],2)));
print("v_2(A_n) at n=100,200,400,800: ", [valuation(a[101],2),valuation(a[201],2),valuation(a[401],2),valuation(a[801],2)]);
for(j=1,5,
  x = compan(R4cf, 6, j, N);
  print("== R4 companion X^(", j, ") ==");
  for(pi=1,5,
    p = [2,3,5,7,13][pi];
    vv = vector(11, i, my(n = 700+i*10-10); my(z = x[n+1]/a[n+1] - x[n]/a[n]); if(z==0,"INF",valuation(z,p)));
    print("   p=", p, "  v_p(t_{n}-t_{n-1}) at n=700,710,...,800: ", vv);
  );
);
}
