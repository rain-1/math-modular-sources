
vp(z, p) = if(z == 0, "INF", valuation(z, p));

rowdiag(nm, cf, r, N, primes) = {
  print("==================== ", nm, "  (order ", r, ", N=", N, ") ====================");
  my(a = genseq(cf, r, [1], N));
  print("  A_0..A_5 = ", vector(6,i,a[i]));
  print("  denexp(A) = ", denexp(a, N));
  for(j = 1, r-1,
    my(x = compan(cf, r, j, N));
    print("  companion X^(", j, "):  x_0..x_", min(6,N), " = ", vector(min(7,N+1),i,x[i]));
    print("    denominator exponent k = ", denexp(x, N));
    for(pi = 1, #primes,
      my(p = primes[pi]);
      \\ slope estimate: v_p(x_{n+1}/a_{n+1} - x_n/a_n)/n for n near N
      my(t1 = x[N+1]/a[N+1] - x[N]/a[N]);
      my(t2 = x[N]/a[N] - x[N-1]/a[N-1]);
      my(v1 = if(t1==0, "INF", valuation(t1,p)), v2 = if(t2==0,"INF",valuation(t2,p)));
      my(qn = round(N/2));
      my(t3 = x[qn+1]/a[qn+1] - x[qn]/a[qn]);
      my(v3 = if(t3==0,"INF",valuation(t3,p)));
      print("    p=", p, ": v_p(dt) at n=", N-1, ": ", v2, "  n=", N, ": ", v1,
            "   -> sigma ~ ", if(type(v1)=="t_INT" && type(v3)=="t_INT", (v1-v3)*1.0/(N-qn), "?"));
    );
  );
};
