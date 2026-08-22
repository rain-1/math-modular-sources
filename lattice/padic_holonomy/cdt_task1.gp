read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");

{
lcmto(n) = my(l=1); for(i=1,n, l=lcm(l,i)); l;
}

/* integrality audit: aa,bb vectors with aa[i]=a_{i-1}; expected exponent e */
{
audit(aa, bb, e, nmax, tag) =
  my(badb, l, ok1, firstfail, failden, den);
  badb = [];
  for(i=1, nmax+1, if(denominator(bb[i]) != 1, badb = concat(badb, [i-1])));
  print(tag, "  b_n in Z for n=0..", nmax, " ? ", if(#badb==0, "YES", Str("NO at n=", badb)));
  ok1 = 1;
  for(n=0, nmax,
    l = lcmto(max(n,1));
    if(denominator(l^e * aa[n+1]) != 1, ok1 = 0; print(tag, "  *** lcm^",e," FAILS at n=",n)));
  print(tag, "  lcm(1..n)^", e, " * a_n in Z for n=0..", nmax, " ? ", if(ok1, "YES", "NO"));
  firstfail = -1; failden = 0;
  for(n=0, nmax,
    l = lcmto(max(n,1));
    den = denominator(l^(e-1) * aa[n+1]);
    if(den != 1 && firstfail < 0, firstfail = n; failden = den));
  if(firstfail < 0,
    print(tag, "  SHARPNESS: lcm^", e-1, " also integral for all n<=", nmax, " -> NOT sharp in range"),
    print(tag, "  SHARPNESS: least n with lcm^", e-1, "*a_n not in Z: n=", firstfail,
          "  residual denominator=", failden, " = ", factor(failden)));
  print("");
}

{
docase(p, k, nmax) =
  my(r, tag);
  r = run0p(p, k, nmax);
  tag = Str("  [X0(", p, ") k=", k, "]");
  print("=== X_0(", p, "), k=", k, "  weight 2k=", 2*k, ", expected tau = 2k+1 = ", 2*k+1, " ===");
  print("  b_0..b_6 = ", vector(min(7,nmax+1), j, r[2][j]));
  print("  a_0..a_6 = ", vector(min(7,nmax+1), j, r[1][j]));
  audit(r[1], r[2], 2*k+1, nmax, tag);
  write(Str("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/data_X0_", p, "_k", k, ".txt"), r);
  r;
}

NMAX = eval(getenv("NMAX"));
if(NMAX == 0, NMAX = 40);
cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
for(i=1, #cases, docase(cases[i][1], cases[i][2], NMAX));
quit;
