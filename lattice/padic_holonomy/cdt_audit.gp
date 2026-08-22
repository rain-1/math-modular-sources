/* shared integrality audit helpers */
{
lcmto(n) = my(l=1); for(i=1,n, l=lcm(l,i)); l;
}
{
audit(aa, bb, e, nmax, tag) =
  my(badb, l, ok1, firstfail, failden, den, bad2);
  badb = [];
  for(i=1, nmax+1, if(denominator(bb[i]) != 1, badb = concat(badb, [i-1])));
  print(tag, "  b_n in Z for n=0..", nmax, " ? ", if(#badb==0, "YES", Str("NO at n=", badb)));
  ok1 = 1; bad2 = [];
  for(n=0, nmax,
    l = lcmto(max(n,1));
    if(denominator(l^e * aa[n+1]) != 1, ok1 = 0; bad2 = concat(bad2,[n])));
  print(tag, "  lcm(1..n)^", e, " * a_n in Z for n=0..", nmax, " ? ",
        if(ok1, "YES", Str("NO at n=", bad2)));
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
