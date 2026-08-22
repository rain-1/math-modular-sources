/* Sharpness detail: for each case list, up to n<=NMAX, the number of n with
   lcm(1..n)^{tau-1} a_n NOT in Z, the least such n, and the primes involved.   */
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
{
lcmto(n) = my(l=1); for(i=1,n,l=lcm(l,i)); l;
}
{
rep(a, tau, nmax, tag) =
  my(cnt, lst, prs, l, den);
  cnt = 0; lst = []; prs = [];
  for(n = 1, nmax,
    l = lcmto(n);
    den = denominator(l^(tau-1)*a[n+1]);
    if(den != 1, cnt++; if(#lst < 6, lst = concat(lst, [n]));
       prs = setunion(prs, Set(factor(den)[,1]~))));
  print(tag, "  tau=", tau, ":  #{n<=", nmax, " : lcm^", tau-1, " a_n not in Z} = ", cnt,
        ";  first such n = ", lst, ";  primes occurring = ", prs);
}
NMAX = 60;
cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
{
for(i=1,#cases, my(p,k,r); p=cases[i][1]; k=cases[i][2];
  r = run0p(p,k,NMAX); rep(r[1], 2*k+1, NMAX, Str("X_0(",p,") k=",k)));
}
rep(run14(NMAX)[1], 2, NMAX, "X_1(4) chi_-4");
rep(run09(NMAX)[1], 2, NMAX, "X_0(9) chi_-3");
rep(run03(NMAX)[1], 2, NMAX, "X_0(3) chi_-3");
quit;
