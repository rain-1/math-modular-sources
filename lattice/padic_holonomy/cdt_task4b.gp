/* Refined archimedean growth: r_n := (|b_{2n}|/|b_n|)^{1/n} -> 1/R_arch.
   (Immune to the phase oscillation caused by complex-conjugate singularities.)  */
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
default(realprecision, 12);
{
growth(r, tag, expect) =
  my(a, b, nn);
  a = r[1]; b = r[2];
  nn = #b - 1;
  print(tag, "   predicted 1/R_arch = ", expect);
  for(i = 1, 5,
    my(n = (nn\2)*i\5);
    if(n >= 10 && b[n+1] != 0 && b[2*n+1] != 0,
      print("      n=", n, "\t (|b_2n|/|b_n|)^(1/n) = ", exp(log(abs(b[2*n+1]*1.0/b[n+1]))/n),
                          "\t (|a_2n|/|a_n|)^(1/n) = ", if(a[n+1]!=0, exp(log(abs(a[2*n+1]*1.0/a[n+1]))/n), 0))));
}
NB = eval(getenv("NBIG")); if(NB == 0, NB = 400);
cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
{
for(i=1, #cases,
  my(p, k);
  p = cases[i][1]; k = cases[i][2];
  growth(run0p(p, k, NB), Str("X_0(", p, ") k=", k), exp(log(p*1.0)*6/(p-1))));
}
growth(run14(NB), "X_1(4) chi_-4 ", 16.0);
growth(run09(NB), "X_0(9) chi_-3 ", 27.0^0.5);
growth(run03(NB), "X_0(3) chi_-3 ", 27.0);
quit;
