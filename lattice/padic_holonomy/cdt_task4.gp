/* Task 4: p-adic slope s_n = v_p(a_n + eta b_n)/n  and archimedean growth |b_n|^{1/n}.
   eta is taken exactly as eta_N = -a_N/b_N with N = NBIG; then
       a_n + eta_N b_n = (a_n b_N - a_N b_n)/b_N
   so v_p(a_n + eta_N b_n) = v_p(a_n b_N - a_N b_n) - v_p(b_N), computed EXACTLY.
   For n << N this equals v_p(a_n + eta_true b_n).
*/
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");

{
slopes(r, p, nlist, nbig, tag, expect) =
  my(a, b, vN, w, vv);
  a = r[1]; b = r[2];
  vN = valuation(b[nbig+1], p);
  print(tag, "   [eta = -a_", nbig, "/b_", nbig, "]   expected slope = ", expect);
  for(i=1, #nlist,
    my(n = nlist[i], vc, dd);
    dd = a[n+1]*b[nbig+1] - a[nbig+1]*b[n+1];
    if(dd == 0, print("      n=", n, "  a_n b_N - a_N b_n = 0 !"),
      vc = valuation(dd, p) - vN;
      print("      n=", n, "\t v_p(a_n+eta b_n) = ", vc, "\t s_n = ", vc/n*1.0)));
}
{
arch(r, nlist, tag, expect) =
  my(a, b);
  a = r[1]; b = r[2];
  print(tag, "   expected lim|b_n|^{1/n} = ", expect);
  for(i=1, #nlist,
    my(n = nlist[i]);
    print("      n=", n,
          "\t |b_n|^{1/n} = ", if(b[n+1]!=0, exp(log(abs(b[n+1]*1.0))/n), 0),
          "\t |b_n/b_{n-1}| = ", if(b[n]!=0, abs(b[n+1]*1.0/b[n]), 0),
          "\t |a_n|^{1/n} = ", if(a[n+1]!=0, exp(log(abs(a[n+1]*1.0))/n), 0)));
}

NB   = eval(getenv("NBIG"));  if(NB   == 0, NB = 160);
NL = [10, 20, 40, 60, 80, 100, 120];
NLA = [10, 20, 40, 60, 80, 100, 120, 140, 160];
default(realprecision, 10);

cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
{
for(i=1, #cases,
  my(p, k, r, tag);
  p = cases[i][1]; k = cases[i][2];
  r = run0p(p, k, NB);
  tag = Str("X_0(", p, ") k=", k);
  print("--- ", tag, " ---");
  slopes(r, p, NL, NB, "  p-adic:", 12/(p-1)*1.0);
  arch(r, NLA, "  arch:  ", exp(log(p*1.0)*6/(p-1)));
  print(""));
}
{
  my(r);
  r = run14(NB);
  print("--- X_1(4), chi_{-4} (2-adic Catalan) ---");
  slopes(r, 2, NL, NB, "  p-adic:", 8.0);
  arch(r, NLA, "  arch:  ", 16.0);
  print("");
  r = run09(NB);
  print("--- X_0(9), chi_{-3} (3-adic) ---");
  slopes(r, 3, NL, NB, "  p-adic:", 3.0);
  arch(r, NLA, "  arch:  ", sqrt(27.0));
  print("");
  r = run03(NB);
  print("--- X_0(3), chi_{-3} (level-3 comparison) ---");
  slopes(r, 3, NL, NB, "  p-adic:", 6.0);
  arch(r, NLA, "  arch:  ", 27.0);
  print("");
}
quit;
