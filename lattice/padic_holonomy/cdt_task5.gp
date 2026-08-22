/* Task 5: Q(x)-linear independence of the CDT inventory
      {1, H, H', ..., H^{(2k)}}   (m = 2k+2 functions),  ' = d/dx,
   with eta specialised to a RATIONAL number (the "zeta_p(2k+1) in Q" hypothesis).
   Rank of { x^j * f_i : i=1..m, j=0..D } inside Q[[x]]/x^NTR, computed mod
   P = 2^61-1.  Full rank = m*(D+1) means no Q(x)-relation with numerator
   degrees <= D.  Also reports the extended inventory {1,H,...,H^{(2k+1)}}
   (m+1 functions), where a rank DROP certifies a linear ODE of order 2k+1.
*/
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
BIGP = 2^61 - 1;

/* coefficient vector of d^i/dx^i of the series with coeffs cc (cc[n+1] = coeff of x^n) */
{
derivvec(cc, i) =
  my(nn = #cc - 1, out);
  out = vector(nn+1-i, n1, my(n = n1-1); cc[n+i+1] * prod(t=1, i, n+t));
  out;
}
/* rank of {x^j f_i} truncated to NTR coefficients */
{
rk(flist, dd, ntr) =
  my(rows, mt);
  rows = List();
  for(i = 1, #flist,
    for(j = 0, dd,
      my(v = vector(ntr, t, if(t-1-j >= 0 && t-j <= #flist[i], flist[i][t-j], 0)));
      listput(rows, v)));
  mt = Mod(matconcat(Col(Vec(rows))), BIGP);
  matrank(mt);
}
{
ranktable(a, b, k, et, tag, dmax, ntr) =
  my(cc, flist, m2, r, need);
  cc = vector(#a, i, a[i] + et*b[i]);
  m2 = 2*k;
  flist = List(); listput(flist, vector(#cc, i, if(i==1, 1, 0)));
  for(i = 0, m2+1, listput(flist, derivvec(cc, i)));
  flist = Vec(flist);   /* [1, H, H', ..., H^{(2k+1)}] : 2k+3 entries */
  print(tag, "   eta = ", et, "    m = 2k+2 = ", 2*k+2, ",  m+1 = ", 2*k+3);
  for(dd = 0, dmax,
    my(f1, f2, r1, r2);
    f1 = vector(2*k+2, i, flist[i]);      /* 1, H, ..., H^{(2k)}   */
    f2 = vector(2*k+3, i, flist[i]);      /* 1, H, ..., H^{(2k+1)} */
    r1 = rk(f1, dd, ntr);
    r2 = rk(f2, dd, ntr);
    print("      D=", dd,
          "\t m=", 2*k+2, ": rank=", r1, "/", (2*k+2)*(dd+1), if(r1==(2*k+2)*(dd+1), "  FULL", "  DEFICIENT"),
          "\t m+1=", 2*k+3, ": rank=", r2, "/", (2*k+3)*(dd+1), if(r2==(2*k+3)*(dd+1), "  FULL", "  DEFICIENT")));
  print("");
}
NN  = eval(getenv("NN"));  if(NN  == 0, NN  = 260);
DMAX = eval(getenv("DMAX")); if(DMAX == 0, DMAX = 5);
NTR = NN - 2*4 - DMAX - 5;
cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
{
for(i=1, #cases,
  my(p, k, r);
  p = cases[i][1]; k = cases[i][2];
  r = run0p(p, k, NN);
  ranktable(r[1], r[2], k, 0,   Str("X_0(", p, ") k=", k), DMAX, NTR);
  ranktable(r[1], r[2], k, 5/7, Str("X_0(", p, ") k=", k), DMAX, NTR));
}
{
  my(r);
  r = run14(NN);
  ranktable(r[1], r[2], 0, 0,   "X_1(4) chi_-4", DMAX, NTR);
  ranktable(r[1], r[2], 0, 5/7, "X_1(4) chi_-4", DMAX, NTR);
  r = run09(NN);
  ranktable(r[1], r[2], 0, 0,   "X_0(9) chi_-3", DMAX, NTR);
  ranktable(r[1], r[2], 0, 5/7, "X_0(9) chi_-3", DMAX, NTR);
}
quit;
