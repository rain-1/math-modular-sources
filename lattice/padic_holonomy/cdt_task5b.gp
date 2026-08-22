/* Task 5b: scan the polynomial-degree bound D to locate the minimal-order ODE.
   ord = expected minimal order of the inhomogeneous linear ODE satisfied by H
         (= 2k+1 for X_0(p),  = 2 for the odd-character/Catalan cases).
   base     = {1, H, H', ..., H^{(ord-1)}}      : ord+1 functions, must stay FULL rank
   extended = {1, H, H', ..., H^{(ord)}}        : ord+2 functions, must DROP
*/
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/cdt_ab.gp");
BIGP = 2^61 - 1;
{
derivvec(cc, i) = vector(#cc - i, n1, my(n = n1-1); cc[n+i+1] * prod(t=1, i, n+t));
}
{
rk(flist, dd, ntr) =
  my(rows, mt);
  rows = List();
  for(i = 1, #flist, for(j = 0, dd,
      listput(rows, vector(ntr, t, if(t-1-j >= 0 && t-j <= #flist[i], flist[i][t-j], 0)))));
  matrank(Mod(matconcat(Col(Vec(rows))), BIGP));
}
{
scan(a, b, ord, et, tag, dmax, ntr) =
  my(cc, fl, base, ext, found);
  cc = vector(#a, i, a[i] + et*b[i]);
  fl = List(); listput(fl, vector(#cc, i, if(i==1, 1, 0)));
  for(i = 0, ord, listput(fl, derivvec(cc, i)));
  fl = Vec(fl);
  base = vector(ord+1, i, fl[i]);
  ext  = vector(ord+2, i, fl[i]);
  found = -1;
  print(tag, "  eta=", et, "  ord=", ord, "  base size=", ord+1, "  ext size=", ord+2);
  for(dd = 0, dmax,
    my(r1, r2, n1, n2);
    n1 = (ord+1)*(dd+1); n2 = (ord+2)*(dd+1);
    if(n2 > ntr - 5, print("      D=", dd, "  (truncation too small, stop)"); break);
    r1 = rk(base, dd, ntr); r2 = rk(ext, dd, ntr);
    if(r1 != n1 || r2 != n2,
      print("      D=", dd, "\t base rank=", r1, "/", n1, if(r1==n1," FULL"," *** DEFICIENT ***"),
                     "\t ext rank=", r2, "/", n2, if(r2==n2," FULL"," *** DEFICIENT ***"))
      ; if(found < 0 && r2 != n2, found = dd));
    if(found >= 0 && dd >= found+1, break));
  if(found < 0, print("      no relation found for D <= ", dmax),
                print("      => minimal-order-", ord, " inhomogeneous ODE with coefficient degrees <= ", found));
  print("");
}
NN = eval(getenv("NN")); if(NN == 0, NN = 420);
DM = eval(getenv("DM")); if(DM == 0, DM = 40);
NTR = NN - 10;
cases = [[2,1],[2,2],[2,3],[3,1],[3,2],[5,1],[5,2],[7,1]];
{
for(i=1, #cases,
  my(p, k, r);
  p = cases[i][1]; k = cases[i][2];
  r = run0p(p, k, NN);
  scan(r[1], r[2], 2*k+1, 5/7, Str("X_0(", p, ") k=", k), DM, NTR));
}
{
  my(r);
  r = run14(NN); scan(r[1], r[2], 2, 5/7, "X_1(4) chi_-4", DM, NTR);
  r = run09(NN); scan(r[1], r[2], 2, 5/7, "X_0(9) chi_-3", DM, NTR);
  r = run03(NN); scan(r[1], r[2], 2, 5/7, "X_0(3) chi_-3", DM, NTR);
}
quit;
