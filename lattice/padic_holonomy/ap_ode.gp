/* ap_rank.gp -- TASK 2: Q(x)-linear independence rank certificate for (p,k)=(5,1).
   f1=1, f2=H, f3=H', f4=H'', f5=H''' ; family {x^j f_i, 0<=j<=D}.
   Rank over F_pr, pr = 2^61-1 (Mersenne prime).  Exact linear algebra over F_pr;
   a full rank over F_pr certifies full rank over Q (hence Q(x)-independence of
   the f_i up to degree D).
*/
default(parisizemax, 12000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/ap_lib.gp");
pr = 2^61 - 1;
print("pr = 2^61-1 prime? ", isprime(pr));
NN = 900;
gettime();
rr = buildmod(NN, pr);
aa = rr[1]; bb = rr[2];
print("built a_n,b_n mod pr up to n=", NN, "  ms=", gettime());
ref = read("/home/ubuntu/code/math-modular-sources/lattice/padic_holonomy/refdata.gp");
okr = 1;
{ for(i=1, 41, if(aa[i] != Mod(ref[1][i], pr) || bb[i] != Mod(ref[2][i], pr), okr = 0)); }
print("mod-pr a_n,b_n match exact rational reference for n<=40: ", okr);

/* fv[i] = coefficient vector of f_i, fv[i][n] = [x^(n-1)] f_i, length len */
{
famvecs(etaval, nfun, len) =
  my(hh, dd, res, ev);
  ev = Mod(etaval, pr);
  hh = vector(len+8, n, aa[n] + ev*bb[n]);
  res = vector(nfun);
  res[1] = vector(len+8, n, if(n==1, Mod(1,pr), Mod(0,pr)));
  dd = hh;
  for(i = 2, nfun,
    res[i] = dd;
    dd = vector(len+8, n, if(n < len+8, n*dd[n+1], Mod(0,pr)));
  );
  res;
}

{
ranktable(etaval, nfun, dmax, ds) =
  my(tt, fv, cols, out, k);
  tt = nfun*(dmax+1) + 200;
  fv = famvecs(etaval, nfun, tt);
  cols = vector(nfun*(dmax+1));
  k = 0;
  for(j = 0, dmax,
    for(i = 1, nfun,
      k++;
      cols[k] = vectorv(tt, m, if(m-1 >= j, fv[i][m-j], Mod(0,pr)));
    );
  );
  out = List();
  for(u = 1, #ds,
    my(dd, need, rk);
    dd = ds[u]; need = nfun*(dd+1);
    rk = matrank(Mat(vector(need, s, cols[s])));
    listput(out, [dd, need, rk, need-rk]);
  );
  Vec(out);
}

/* inspect the unique relation at D=8 in the 5-function family: is it really
   an inhomogeneous ODE of order 3 (coefficient of H''' nonzero)?          */
{
odeinfo(etaval) =
  my(nfun, dmax, tt, fv, cols, mt, kr, v, cf, dg);
  nfun = 5; dmax = 8; tt = 400;
  fv = famvecs(etaval, nfun, tt);
  cols = vector(nfun*(dmax+1));
  my(k = 0);
  for(j = 0, dmax, for(i = 1, nfun, k++;
    cols[k] = vectorv(tt, m, if(m-1 >= j, fv[i][m-j], Mod(0,pr)))));
  mt = Mat(cols);
  kr = matker(mt);
  print("  kernel dimension at D=8: ", #kr);
  if(#kr == 1,
    v = kr[,1];
    /* coefficient polynomial of f_i is sum_j v[j*5+i] x^j */
    cf = vector(nfun, i, sum(j = 0, dmax, if(v[j*nfun+i] == 0, 0, 1)*x^j));
    dg = vector(nfun, i, my(dz = -1); for(j = 0, dmax, if(v[j*nfun+i] != 0, dz = j)); dz);
    print("  deg of the polynomial coefficient of [1, H, H', H'', H''']: ", dg);
  );
}
print(""); print("--- ODE structure, eta = 5/7 ---"); odeinfo(5/7);
print("--- ODE structure, eta = 0 ---");   odeinfo(0);
quit;
