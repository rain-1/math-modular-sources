\\ Verification of the two headline second-order rows and of the CDT inputs that
\\ this note uses:  a_n in Z, d_n^2 b_n in Z (k = 2 sharp), the characteristic
\\ roots, x_2 = 1/lam_2, and the radius of convergence of H = xi*A - B.
\\ (PARI builtin names psi, M, Phi, S, cmp deliberately avoided.)
default(realprecision, 60);

mkrow(alv, bev, gav, dev, epv, zev, NN, a0, a1) =
{ my(v = vector(NN+1)); v[1] = a0; v[2] = a1;
  for(n = 1, NN-1,
    v[n+2] = ((alv*n^2+bev*n+gav)*v[n+1] - (dev*n^2+epv*n+zev)*v[n]) / (n+1)^2);
  v; }

dnv(NN) = { my(w = vector(NN+1)); w[1] = 1;
            for(n = 1, NN, w[n+1] = lcm(w[n], n)); w; }

checkrow(nm, alv, bev, gav, dev, epv, zev, NN) =
{ my(av, bv, dn, kk, ok, l1, l2, rts);
  av = mkrow(alv,bev,gav,dev,epv,zev,NN,1,gav);
  bv = mkrow(alv,bev,gav,dev,epv,zev,NN,0,1);
  dn = dnv(NN);
  ok = 1; for(n=0,NN, if(denominator(av[n+1]) != 1, ok = 0));
  kk = 0; while(kk <= 6,
    if(sum(n=0,NN, if(denominator(dn[n+1]^kk*bv[n+1]) != 1, 1, 0)) == 0, break); kk++);
  rts = polroots(x^2 - alv*x + dev);
  l1 = rts[1]; l2 = rts[2];
  if(abs(l1) < abs(l2), [l1,l2] = [l2,l1]);
  printf("%-32s a_n in Z (n<=%d): %d   k = %d (sharp: d_n^%d b_n not in Z: %d)\n",
         nm, NN, ok, kk, kk-1,
         sum(n=0,NN, if(denominator(dn[n+1]^(kk-1)*bv[n+1]) != 1, 1, 0)) > 0);
  printf("      lam_1 = %.12f   lam_2 = %.12f   lam_1*lam_2 = %d\n", real(l1), real(l2), dev);
  printf("      x_2 = 1/lam_2 = %.12f   4|x_2| = %.12f   log 4|x_2| = %.10f\n",
         real(1/l2), abs(4/l2), log(abs(4/l2)));
  printf("      score = log(1/|lam_2|) - 2 = %+.6f\n", log(1/abs(l2)) - 2);
  \\ the Apery limit and the decay rate of the linear form
  my(xi = bv[NN+1]/av[NN+1], rr = vector(20));
  for(j=1,20, rr[j] = abs(av[NN-20+j]*xi - bv[NN-20+j])^(1.0/(NN-21+j)));
  printf("      xi = %.30f\n", xi);
  printf("      |a_n xi - b_n|^(1/n) at n=%d : %.8f   (should be |lam_2| = %.8f)\n\n",
         NN-1, rr[19], abs(l2));
}

print("=== second-order rows: CDT inputs ===");
checkrow("Beukers sqrt-Apery",       136,68,10,16,-16,4, 160);
checkrow("level-5 Gamma_0(5)+5",      88,44, 6,-64,64,-12, 160);
checkrow("sqrt-T",                    24,12, 2,16,-16,4, 160);
checkrow("sqrt-Domb",                 20,10, 2,64,-64,16, 160);
checkrow("sqrt-Cooper s_7",           26,13, 2,-27,27,-6, 160);
checkrow("Apery zeta(2) (calib)",     11,11, 3,-1,0,0, 160);
checkrow("Zagier C = CDT (calib)",    10,10, 3,9,0,0, 160);
checkrow("Zagier E = Catalan (calib)",12,12, 4,32,0,0, 160);

print("=== geometric denominator: the level-5 modular coordinate ===");
print("t = 2x, so |t_2| = 2|x_2|; N1: score = log|t_2| - 2 - log 2 = log|x_2| - 2.");
{ my(x2 = 1/(44-20*sqrt(5)));
  printf("  |x_2| = %.12f, |t_2| = %.12f, log|t_2|-2-log(2) = %+.6f\n",
         abs(x2), 2*abs(x2), log(2*abs(x2))-2-log(2)); }
{ my(x2 = 1/(4*(17-12*sqrt(2))));
  printf("  Beukers: |x_2| = %.12f, |t_2| = 4|x_2| = %.12f = (1+sqrt2)^4 = %.12f\n",
         abs(x2), 4*abs(x2), (1+sqrt(2))^4);
  printf("           score = log|t_2|-2-log(4) = %+.6f\n", log(4*abs(x2))-2-log(4)); }
