/* 10_domb_ident.gp -- high-precision confirmation of the two new identifications
   found in the sqrt(Domb) cusp-move orbit:
      sqrtDomb.5 = (20,10,64), P = 20n^2+26n+10, Q = 64n^2 :
          xi = (15 L(2,chi_-3) - 6 zeta(2))/16
      sqrtDomb.7 = (20, 4,64), P = 20n^2+14n+4,  Q = 64n^2 :
          xi = ( 6 zeta(2) + 15 L(2,chi_-3))/32                                   */
read("lib.gp");
default(realprecision, 260);
zt2 = zeta(2); lch = lfun(-3,2);
{ tst(nm, pp, qq, target, nmax) =
  my(av = seqA([pp,qq], nmax), bv = seqB([pp,qq], nmax), x, kk, dn, ok);
  x = bv[nmax+1]/av[nmax+1]*1.0;
  print("  ", nm, "  P=", pp, "  Q=", qq);
  print("     xi (n=", nmax, ") = ", strprintf("%.60f", x));
  print("     target            = ", strprintf("%.60f", target));
  print("     |xi - target|     = ", strprintf("%.3e", abs(x-target)));
  print("     predicted decay |lam2/lam1|^n = ", strprintf("%.3e", (4.0/16)^nmax));
  dn = dnvec(nmax); ok = 1;
  for(n=0,nmax, if(denominator(av[n+1])!=1, ok=0; break));
  print("     a_n integral to n=", nmax, ": ", ok, "   k = ", denomexp(bv, nmax, 5), " (sharp)");
}
print("=== high-precision identification, sqrt(Domb) orbit ===");
tst("sqrtDomb.5", 20*nv^2+26*nv+10, 64*nv^2, (15*lch - 6*zt2)/16, 200);
tst("sqrtDomb.7", 20*nv^2+14*nv+4,  64*nv^2, ( 6*zt2 + 15*lch)/32, 200);
print("");
print("  parent sqrtDomb.1 = (20,2,64): xi = L(f_12,2), not an Eisenstein value:");
tst("sqrtDomb.1", 20*nv^2+10*nv+2, 64*nv^2-64*nv+16, 0.0, 200);
print("");
print("=== the two forms in 1, zeta(2), L(2,chi-3) they produce ===");
print(" row .5 : 16 d_n^2 (a_n xi - b_n) = d_n^2 a_n (15 L - 6 zeta(2)) - 16 d_n^2 b_n");
print(" row .7 : 32 d_n^2 (a_n xi - b_n) = d_n^2 a_n ( 6 zeta(2) + 15 L) - 32 d_n^2 b_n");
print(" coefficient vectors in the (zeta(2), L)-plane: (-6,15) and (6,15) -- independent.");
{ my(sg = 2+log(16), dl = -2-log(4));
  print(" sigma = k + log lam1 = ", sg, "   delta = -k - log|lam2| = ", dl);
  print(" tau = delta/sigma = ", dl/sg, "   Nesterenko dim >= ", 1+dl/sg);
  print(" deficit to dim>=2 : ", -dl, " nats/step;  to dim>=3 : ", 2*sg-dl, " nats/step"); }
quit;
