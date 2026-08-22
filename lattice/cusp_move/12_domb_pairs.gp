/* 12_domb_pairs.gp -- are the archimedean periods of the p-adically related
   sqrt(Domb) placements rationally proportional?  (they are not)             */
read("lib.gp");
default(realprecision, 300);
{ xi(pp, qq, nmax) = my(av=seqA([pp,qq],nmax), bv=seqB([pp,qq],nmax)); bv[nmax+1]/av[nmax+1]*1.0; }
x1 = xi(20*nv^2+10*nv+2,   64*nv^2-64*nv+16, 400);
x2 = xi(-28*nv^2-22*nv-6,  192*nv^2-96*nv,   700);
x3 = xi(8*nv^2+2*nv,       -48*nv^2+24*nv,   400);
x4 = xi(8*nv^2+14*nv+6,    -48*nv^2-24*nv,   400);
x5 = xi(20*nv^2+26*nv+10,  64*nv^2,          400);
x6 = xi(-28*nv^2-34*nv-12, 192*nv^2+96*nv,   700);
x7 = xi(20*nv^2+14*nv+4,   64*nv^2,          400);
{ tst(nm, u, w) = my(r = u/w, q = bestappr(r, 10^15));
  print("  ", nm, " = ", strprintf("%.40f", r), "   best rational (height<=10^15) ", q,
        "   err ", strprintf("%.3e", abs(r-q*1.0))); }
print("=== archimedean ratios of the p-adically related sqrt(Domb) placements ===");
print("  (p-adic:  xi_p(2) = -2 xi_p(3) at p=2 and p=3;  xi_p(4) = xi_p(6);  xi_2(5) = 2 xi_2(7))");
tst("xi(2)/xi(3)", x2, x3);
tst("xi(4)/xi(6)", x4, x6);
tst("xi(5)/xi(7)", x5, x7);
print("  for reference, xi(5) = (15 L(2,chi-3) - 6 zeta(2))/16, xi(7) = (6 zeta(2) + 15 L(2,chi-3))/32,");
print("  so xi(5)/xi(7) = 2(15L-6z)/(6z+15L) is irrational -- consistent.");
quit;
