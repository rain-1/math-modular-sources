\\ /home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/lib12.gp
\\ Conductor-12 well-poised rows for L(2,chi_{-4}) = Catalan's constant G.
\\ Exact rational arithmetic.  Names avoid PARI builtins.
\\
\\ Poles of order 2 at t = -(j + r/12), r in {1,5,7,11}, j = 0..b:  D_b of degree 8b+8.
\\ Numerator  N(t) = (2t+b+1) * prod_{j=1..a}(t-j+1)(t+b+j) * (u^4 + c1 u^2 + c0),
\\ u = t+(b+1)/2, so N is odd in u, i.e. N(-t-b-1) = -N(t)  (antisymmetry).
\\ The involution sigma(t) = -t-b-1 maps the class r to the class 12-r, so
\\   sum A^{(11)} = -sum A^{(1)},  sum A^{(7)} = -sum A^{(5)},
\\   sum B^{(11)} = +sum B^{(1)},  sum B^{(7)} = +sum B^{(5)}.
\\ deg N <= deg D - 2 forces sum_i B_i = 0, hence sum B^{(5)} = -sum B^{(1)}.
\\ Two free parameters (c1,c0), two linear conditions:
\\   (C1)  sum B^{(1)} = 0            (kills all psi(r/12), i.e. all logs and pi cot)
\\   (C2)  sum A^{(1)} - sum A^{(5)} = 0   (kills L(2,chi_{-3}))
\\ Then, with  zeta(2,1/12)-zeta(2,11/12) = 72(L_12(2,chi_-4)+L_12(2,chi_-3)),
\\             zeta(2,5/12)-zeta(2,7/12) = 72(L_12(2,chi_-4)-L_12(2,chi_-3)),
\\             L_12(2,chi_-4) = (1-chi_-4(3)/9) G = (10/9) G,
\\   sum_{t>=0} R(t) = 144 * (10/9) * sum A^{(1)} * G - P = 160 * sum A^{(1)} * G - P.

\\ raw data for a given (a,b) and a given even quartic q(u)=u^4+c1 u^2+c0:
\\ returns [avec, bvec] over the nn = 4(b+1) poles, ordered class-major r=1,5,7,11.
c12raw(aa, bb, c1, c0) =
{ my(rs = [1,5,7,11], hh = 4, nn = 4*(bb+1), po = vector(nn),
     av = vector(nn), bv = vector(nn), hb = (bb+1)/2);
  for(c = 1, hh, for(j = 0, bb, po[(c-1)*(bb+1)+j+1] = -(j + rs[c]/12)));
  for(i = 1, nn,
    my(p = po[i], u = p + hb, qq, qd, val, dlg, pr = 1, sm = 0);
    qq = u^4 + c1*u^2 + c0;
    qd = 4*u^3 + 2*c1*u;
    val = (2*p + bb + 1) * qq;
    dlg = 2/(2*p + bb + 1) + qd/qq;
    for(j = 1, aa, my(f1 = p-j+1, f2 = p+bb+j); val *= f1*f2; dlg += 1/f1 + 1/f2);
    for(i2 = 1, nn, if(i2 != i, my(d = p - po[i2]); pr *= d^2; sm += 1/d));
    av[i] = val/pr;
    bv[i] = av[i]*(dlg - 2*sm));
  [av, bv];
};

\\ the two linear functionals, as functions of (c1,c0), are affine; solve exactly.
\\ Returns [c1, c0].
c12cond(aa, bb, c1, c0) =
{ my(z = c12raw(aa,bb,c1,c0), av = z[1], bv = z[2], nb = bb+1);
  [ sum(j=1,nb, bv[j]), sum(j=1,nb, av[j]) - sum(j=nb+1,2*nb, av[j]) ];
};

c12solve(aa, bb) =
{ my(v0, v1, v2, mm, rhs, so);
  v0 = c12cond(aa,bb,0,0);
  v1 = c12cond(aa,bb,1,0) - v0;
  v2 = c12cond(aa,bb,0,1) - v0;
  mm = [v1[1], v2[1]; v1[2], v2[2]];
  rhs = [-v0[1], -v0[2]]~;
  if(matdet(mm) == 0, error("c12solve: degenerate at a=",aa," b=",bb));
  so = matsolve(mm, rhs);
  [so[1], so[2]];
};

\\ full row: returns [Q, P, c1, c0] with sum_{t>=0} R(t) = Q*G - P.
c12row(aa, bb) =
{ my(cc = c12solve(aa,bb), z, av, bv, nb = bb+1, nn = 4*nb, qc, pc, rs = [1,5,7,11]);
  z = c12raw(aa, bb, cc[1], cc[2]); av = z[1]; bv = z[2];
  if(sum(i=1,nn, bv[i]) != 0, error("c12: total sum B nonzero"));
  if(sum(j=1,nb, bv[j]) != 0, error("c12: sum B^(1) nonzero"));
  if(sum(j=nb+1,2*nb, bv[j]) != 0, error("c12: sum B^(5) nonzero"));
  if(sum(i=1,nn, av[i]) != 0, error("c12: total sum A nonzero"));
  if(sum(j=1,nb, av[j]) != sum(j=nb+1,2*nb, av[j]), error("c12: chi_-3 part nonzero"));
  qc = 160*sum(j=1,nb, av[j]);
  pc = 0;
  for(c = 1, 4,
    my(base = rs[c]/12, t1 = 0, t2 = 0);
    for(j = 0, bb,
      my(i = (c-1)*nb + j + 1);
      pc += av[i]*t2 + bv[i]*t1;
      t2 += 1/(j+base)^2; t1 += 1/(j+base)));
  [qc, pc, cc[1], cc[2]];
};

\\ direct numerical check of sum_{t>=0} R(t)
c12direct(aa, bb, c1, c0, tmax) =
{ my(rs = [1,5,7,11], hb = (bb+1)/2, tot = 0.0);
  for(t = 0, tmax,
    my(u = t + hb, nu = (2*t+bb+1)*(u^4 + c1*u^2 + c0), de = 1.0);
    for(j = 1, aa, nu *= (t-j+1)*(t+bb+j));
    for(c = 1, 4, for(j = 0, bb, de *= (t + j + rs[c]/12)^2));
    tot += nu/de);
  tot;
};

\\ ===================================================================
\\ CLEANER DESIGN (no linear solve): split the four pole classes into the
\\ two sigma-pairs {1/12,11/12} and {5/12,7/12}.  Each half is a well-poised
\\ antisymmetric row on its own, so sum B vanishes in each class automatically;
\\ the only remaining condition (kill L(2,chi_-3)) is met by cross-scaling.
\\ half(aa,bb,r) : poles at -(j+r/12), -(j+(12-r)/12), j=0..bb; numerator
\\   N(t) = (2t+bb+1) prod_{j=1..aa}(t-j+1)(t+bb+j).
\\ returns [sumA_firstclass, Ppart] with
\\   sum_{t>=0} R(t) = sumA * (zeta(2,r/12)-zeta(2,(12-r)/12)) - Ppart.
c12half(aa, bb, rr) =
{ my(nn = 2*bb+2, po = vector(nn), av = vector(nn), bv = vector(nn),
     sa1, sa2, sb1, pc, t1, t2, base);
  for(j = 0, bb, po[j+1] = -(j + rr/12); po[bb+2+j] = -(j + (12-rr)/12));
  for(i = 1, nn,
    my(p = po[i], val = 2*p + bb + 1, dlg = 2/(2*p + bb + 1), pr = 1, sm = 0);
    for(j = 1, aa, my(f1 = p-j+1, f2 = p+bb+j); val *= f1*f2; dlg += 1/f1 + 1/f2);
    for(i2 = 1, nn, if(i2 != i, my(d = p - po[i2]); pr *= d^2; sm += 1/d));
    av[i] = val/pr;
    bv[i] = av[i]*(dlg - 2*sm));
  if(sum(i=1,nn, bv[i]) != 0, error("c12half: sum B nonzero"));
  sa1 = sum(j=1,bb+1, av[j]); sa2 = sum(j=bb+2,nn, av[j]);
  sb1 = sum(j=1,bb+1, bv[j]);
  if(sa1 + sa2 != 0, error("c12half: symmetric A part nonzero"));
  if(sb1 != 0, error("c12half: psi coefficient nonzero"));
  pc = 0;
  for(cl = 0, 1,
    base = if(cl == 0, rr/12, (12-rr)/12); t1 = 0; t2 = 0;
    for(j = 0, bb,
      my(i = cl*(bb+1) + j + 1);
      pc += av[i]*t2 + bv[i]*t1;
      t2 += 1/(j+base)^2; t1 += 1/(j+base)));
  [sa1, pc];
};

\\ full conductor-12 Catalan row, product normalisation (no denominators created):
\\   Q = 160 * s1 * s5,   P = s5*P1 + s1*P5,   sum R = Q*G - P.
c12rowB(aa1, aa5, bb) =
{ my(h1 = c12half(aa1,bb,1), h5 = c12half(aa5,bb,5));
  [ 160*h1[1]*h5[1], h5[1]*h1[2] + h1[1]*h5[2] ];
};
