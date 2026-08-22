default(parisizemax, 4000000000);
default(parisize, 1000000000);
default(realprecision, 250);

\\ ---------------------------------------------------------------
\\ Fold constant xi = C_minus/A_minus for second-order Apery-like rows
\\ (n+1)^2 u_{n+1} = P(n) u_n - Q(n) u_{n-1},  P=a n^2+b n+c, Q=d n^2+e n+f
\\ ---------------------------------------------------------------

\\ Frobenius recursion solver.
\\ al[j+1]=alpha_j (j=0..4), be[j+1]=beta_j (j=0..3), ga[j+1]=gamma_j (j=0..2)
\\ rr[m+1] = coefficient of s^(sig+m-1) in the RHS, m=0..MM
\\ returns y[m+1] = y_m , with y_0 = y0 imposed.
frob(al,be,ga, sig, y0, rr, MM) = {
  my(y = vector(MM+1), kn, dd);
  y[1] = y0;
  for(m=1, MM,
    kn = 0;
    for(j=2,4, if(m+1-j>=0, kn += al[j+1]*(sig+m+1-j)*(sig+m-j)*y[m+1-j+1]));
    for(j=1,3, if(m-j>=0,   kn += be[j+1]*(sig+m-j)*y[m-j+1]));
    for(j=0,2, if(m-1-j>=0, kn += ga[j+1]*y[m-1-j+1]));
    dd = al[2]*(sig+m)*(sig+m-1) + be[1]*(sig+m);
    y[m+1] = (rr[m+1] - kn)/dd;
  );
  y;
}

\\ evaluate sum y_m s^m and its s-derivative, plus 2nd derivative of s^sig*sum
evser(y, sig, s, MM) = {
  my(v0=0., v1=0., v2=0., sp=1.);
  for(m=0, MM,
    v0 += y[m+1]*sp;
    v1 += (sig+m)*y[m+1]*sp;
    v2 += (sig+m)*(sig+m-1)*y[m+1]*sp;
    sp *= s;
  );
  \\ f = s^sig*v0 ; f' = s^(sig-1)*v1 ; f'' = s^(sig-2)*v2
  [s^sig*v0, s^(sig-1)*v1, s^(sig-2)*v2];
}

dofold(nam, a,b,c,d,e,f, tc, NA, MM, NEST) = {
  my(x='x, Alp,Bep,Gap, al,be,ga, rho, hh, sm, tm, uu,vv,pp,ww, rr, U,V,PP, Ap,Am,Cp,Cm,xi);
  print("==================================================================");
  print("ROW: ", nam);
  print("  a,b,c,d,e,f = ", [a,b,c,d,e,f]);
  Alp = subst(x^2 - a*x^3 + d*x^4, x, tc+x);
  Bep = subst(x - (a+b)*x^2 + (3*d+e)*x^3, x, tc+x);
  Gap = subst(-c*x + (d+e+f)*x^2, x, tc+x);
  al = vector(5, j, polcoef(Alp, j-1));
  be = vector(4, j, polcoef(Bep, j-1));
  ga = vector(3, j, polcoef(Gap, j-1));
  print("  |alpha_0| (should be ~0): ", abs(al[1]));
  al[1] = 0;
  rho = 1 - be[1]/al[2];
  print("  t_c   = ", tc);
  print("  |t_c| = ", abs(tc), "   1/|t_c| = ", 1/abs(tc));
  print("  rho   = ", rho);
  print("  alpha_1 = ", al[2]);

  \\ --- convergence radii ---
  my(rts = polroots(d*x^2 - a*x + 1), other = 0);
  for(i=1,#rts, if(abs(rts[i]-tc) > 1e-100, other = rts[i]));
  my(radloc = min(abs(tc), abs(tc-other)));
  hh = 0.4*abs(tc);
  sm = -hh*tc/abs(tc);
  tm = tc + sm;
  print("  other root = ", other, "   local radius = ", radloc);
  print("  h = ", hh, "  ratio_origin = ", abs(tm)/abs(tc), "  ratio_local = ", hh/radloc);
  print("  ratio_local^MM = ", (hh/radloc)^MM, "   ratio_origin^NA = ", (abs(tm)/abs(tc))^NA);

  \\ --- local solutions ---
  rr = vector(MM+1);
  uu = frob(al,be,ga, 0, 1, rr, MM);
  my(islog = (abs(rho) < 1e-150));
  if(islog,
    print("  *** LOGARITHMIC CASE (rho = 0) ***");
    \\ v = u*log(s) + w ; L[w] = -( alpha*(2u'/s - u/s^2) + beta*u/s )
    rr = vector(MM+1);
    for(m=0, MM,
      my(t0=0);
      for(j=1,4, if(m+1-j>=0, t0 += al[j+1]*(2*(m+1-j)-1)*uu[m+1-j+1]));
      for(j=0,3, if(m-j>=0,   t0 += be[j+1]*uu[m-j+1]));
      rr[m+1] = -t0;
    );
    print("  rr_w[0] (should be ~0): ", abs(rr[1]));
    rr[1] = 0;
    ww = frob(al,be,ga, 0, 0, rr, MM);
  ,
    rr = vector(MM+1);
    vv = frob(al,be,ga, rho, 1, rr, MM);
  );
  \\ particular solution p, L[p] = t = tc + s
  rr = vector(MM+1);
  if(MM>=1, rr[2] = tc);
  if(MM>=2, rr[3] = 1);
  pp = frob(al,be,ga, 0, 0, rr, MM);

  \\ --- evaluate at s = sm ---
  my(lg = log(abs(sm)) + I*arg(sm*(1+0.*I)));
  my(evalU(s) = evser(uu,0,s,MM));
  my(evalP(s) = evser(pp,0,s,MM));
  my(evalV(s) = if(islog,
        my(EU=evser(uu,0,s,MM), EW=evser(ww,0,s,MM), L=log(abs(s))+I*arg(s*(1+0.*I)));
        [EU[1]*L + EW[1], EU[2]*L + EU[1]/s + EW[2], EU[3]*L + 2*EU[2]/s - EU[1]/s^2 + EW[3]]
      , evser(vv,rho,s,MM)));

  \\ residual checks
  my(Lop(EE, tt) = (tt^2*(1-a*tt+d*tt^2))*EE[3] + (tt*(1-(a+b)*tt+(3*d+e)*tt^2))*EE[2] + (-c*tt+(d+e+f)*tt^2)*EE[1]);
  my(EU=evalU(sm), EV=evalV(sm), EP=evalP(sm));
  print("  residual L[u]     = ", abs(Lop(EU,tm)));
  print("  residual L[v]     = ", abs(Lop(EV,tm)));
  print("  residual L[p] - t = ", abs(Lop(EP,tm) - tm));

  \\ --- global series ---
  my(A = vector(NA+1), B = vector(NA+1));  \\ A[n+1]=a_n
  A[1] = 1; A[2] = c; B[1] = 0; B[2] = 1;
  for(n=1, NA-1,
    my(Pn = a*n^2+b*n+c, Qn = d*n^2+e*n+f);
    A[n+2] = (Pn*A[n+1] - Qn*A[n])/(n+1)^2;
    B[n+2] = (Pn*B[n+1] - Qn*B[n])/(n+1)^2;
  );
  print("  a_n, n=0..19: ", vector(20,i,A[i]));
  print("  b_n, n=0..19: ", vector(20,i,B[i]));

  my(y0=0., y0d=0., yB=0., yBd=0., pw=1.);
  for(n=0, NA,
    my(an = A[n+1]*1., bn = B[n+1]*1.);
    y0 += an*pw; yB += bn*pw;
    if(n>=1, y0d += n*an*pw/tm; yBd += n*bn*pw/tm);
    pw *= tm;
  );
  print("  last term magnitude: ", abs(A[NA+1]*1.*tm^NA));

  \\ --- solve ---
  my(det = EU[1]*EV[2] - EU[2]*EV[1]);
  Ap = ( y0*EV[2] - y0d*EV[1] )/det;
  Am = ( EU[1]*y0d - EU[2]*y0 )/det;
  my(g = yB - EP[1], gd = yBd - EP[2]);
  Cp = ( g*EV[2] - gd*EV[1] )/det;
  Cm = ( EU[1]*gd - EU[2]*g )/det;
  xi = Cm/Am;
  write("conn.txt", Ap); write("conn.txt", Am); write("conn.txt", Cp); write("conn.txt", Cm);
  print("  A_plus  = ", Ap);
  print("  A_minus = ", Am);
  print("  C_plus  = ", Cp);
  print("  C_minus = ", Cm);

  \\ --- independent verification at a second point ---
  my(s2 = -0.28*abs(tc)*tc/abs(tc), t2 = tc + s2);
  my(FU=evalU(s2), FV=evalV(s2), FP=evalP(s2));
  my(z0=0., zB=0., pw2=1.);
  for(n=0, NA, z0 += A[n+1]*1.*pw2; zB += B[n+1]*1.*pw2; pw2 *= t2);
  print("  CHECK y_0 at 2nd pt, rel err = ", abs((Ap*FU[1]+Am*FV[1] - z0)/z0));
  print("  CHECK y_B at 2nd pt, rel err = ", abs((Cp*FU[1]+Cm*FV[1]+FP[1] - zB)/zB));

  print("  THETA* = C_plus/A_plus = ", Cp/Ap);
  print("  XI = ", xi);
  print("  Re(xi) = ", real(xi));
  print("  Im(xi) = ", imag(xi));

  \\ --- cross-check estimator ---
  if(NEST > 0,
    my(lam = 1/tc, cl = conj(lam), acc = 0., cnt = 0);
    for(n=NEST\2, NEST,
      my(num = B[n+2]*1. - cl*B[n+1]*1., den = A[n+2]*1. - cl*A[n+1]*1.);
      acc += num/den; cnt++;
    );
    my(est = acc/cnt);
    print("  estimator (avg n=",NEST\2,"..",NEST,") = ", est);
    print("    |est - xi|       = ", abs(est - xi));
    print("    |est - conj(xi)| = ", abs(est - conj(xi)));
  );
  xi;
}
