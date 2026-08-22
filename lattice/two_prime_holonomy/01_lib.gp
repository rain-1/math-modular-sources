\\ NB: parisizemax is set by the *calling* scripts (01_*.gp); putting the
\\ default(parisizemax,...) call inside a read() file aborts parsing of the rest.
\\ Exact-rational rebuild of CDT's 14 functions (port of
\\ /home/ubuntu/code/math-modular-sources/lattice/cdt_finder/indep_check2.py
\\ from mod 2^61-1 to exact Q, so that p-adic valuations survive).

\\ ---- CDT Prop 11.1.4 inhomogeneous ODE, exact port of solve_ode ----
\\ x(1-x)(1-9x) y'' + (1-20x+27x^2) y' + 3(3x-1) y = bb + cc/(1-x)
solveode(bb, cc, pr) = {
  my(ys = vector(pr+5), ss, kk, cf, rr);
  ys[1] = if(bb==0 && cc==0, 1, 0);
  for(nn = 0, pr-1,
    ss = 0;
    if(nn>=1, kk=nn-1; ss +=   1*(kk+2)*(kk+1)*ys[kk+3]);
    if(nn>=2, kk=nn-2; ss += -10*(kk+2)*(kk+1)*ys[kk+3]);
    if(nn>=3, kk=nn-3; ss +=   9*(kk+2)*(kk+1)*ys[kk+3]);
    kk=nn;              ss +=   1*(kk+1)*ys[kk+2];
    if(nn>=1, kk=nn-1; ss += -20*(kk+1)*ys[kk+2]);
    if(nn>=2, kk=nn-2; ss +=  27*(kk+1)*ys[kk+2]);
    ss += -3*ys[nn+1];
    if(nn>=1, ss += 9*ys[nn]);
    cf = (nn+1)^2;
    rr = if(nn==0, bb+cc, cc);
    ys[nn+2] = (rr - ss)/cf;
  );
  vector(pr, ii, ys[ii]);
}

\\ ---- composition with w(x) = x/(x-1) ----
\\ w^n = (-1)^n x^n (1-x)^{-n}; [x^k] w^n = (-1)^n binom(k-1,n-1)
\\ f(w(x))[k] = sum_{n=1..k} f[n] (-1)^n binom(k-1,n-1);  [x^0] = f[0]
composew(fv, nx) = {
  my(rv = vector(nx), tt, bc);
  rv[1] = fv[1];
  for(kk = 1, nx-1,
    tt = 0; bc = 1;
    for(nn = 1, kk,
      if(fv[nn+1] != 0, tt += if(nn%2, -fv[nn+1], fv[nn+1])*bc);
      bc = bc*(kk-nn)/nn;
    );
    rv[kk+1] = tt;
  );
  rv;
}

\\ ---- re-expansion in y = x^2/(x-1) (port of to_y) ----
\\ Y^n = (-1)^n x^{2n}(1-x)^{-n};  [x^k] Y^n = (-1)^n binom(k-n-1,n-1), k>=2n
toy(fv, nx, ny) = {
  my(cur = fv, gv = vector(ny), bc, sg, gn);
  gv[1] = cur[1]; cur[1] = 0;
  for(nn = 1, ny-1,
    if(2*nn > nx-1, break);
    sg = if(nn%2, -1, 1);
    gn = cur[2*nn+1]*sg;
    gv[nn+1] = gn;
    if(gn != 0,
      bc = 1;
      for(kk = 2*nn, nx-1,
        cur[kk+1] -= gn*sg*bc;
        bc = bc*(kk-nn)/(kk-2*nn+1);
      );
    );
  );
  gv;
}

\\ ---- multiply by anti(x) = x - w(x) = 2x + x^2 + x^3 + ... ----
antimul(uv, nx) = {
  my(rv = vector(nx), ps = 0);
  \\ rv[k] = 2*u[k-1] + sum_{m=0}^{k-2} u[m]
  for(kk = 1, nx-1,
    if(kk >= 2, ps += uv[kk-1]);
    rv[kk+1] = 2*uv[kk] + ps;
  );
  rv;
}

\\ ---- d/dy and antiderivative, literal ports of d_ and i_ ----
dy(fv) = vector(#fv-1, ii, ii*fv[ii+1]);
iy(fv) = concat([0], vector(#fv-1, ii, fv[ii]/ii));

\\ ---- valuation helpers ----
vp(z, pp) = if(z == 0, "oo", valuation(z, pp));

\\ least squares slope of (n, v_p(c_n)) over n in [n0,n1], skipping zeros
lsslope(vv, n0, n1) = {
  my(sx=0, sy=0, sxx=0, sxy=0, cnt=0, t);
  for(nn = n0, n1,
    t = vv[nn+1];
    if(type(t) == "t_INT" || type(t)=="t_REAL",
      sx += nn; sy += t; sxx += nn^2; sxy += nn*t; cnt++);
  );
  if(cnt < 2, return("n/a"));
  (cnt*sxy - sx*sy)/(cnt*sxx - sx^2)*1.0;
}
