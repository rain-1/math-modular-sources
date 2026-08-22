\\ sc_rows.gp -- common row definitions for the supercongruence study.
\\ Rows R1..R5, A-sequence and companions X^(j), all exact over Q.
\\ (parisizemax set by the driver script)

\\ ---------------------------------------------------------------
\\ Generic driver.  A row is given by an order-r recurrence
\\      sum_{i=0}^{r} Cf[i+1](n) * u_{n-i} = 0      (valid for n >= 1)
\\ with Cf[1](n) the leading coefficient (of u_n), nonvanishing for n>=1.
\\ genseq(Cf, r, init, N) : init = vector of the first #init values
\\   u_0..u_{#init-1}; the recurrence is used for n >= #init.
\\ ---------------------------------------------------------------
genseq(Cf, r, init, N) = {
  my(u = vector(N+1));
  for(i=1, #init, u[i] = init[i]);
  for(n = #init, N,
    my(s = 0);
    for(i = 1, r, if(n-i >= 0, s += subst(Cf[i+1], 'n, n) * u[n-i+1]));
    u[n+1] = -s / subst(Cf[1], 'n, n);
  );
  u;
};

\\ companion X^(j): x_i = 0 for i<j, x_j = 1, recurrence for n > j
compan(Cf, r, j, N) = {
  my(init = vector(j+1)); init[j+1] = 1;
  genseq(Cf, r, init, N);
};

\\ ---------------------------------------------------------------
\\ R1: Apery zeta(3).  (n+1)^3 u_{n+1} = (2n+1)(17n^2+17n+5) u_n - n^3 u_{n-1}
\\ shift n -> n-1 :  n^3 u_n - (2n-1)(17n^2-17n+5) u_{n-1} + (n-1)^3 u_{n-2} = 0
\\ ---------------------------------------------------------------
R1cf = [ 'n^3, -(2*'n-1)*(17*'n^2-17*'n+5), ('n-1)^3 ];
R1r  = 2; R1k = 3;

\\ R2: Zagier C (10,3,9): (n+1)^2 u_{n+1} = (10n^2+10n+3)u_n - 9n^2 u_{n-1}
R2cf = [ 'n^2, -(10*'n^2-10*'n+3), 9*('n-1)^2 ];
R2r  = 2; R2k = 2;

\\ R3: AZ eta (11,5,125): (n+1)^3 u_{n+1} = (2n+1)(11n^2+11n+5)u_n - 125 n^3 u_{n-1}
R3cf = [ 'n^3, -(2*'n-1)*(11*'n^2-11*'n+5), 125*('n-1)^3 ];
R3r  = 2; R3k = 3;

\\ R5: AESZ 207.  L = sum_i z^i Q_i(theta); recurrence sum_i Q_i(n-i) u_{n-i} = 0
{ R5Q = [ 'n^4,
        -2^4*(1072*'n^4 - 17824*'n^3 - 10888*'n^2 - 1976*'n - 145),
        -2^17*(51088*'n^4 + 116368*'n^3 - 45264*'n^2 - 14228*'n - 1397),
        2^28*13*(73104*'n^4 + 1536*'n^3 - 488*'n^2 + 384*'n + 97),
        -2^44*13^2*(2*'n+1)^4 ]; }
R5cf = vector(5, i, subst(R5Q[i], 'n, 'n-(i-1)));
R5r  = 4; R5k = 4;

\\ R4: Sym^3 of Zagier E.  A_n = [t^n] F(t)^3, F from (n+1)^2 e_{n+1}=(12n^2+12n+4)e_n-32n^2 e_{n-1}
Ecf = [ 'n^2, -(12*'n^2-12*'n+4), 32*('n-1)^2 ];
r4A(N) = {
  my(e = genseq(Ecf, 2, [1,4], N+2));
  my(F = sum(i=0, N, e[i+1]*'t^i) + O('t^(N+1)));
  my(G = F^3);
  vector(N+1, i, polcoeff(G, i-1, 't));
};

\\ fit a recurrence sum_{i=0}^{ord} P_i(n) u_{n-i} = 0 of poly degree d
scfit(A, ord, d) = {
  my(rows = List(), nmax = min(#A-1, ord + (ord+1)*(d+1) + 12));
  for(n = ord, nmax,
    my(rr = List());
    for(i = 0, ord, for(e = 0, d, listput(rr, n^e * A[n-i+1])));
    listput(rows, Vec(rr));
  );
  my(mx = matconcat(vector(#rows, i, Mat(rows[i]))~));
  my(kk = matker(mx));
  if(#kk == 0, return(0));
  if(#kk > 1, return(-1));
  my(v = kk[,1]); v = v/content(v);
  vector(ord+1, j, Polrev(vector(d+1, e, v[(j-1)*(d+1)+e]), 'n));
};

\\ R4: minimal recurrence for A_n = [t^n] F(t)^3 (fitted here, order 6, degree 4;
\\ verified 0 failures n=6..148, and no order<=5 / degree<=14 relation exists).
{ R4cf = [ 'n^4,
  -36*'n^4 + 72*'n^3 - 88*'n^2 + 52*'n - 12,
  528*'n^4 - 2112*'n^3 + 3952*'n^2 - 3680*'n + 1392,
  -4032*'n^4 + 24192*'n^3 - 61120*'n^2 + 74496*'n - 36288,
  16896*'n^4 - 135168*'n^3 + 430592*'n^2 - 641024*'n + 373248,
  -36864*'n^4 + 368640*'n^3 - 1417216*'n^2 + 2478080*'n - 1658880,
  32768*('n-3)^4 ]; }
R4r = 6;

\\ measure the sharp denominator exponent k for a rational sequence x
denexp(x, N) = {
  my(k = 0, L = 1);
  for(n = 1, N,
    L = lcm(L, n);
    my(d = denominator(x[n+1]));
    if(d > 1,
      my(kk = 0, dd = d);
      while(dd > 1, kk++; dd = dd/gcd(dd, L); if(kk > 20, break));
      if(kk > k, k = kk);
    );
  );
  k;
};
