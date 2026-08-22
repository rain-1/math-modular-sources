/* 15_k3jtest.gp -- the K3 window of the J-map test.
 *
 * Companion to 05_jtest.gp (deg J <= 12, rational elliptic surfaces).
 * Here the window is raised to deg J <= 24, i.e. the elliptic-K3 case
 * (sum e = 24), which is open item 4 of HERFURTNER_CLASSIFICATION.md sec. 8.
 *
 * Test.  Let q(t) = t exp(g/y0) be the canonical nome of
 *   L = theta^2 - t P(theta) + t^2 Q(theta+1),
 *   P(n) = A(n^2 + (2*mm-j1-j2)/(2*mm) n) + B,  Q(n) = C (mm n - j1)(mm n - j2).
 * L is (projectively) the Picard-Fuchs operator of an elliptic surface over
 * P^1_t with a cusp of width h at t=0 iff for some h >= 1 and some constant
 * gam the composite  Jcal(t) := J(gam q(t)^h)  is a rational function U/V of t;
 * deg Jcal = sum of the n_i over the I_n, I_n^* fibres.  Writing
 *   W := x^h * Jcal   (a power series of valuation 0 in x = t)
 * the condition is  V*W = x^h*U  with deg U, deg V <= dmax.
 *
 * Implementation.  The linear system for V alone (eliminating U: the
 * coefficient of x^e in V*W must vanish for every e outside [h, h+dmax]) is
 * solved MODULO a 61-bit prime as a filter -- reduction mod p is a ring
 * homomorphism on the p-integral coefficients involved, so the filter has no
 * false negatives -- and every surviving (h, gam) is then redone exactly over Q
 * and certified: V is fitted from 2*dmax+2 coefficients only and the identity
 * V*W - x^h*U = O(x^NTERM) is checked on all the remaining ones.
 *
 * Variable names avoid the PARI builtins psi, M, Phi, S, cmp.
 */

default(realprecision, 38);

KPRIME = 2^61 - 1;      /* Mersenne prime; larger than any denominator here */
NTERM  = 140;           /* series length; >= 2*24 + 10 by a wide margin      */
DMAX   = 24;            /* deg J window: elliptic K3                          */
HMAX   = 24;            /* cusp width at the MUM point; h <= deg J <= 24      */
GLIM   = 4096;          /* gam = +- m^(+-h), m running over 7-smooth m <= GLIM */

/* ---------- the gam candidate list -------------------------------------- */

smoothbase(lim) =
{ my(v = List());
  for(m = 1, lim,
    my(tt = m);
    foreach([2,3,5,7], pr, while(tt % pr == 0, tt = tt/pr));
    if(tt == 1, listput(v, m)));
  Vec(v);
}
GBASE = smoothbase(GLIM);

/* ---------- the classical j-function q-expansion ------------------------ */

jvec(n) =
{ my(e4, dl, jj);
  e4 = 1 + 240*sum(k = 1, n, sigma(k,3)*x^k) + O(x^(n+1));
  dl = x*prod(k = 1, n, (1 - x^k + O(x^(n+1)))^24);
  jj = x*e4^3/dl;                       /* 1 + 744 x + 196884 x^2 + ...     */
  vector(n+1, i, polcoeff(jj, i-1));
}
JV = jvec(NTERM);

/* ---------- the canonical nome ------------------------------------------ */
/* onec = 1 (or Mod(1,p)); returns q = x + ... to O(x^(N+1)).                */

nomegen(mm, j1, j2, A, B, C, N, onec) =
{ my(be, de, ep, ze, u, g, y0, gg);
  be = onec*A*(2*mm - j1 - j2)/(2*mm);
  de = onec*C*mm^2;  ep = -onec*C*mm*(j1 + j2);  ze = onec*C*j1*j2;
  u = vector(N+2); g = vector(N+2);
  u[1] = onec; g[1] = 0*onec;
  for(n = 1, N+1,
    my(P0 = onec*A*(n-1)^2 + be*(n-1) + onec*B,
       Q0 = de*(n-1)^2 + ep*(n-1) + ze,
       dP = 2*onec*A*(n-1) + be,
       dQ = 2*de*(n-1) + ep);
    u[n+1] = (P0*u[n] - if(n >= 2, Q0*u[n-1], 0))/n^2;
    g[n+1] = (dP*u[n] + P0*g[n] - if(n >= 2, dQ*u[n-1] + Q0*g[n-1], 0)
              - 2*n*u[n+1])/n^2;
  );
  y0 = sum(n = 0, N, u[n+1]*x^n) + O(x^(N+1));
  gg = sum(n = 0, N, g[n+1]*x^n) + O(x^(N+1));
  x*exp(gg/y0);
}

/* ---------- W = x^h * J(gam * q^h) as a coefficient vector -------------- */
/* zp[m+1] = (q^h)^m ; returns vector of coefficients of W, indices 0..N-1.  */

wcoeffs(zp, mmax, base, gam, N) =
{ my(acc, W);
  acc = sum(m = 0, mmax, JV[m+1]*gam^m*zp[m+1]);
  W = base*acc/gam;
  vector(N, i, polcoeff(W, i-1));
}

/* the e's for which the coefficient of x^e in V*W must vanish */
eqlist(h, dmax, N) =
{ my(v = List());
  for(e = 0, h-1, listput(v, e));
  for(e = h + dmax + 1, N-1, listput(v, e));
  Vec(v);
}

vmat(wv, es, dmax) =
  matrix(#es, dmax+1, i, j, if(es[i] >= j-1, wv[es[i] - (j-1) + 1], 0));

/* ---------- exact certification of one (h, gam) ------------------------- */
/* returns [degJ, nused, nextra, U, V] or 0.                                */

certify(mm, j1, j2, A, B, C, h, gam, dmax, N) =
{ my(qq, zz, base, mmax, zp, wv, es, nused, ker, vv, uu, gg, R, nextra, dg,
     wpol, wser);
  qq   = nomegen(mm, j1, j2, A, B, C, N, 1);
  zz   = qq^h;
  base = x^h/zz;
  mmax = (N-1)\h;
  zp = vector(mmax+1); zp[1] = 1 + O(x^(N+1));
  for(m = 1, mmax, zp[m+1] = zp[m]*zz);
  wv = wcoeffs(zp, mmax, base, gam, N);
  es = eqlist(h, dmax, N);
  nused = 2*dmax + 2;  if(nused > #es, nused = #es);
  ker = matker(vmat(wv, vector(nused, i, es[i]), dmax));
  if(#ker == 0, return(0));
  wpol = sum(i = 0, N-1, wv[i+1]*x^i);
  wser = wpol + O(x^N);
  /* the kernel is Vmin * {polys of degree <= dmax - deg Vmin}: gcd it out */
  vv = 0;
  for(c = 1, #ker, my(pz = sum(k = 0, dmax, ker[k+1, c]*x^k));
      vv = if(vv == 0, pz, gcd(vv, pz)));
  if(vv == 0, return(0));
  uu = sum(k = 0, dmax, polcoeff(vv*wpol, h+k)*x^k);
  gg = gcd(uu, vv); if(gg != 0 && poldegree(gg) > 0, uu = uu/gg; vv = vv/gg);
  /* full check on ALL coefficients up to x^(N-1), incl. those not fitted */
  R = vv*wser - x^h*uu;
  for(e = 0, N-1, if(polcoeff(R, e) != 0, return(0)));
  nextra = #es - nused;
  dg = max(poldegree(uu), poldegree(vv));
  [dg, nused, nextra, uu, vv];
}

/* ---------- the test ---------------------------------------------------- */
/* returns [h, degJ, gam, nused, nextra] or 0.                              */

k3test(mm, j1, j2, A, B, C, dmax = DMAX, hmax = HMAX, N = NTERM, p = KPRIME) =
{ my(onec, qq, glist, nq1);
  onec = Mod(1, p);
  qq = nomegen(mm, j1, j2, A, B, C, N, onec);
  nq1 = dmax + 3;
  for(h = 1, hmax,
    my(zz, base, mmax, zp, es, es1);
    zz   = qq^h;
    base = x^h/zz;
    mmax = (N-1)\h;
    zp = vector(mmax+1); zp[1] = onec + O(x^(N+1));
    for(m = 1, mmax, zp[m+1] = zp[m]*zz);
    es  = eqlist(h, dmax, N);
    es1 = vector(nq1, i, es[#es - nq1 + i]);            /* tail block */
    glist = List();
    for(i = 1, #GBASE,
      my(b = onec*GBASE[i]^h);
      listput(glist, b); listput(glist, -b);
      listput(glist, 1/b); listput(glist, -1/b));
    for(gi = 1, #glist,
      my(gam = glist[gi], wv);
      wv = wcoeffs(zp, mmax, base, gam, N);
      if(#matker(vmat(wv, es1, dmax)) == 0, next);
      if(#matker(vmat(wv, es, dmax)) == 0, next);
      /* candidate: redo exactly over Q */
      /* recover the exact rational gam from the index gi */
      my(gq = (-1)^(gi+1) * GBASE[(gi-1)\4 + 1]^(if((gi-1)%4 < 2, h, -h)), cert);
      cert = certify(mm, j1, j2, A, B, C, h, gq, dmax, N);
      if(cert != 0, return([h, cert[1], gq, cert[2], cert[3]]))
    )
  );
  0;
}
