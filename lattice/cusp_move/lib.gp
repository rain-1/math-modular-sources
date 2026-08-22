/* lib.gp -- cusp-move library for second-order Apery-like rows.
   Claude (Opus 5), 2026-08-22.  PARI/GP 2.15.4.
   Row convention:  (n+1)^2 u_{n+1} = Pp(n) u_n - Qq(n) u_{n-1},  u_0=1, u_{-1}=0.
   Rows are stored as [Pp, Qq, name] with Pp,Qq polynomials in the variable nv.  */

nv = 'nv;

/* ---------- basic row data ---------- */
rowa(rw) = polcoef(rw[1], 2, nv);            \\ leading coeff of P  = a
rowd(rw) = polcoef(rw[2], 2, nv);            \\ leading coeff of Q  = d
rowb(rw) = subst(rw[1], nv, 0);              \\ u_1 = P(0) = b
/* characteristic roots: roots of x^2 - a x + d */
charpol_row(rw) = 'x^2 - rowa(rw)*'x + rowd(rw);
charroots(rw) = polroots(charpol_row(rw));
/* exponents at infinity are 1-r_i, r_i = roots of Q(n) */
qroots(rw) = if(poldegree(rw[2],nv)==2, polroots(subst(rw[2],nv,'x)), [0,0]~);

/* ---------- sequences ---------- */
/* A-row: u_0=1, u_1=P(0) */
seqA(rw, nmax) =
{ my(v = vector(nmax+1), pp = rw[1], qq = rw[2]);
  v[1] = 1;
  for(n = 0, nmax-1,
    my(pn = subst(pp,nv,n), qn = subst(qq,nv,n),
       prev = if(n==0, 0, v[n]));
    v[n+2] = (pn*v[n+1] - qn*prev)/(n+1)^2);
  v; }

/* companion: b_0=0, b_1=1 ; equivalently L B = t */
seqB(rw, nmax) =
{ my(v = vector(nmax+1), pp = rw[1], qq = rw[2]);
  v[1] = 0;
  for(n = 0, nmax-1,
    my(pn = subst(pp,nv,n), qn = subst(qq,nv,n),
       prev = if(n==0, 0, v[n]), inh = if(n==0, 1, 0));
    v[n+2] = (pn*v[n+1] - qn*prev + inh)/(n+1)^2);
  v; }

/* generalised companion: L Bh = g(t) = sum_{k>=1} gc[k] t^k   (gc a closure k->coef) */
seqBg(rw, nmax, gc) =
{ my(v = vector(nmax+1), pp = rw[1], qq = rw[2]);
  v[1] = 0;
  for(n = 0, nmax-1,
    my(pn = subst(pp,nv,n), qn = subst(qq,nv,n),
       prev = if(n==0, 0, v[n]));
    v[n+2] = (pn*v[n+1] - qn*prev + gc(n+1))/(n+1)^2);
  v; }

/* ---------- the cusp move on coefficient sequences ----------
   A^#(s) = (1+lam s)^{-al} A(t),  t = s/(1+lam s)
   a^#_n  = sum_m a_m (-lam)^{n-m} binomial(n+al-1, n-m)                  */
cmove(v, lam, al) =
{ my(nmax = #v-1, w = vector(nmax+1));
  for(n = 0, nmax,
    my(acc = 0);
    for(m = 0, n, acc += v[m+1]*(-lam)^(n-m)*binomial(n+al-1, n-m));
    w[n+1] = acc);
  w; }

/* inverse move: lam -> -lam, al -> -al?  Not quite; the inverse of
   (s = t/(1-lam t), z = (1-lam t)^al y) is (t = s/(1+lam s), y = (1+lam s)^{-al} z),
   i.e. cmove with lam -> -lam and al -> -al is NOT it.  Direct: use cmove(w,-lam,-al)?
   Check: applying cmove twice with (lam,al) then (-lam,-al):
   z(s)=(1-lam t)^al y ; then (1+lam s)^{-(-al)}... we verify numerically in 00_selftest.gp */
cmoveinv(v, lam, al) = cmove(v, -lam, -al);

/* ---------- recurrence fitting ----------
   look for Pp,Qq of degree <= dg with (n+1)^2 v_{n+1} = Pp(n) v_n - Qq(n) v_{n-1}.
   Returns [Pp,Qq] normalised, or 0.  */
fitrow(v, dg = 2) =
{ my(nun = 2*(dg+1), n0 = 1, ne = nun + 6, mt, rh, sol);
  if(#v - 2 < n0 + ne - 1, return(0));
  mt = matrix(ne, nun, i, j,
        my(n = n0 + i - 1);
        if(j <= dg+1, n^(j-1)*v[n+1], -n^(j-dg-2)*v[n]));
  rh = vectorv(ne, i, my(n = n0 + i - 1); (n+1)^2*v[n+2]);
  if(matrank(mt[1..nun,]) < nun, return(0));
  sol = matsolve(mt[1..nun,], rh[1..nun]);
  for(i = nun+1, ne, if(mt[i,]*sol != rh[i], return(0)));
  my(pp = sum(j=0,dg, sol[j+1]*nv^j), qq = sum(j=0,dg, sol[dg+2+j]*nv^j));
  [pp, qq]; }

/* full check of a fitted row against the sequence */
checkrow(rw, v) =
{ my(ok = 1);
  for(n = 1, #v-2,
    if((n+1)^2*v[n+2] != subst(rw[1],nv,n)*v[n+1] - subst(rw[2],nv,n)*v[n],
       ok = 0; break));
  ok; }

/* ---------- integrality ---------- */
/* minimal c>0 (a positive rational, we only search integers and 1/integers) with
   c^n v_n in Z for all n <= nmax.  We return the minimal positive integer c that works,
   or 0 if none up to bound. */
minscale(v, bound = 4096) =
{ my(nmax = #v-1, den = 1);
  for(n = 0, nmax, den = lcm(den, denominator(v[n+1])));
  if(den == 1, return(1));
  my(f = factor(den), c = 1);
  for(i = 1, #f~,
    my(p = f[i,1], e = 0, need = 0);
    \\ smallest e with p^(e n) killing denominator of v_n for all n
    for(n = 1, nmax,
      my(w = valuation(denominator(v[n+1]), p));
      if(w > 0, need = max(need, ceil(w/n))));
    c *= p^need);
  if(c > bound, return(0));
  \\ verify
  for(n = 0, nmax, if(denominator(c^n*v[n+1]) != 1, return(0)));
  c; }

rescale(v, c) = vector(#v, i, c^(i-1)*v[i]);
rescalerow(rw, c) = [c*subst(rw[1],nv,nv), c^2*rw[2]];

/* ---------- denominator exponent k ---------- */
dnvec(nmax) = { my(v = vector(nmax+1)); v[1]=1; for(n=1,nmax, v[n+1]=lcm(v[n],n)); v; }

/* smallest k>=0 with d_n^k b_n in Z for n<=nmax ; returns [k, sharp] */
denomexp(bv, nmax, kmax = 6) =
{ my(dn = dnvec(nmax), k = -1);
  for(kk = 0, kmax,
    my(ok = 1);
    for(n = 0, nmax, if(denominator(dn[n+1]^kk * bv[n+1]) != 1, ok = 0; break));
    if(ok, k = kk; break));
  k; }

/* ---------- score ---------- */
scorerow(rw, k) =
{ my(rt = charroots(rw), l1, l2, m1, m2);
  m1 = abs(rt[1]); m2 = abs(rt[2]);
  if(m1 < m2, my(tmp=m1); m1=m2; m2=tmp);
  [m1, m2, log(1/m2) - k]; }

/* ---------- p-adic Apery limit diagnostics ---------- */
/* v_p(b_n/a_n - b_{n-1}/a_{n-1}) */
padicdiff(av, bv, p, n) = valuation(bv[n+1]/av[n+1] - bv[n]/av[n], p);
