/* ap_lib.gp -- 5-adic (Z/5^MM) computation of the Calegari coefficients
   a_n, b_n, e_n for (p,k)=(5,1) on X_0(5).

     x(q)    = q * prod_{n>=1} ((1-q^{5n})/(1-q^n))^6         in q + q^2 Z[[q]]
     Estar   = E_2(tau) - 5 E_2(5 tau),  E_2 = 1-24 sum sigma_1(n) q^n
               (UNDIVIDED: constant term -4)
     Eprime  = sum_{n>=1} ( sum_{d|n, 5 !| d} d^{-3} ) q^n
     H_eta   = Estar*(Eprime + eta) = sum_n (a_n + eta b_n) x^n
   so   b_n = [x^n] Estar,  e_n = [x^n] Eprime,  a_n = sum_i b_i e_{n-i}.

   Everything lies in Z_5, so we work in Z/5^MM.  All ring operations only
   (no division by 5): the residues are EXACT reductions of the true values.

   NB: builtin names psi, M, Phi, S, cmp are avoided.
*/

/* prod_{n>=1} (1-q^n) by the pentagonal number theorem (exact, sparse) */
{
eulerP(nn) =
  my(ss, k, t1, t2);
  ss = 1 + O(q^(nn+1));
  k = 1;
  while(1,
    t1 = k*(3*k-1)/2; t2 = k*(3*k+1)/2;
    if(t1 > nn && t2 > nn, break);
    if(t1 <= nn, ss += (-1)^k*q^t1);
    if(t2 <= nn, ss += (-1)^k*q^t2);
    k++;
  );
  ss;
}

/* "peel" F(q) into the x-coordinate:  F = sum_n res[n+1] * x^n,  x = q*wq.
   winv = 1/wq.  nn = highest index wanted.  F must be known to O(q^(nn+1)). */
{
peel(F, winv, nn) =
  my(res, ff, cc);
  res = vector(nn+1);
  ff = F;
  for(n = 0, nn,
    cc = polcoeff(ff, 0);
    res[n+1] = cc;
    if(n < nn,
      ff = ((ff - cc)*(winv + O(q^(nn+1-n))))/q;
    );
  );
  res;
}

/* build everything.  returns [avec, bvec, evec], index i <-> n=i-1 */
{
buildall(nn, mm) =
  my(mo, one, P1, P5, wq, winv, E2, Es, ep, cf, iv, bb, ee, aa, hb, he, ha);
  mo = 5^mm; one = Mod(1, mo);
  P1 = eulerP(nn+1);
  P5 = subst(P1, q, q^5) + O(q^(nn+2));
  wq = one*(P5^6/P1^6) + O(q^(nn+2));
  winv = 1/wq;
  print("  [built wq, winv]  t=", gettime());
  E2 = 1 - 24*sum(n=1, nn+1, sigma(n)*q^n) + O(q^(nn+2));
  Es = one*(E2 - 5*subst(E2, q, q^5)) + O(q^(nn+2));
  /* Eprime */
  cf = vector(nn+1, i, Mod(0, mo));   /* cf[i] = coeff of q^i, i=1..nn ; slot nn+1 unused */
  for(d = 1, nn,
    if(d % 5,
      iv = Mod(d, mo)^(-3);
      forstep(n = d, nn, d, cf[n] += iv);
    );
  );
  ep = sum(n = 1, nn, cf[n]*q^n) + O(q^(nn+2));
  print("  [built Es, Eprime]  t=", gettime());
  bb = peel(Es, winv, nn);
  print("  [peeled Estar]  t=", gettime());
  ee = peel(ep, winv, nn);
  print("  [peeled Eprime]  t=", gettime());
  /* a = b*e as series in x */
  hb = sum(n = 0, nn, bb[n+1]*x^n) + O(x^(nn+1));
  he = sum(n = 0, nn, ee[n+1]*x^n) + O(x^(nn+1));
  ha = hb*he;
  aa = vector(nn+1, i, polcoeff(ha, i-1));
  print("  [convolved a]  t=", gettime());
  [aa, bb, ee];
}

/* 5-adic valuation of a t_INTMOD residue mod 5^mm (returns mm if 0) */
{
v5(z, mm) = my(u); u = lift(z); if(u == 0, mm, valuation(u, 5));
}

/* same as buildall but over an arbitrary modulus mo (e.g. a big prime).
   valid because a_n,b_n,e_n only need 1/d^3 for 5 !| d ; mo must be prime to
   all d<=nn with 5!|d (true for mo = 2^61-1). */
{
buildmod(nn, mo) =
  my(one, P1, P5, wq, winv, E2, Es, ep, cf, iv, bb, ee, aa, hb, he, ha);
  one = Mod(1, mo);
  P1 = eulerP(nn+1);
  P5 = subst(P1, q, q^5) + O(q^(nn+2));
  wq = one*(P5^6/P1^6) + O(q^(nn+2));
  winv = 1/wq;
  E2 = 1 - 24*sum(n=1, nn+1, sigma(n)*q^n) + O(q^(nn+2));
  Es = one*(E2 - 5*subst(E2, q, q^5)) + O(q^(nn+2));
  cf = vector(nn+1, i, Mod(0, mo));
  for(d = 1, nn,
    if(d % 5,
      iv = Mod(d, mo)^(-3);
      forstep(n = d, nn, d, cf[n] += iv);
    );
  );
  ep = sum(n = 1, nn, cf[n]*q^n) + O(q^(nn+2));
  bb = peel(Es, winv, nn);
  ee = peel(ep, winv, nn);
  hb = sum(n = 0, nn, bb[n+1]*x^n) + O(x^(nn+1));
  he = sum(n = 0, nn, ee[n+1]*x^n) + O(x^(nn+1));
  ha = hb*he;
  aa = vector(nn+1, i, polcoeff(ha, i-1));
  [aa, bb, ee];
}
