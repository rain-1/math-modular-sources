/* cdt_ab.gp : exact a_n,b_n for Calegari / CDT p-adic Apery constructions.
   Normalisation:
     x(q)   = q prod ((1-q^{p n})/(1-q^n))^{24/(p-1)}          (X_0(p), p-1|24)
     Estar  = E_{2k}(tau) - p^{2k-1} E_{2k}(p tau)   in Z[[q]], constant term 1-p^{2k-1}
             (nrm=1); with nrm=0 it is divided by (1-p^{2k-1}) to be 1+O(q)  [Calegari's
             normalisation, which for (p,k)=(2,1) is exactly multiplication by -1].
     Eprime = sum_{n>=1} ( sum_{d|n, p!|d} d^{-(2k+1)} ) q^n   (constant term omitted)
     H      = Estar*(Eprime + eta) = sum_n (a_n + eta*b_n) x^n
   so b_n = [x^n] Estar,  a_n = [x^n] (Estar*Eprime).
   NB: avoid builtin names psi,M,Phi,S,cmp.
*/

/* classical E_{2k} = 1 - (4k/B_{2k}) sum sigma_{2k-1}(n) q^n */
{
eis(k, nn) =
  my(c, ss);
  c = -4*k/bernfrac(2*k);
  ss = 1 + O(q^(nn+1));
  for(n = 1, nn, ss += c*sigma(n, 2*k-1)*q^n);
  ss;
}

/* q prod_{n>=1} ((1-q^{Nn})/(1-q^n))^s  ... eta-quotient hauptmodul, = q + ... */
{
hauptser(nn, N, s) =
  my(lg, tt);
  /* use log of the product: -s*sum_m ( sigma-ish ) ; simpler: direct product */
  lg = O(q^(nn+1));
  for(n = 1, nn,
    lg += s*log(1 - q^(N*n) + O(q^(nn+1))) - s*log(1 - q^n + O(q^(nn+1)));
  );
  tt = exp(lg);
  q*tt + O(q^(nn+2));
}

/* Eprime: sum_{n>=1} (sum_{d|n, p!|d} d^{-w}) q^n */
{
eprimeser(nn, p, w) =
  my(ss);
  ss = O(q^(nn+1));
  for(n = 1, nn, ss += sumdiv(n, d, if(d%p, 1/d^w, 0))*q^n);
  ss;
}

/* main driver for X_0(p) : returns [avec, bvec] indexed a[i]=a_{i-1} */
{
run0p(p, k, nn, nrm=1) =
  my(s, xq, qofx, ee, ep, ha, hb, aa, bb, mm);
  mm = nn+2;
  s = 24/(p-1);
  xq   = hauptser(mm, p, s);
  qofx = serreverse(xq);
  ee   = eis(k, mm) - p^(2*k-1)*subst(eis(k, mm), q, q^p + O(q^(mm+1)));
  if(nrm == 0, ee = ee/(1 - p^(2*k-1)));
  ep   = eprimeser(mm, p, 2*k+1);
  hb   = subst(ee, q, qofx);
  ha   = subst(ee*ep, q, qofx);
  aa = vector(nn+1, i, polcoeff(ha, i-1));
  bb = vector(nn+1, i, polcoeff(hb, i-1));
  [aa, bb];
}

/* ---------------- odd-character (Catalan-type) cases -------------------------
   chi = kronecker(D,.) odd quadratic, conductor |D|.
     G1  = 1 + c * sum_{n>=1} ( sum_{d|n} chi(d) ) q^n        (c = 2/L(0,chi))
     G'  = sum_{n>=1} ( sum_{d|n} chi(d) d^{-2} ) q^n         (constant term omitted)
   For D=-4 (p=2) : L(0,chi)=1/2, c=4, F_1 = (1/4)*G1 in Calegari's normalisation,
     and G' = Calegari's F'_{-1} = sum_{n>=0} q^{2n+1}(-1)^n/((2n+1)^2(1-q^{2n+1})).
   For D=-3 (p=3) : L(0,chi)=1/3, c=6.
   No p-stabilisation is needed: chi(d)=0 whenever p|d, so sum_{d|n} = sum_{d|n,p!|d}
   already.  (This is the exact analogue of Calegari's p=2 / chi_{-4} situation.)
*/
{
chieis1(nn, D, c) =
  my(ss);
  ss = 1 + O(q^(nn+1));
  for(n = 1, nn, ss += c*sumdiv(n, d, kronecker(D, d))*q^n);
  ss;
}
{
chieisprime(nn, D, w) =
  my(ss);
  ss = O(q^(nn+1));
  for(n = 1, nn, ss += sumdiv(n, d, kronecker(D, d)/d^w)*q^n);
  ss;
}

/* generic driver: level N hauptmodul q*prod((1-q^{Nn})/(1-q^n))^s, character D */
{
runchi(nn, N, s, D, c) =
  my(xq, qofx, g1, gp, ha, hb, mm);
  mm = nn+2;
  xq   = hauptser(mm, N, s);
  qofx = serreverse(xq);
  g1   = chieis1(mm, D, c);
  gp   = chieisprime(mm, D, 2);
  hb = subst(g1, q, qofx);
  ha = subst(g1*gp, q, qofx);
  [vector(nn+1, i, polcoeff(ha, i-1)), vector(nn+1, i, polcoeff(hb, i-1))];
}

run14(nn) = runchi(nn, 4, 8, -4, 4);   /* X_1(4)=X_0(4), chi_{-4}, Catalan */
run09(nn) = runchi(nn, 9, 3, -3, 6);   /* X_0(9), chi_{-3}                 */
run03(nn) = runchi(nn, 3, 12, -3, 6);  /* X_0(3), chi_{-3}  (comparison)   */
