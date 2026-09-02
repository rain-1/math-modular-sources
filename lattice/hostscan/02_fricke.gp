/* 02_fricke.gp -- EXHAUSTIVE sweep of Fricke hosts in the Theorem 3.4 shape
     u = prod_d eta(d tau)^{r_d},  r_{N/d} = -r_d  (=> sum r_d = 0, u|W_N = 1/(C u)),
     ord_inf(u) = 1,   C^{-2} = prod_d (N/d)^{r_d},   C in Z_{>0},
     x = u/(1 + B u + C u^2),  F = D log u in M_2(Gamma_0(N)),  a_n = [x^n] F,
     lambda_{1,2} = B +- 2 sqrt(C),  c = lam1 lam2 = B^2 - 4C.
   Enumerated by the ANTISYMMETRIC cusp-order vector with deg(u) <= DEGMAX,
   using the Ligozat matrix  A[c,d] = N/(24 gcd(c,N/c) c) * gcd(c,d)^2/d.
   Output lines:  N | ord-vector | r-vector | C | deg | (then B-scan done in python)  */
default(parisizemax, 4000000000);
DEGMAX = 4;
NMAX   = 120;

{ligozat(n) = my(dv = divisors(n), k = #dv);
  matrix(k, k, i, j, my(c = dv[i], d = dv[j]);
    n/(24*gcd(c, n/c)*c) * gcd(c,d)^2/d );}

{epsc(n,c) = eulerphi(gcd(c, n/c));}

print("N|ord|r|C|deg");
{for(nn=2, NMAX,
  my(dv = divisors(nn), k = #dv);
  my(A = ligozat(nn));
  if(matrank(A) < k, next);
  my(Ai = A^(-1));
  /* pairs (c, nn/c) with c < nn/c ; index into dv */
  my(pairs = List());
  for(i=1,k, my(c=dv[i]); if(c*c < nn, my(j = 0);
      for(t=1,k, if(dv[t]==nn/c, j=t; break)); listput(pairs, [i,j])));
  my(np = #pairs);
  /* pair 1 is always (1, nn) since dv[1]=1 ; its order is fixed: ord_1=-1, ord_nn=+1 */
  /* budget: sum over pairs eps_c * |ord| <= DEGMAX, with the (1,nn) pair using 1 */
  my(budget = DEGMAX - 1);
  /* enumerate ord values on pairs 2..np, each |o| <= budget/eps */
  my(lim = vector(np, t, if(t==1, 0, floor(budget/max(1,epsc(nn, dv[pairs[t][1]]))))));
  my(tot = 1); for(t=2,np, tot *= (2*lim[t]+1));
  if(tot > 3000000, print("SKIP ",nn," tot=",tot); next);
  for(code = 0, tot-1,
    my(cc = code, ord = vector(k), ok = 1, used = 1);
    ord[1] = -1;  /* c = 1  (cusp 0) */
    my(inn = 0); for(t=1,k, if(dv[t]==nn, inn=t; break));
    ord[inn] = 1;
    for(t=2, np,
      my(m = 2*lim[t]+1, o = (cc % m) - lim[t]); cc = cc \ m;
      my(i = pairs[t][1], j = pairs[t][2]);
      ord[i] = o; ord[j] = -o;
      used += epsc(nn, dv[i]) * abs(o);
    );
    if(used > DEGMAX, next);
    my(r = Ai * ord~);
    for(t=1,k, if(denominator(r[t]) != 1, ok = 0; break));
    if(!ok, next);
    /* C^{-2} = prod_d (N/d)^{r_d} */
    my(Cm2 = prod(t=1,k, (nn/dv[t])^r[t]));
    my(C2 = 1/Cm2);
    if(denominator(C2) != 1, next);
    if(!issquare(C2), next);
    my(C = sqrtint(C2));
    if(C < 1, next);
    print(nn,"|",ord,"|",Vec(r~),"|",C,"|",used);
  );
);}
quit;
