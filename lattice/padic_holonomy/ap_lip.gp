/* ap_lip.gp -- TASK 3: archimedean Lipschitz data for x(q) on X_0(5).
   Lam(q) = | q d/dq log x(q) |.
   From x = q prod_{n>=1} ((1-q^{5n})/(1-q^n))^6 :
     q dlog x/dq = 1 + 6 sum_{n>=1} n ( q^n/(1-q^n) - 5 q^{5n}/(1-q^{5n}) )
                 = 1 + 6 sum_{m>=1} sigma_1(m) q^m - 30 sum_{m>=1} sigma_1(m) q^{5m}
                 = -(E_2(q) - 5 E_2(q^5))/4 .
   FLOATING POINT (\p 60).
*/
default(parisizemax, 12000000000);
default(realprecision, 60);

/* --- symbolic verification of the formula, exact, to O(q^60) --- */
{
verify() =
  my(nn, P1, P5, xq, lhs, rhs, rhs2, ok);
  nn = 60;
  P1 = 1 + O(q^(nn+1)); 
  for(n=1, nn, P1 = P1*(1 - q^n + O(q^(nn+1))));
  P5 = 1 + O(q^(nn+1));
  for(n=1, nn\5, P5 = P5*(1 - q^(5*n) + O(q^(nn+1))));
  xq = q*P5^6/P1^6 + O(q^(nn+1));
  lhs = q*deriv(xq)/xq;                              /* q x'/x */
  rhs = 1 + 6*sum(n=1, nn, n*( q^n/(1-q^n+O(q^(nn+1))) - 5*q^(5*n)/(1-q^(5*n)+O(q^(nn+1))) )) + O(q^(nn+1));
  rhs2 = -( (1 - 24*sum(m=1,nn,sigma(m)*q^m)) - 5*(1 - 24*sum(m=1,nn\5,sigma(m)*q^(5*m))) )/4 + O(q^(nn+1));
  ok = [lhs - rhs == 0, lhs - rhs2 == 0];
  print("formula verification (exact series to O(q^61)): [Lambert form, E2* form] = ", ok);
  print("  q dlogx/dq = ", lhs + O(q^9));
}
verify();

/* coefficients c_m of  q dlog x/dq = 1 + sum_{m>=1} c_m q^m   (exact integers) */
NMAX = 4000;
cc = vector(NMAX);
{ for(m = 1, NMAX, cc[m] = 6*sigma(m) - if(m%5==0, 30*sigma(m/5), 0)); }

/* evaluate on circle |q|=r at K points; return [maxmod, minmod, argmax, majorant] */
{
circ(r, kk) =
  my(dm, mx, mn, amx, val, z, tot, maj, lim);
  /* truncation: need c_m r^m < 1e-70 */
  lim = NMAX;
  while(lim > 10 && (lim^3)*r^lim < 1e-75, lim--);
  dm = vector(lim, m, cc[m]*r^m);
  mx = 0.; mn = 1e100; amx = 0.;
  for(k = 0, kk\2,
    z = exp(2*Pi*I*k/kk);
    tot = 1.;
    my(zp); zp = 1.;
    for(m = 1, lim, zp = zp*z; tot += dm[m]*zp);
    val = abs(tot);
    if(val > mx, mx = val; amx = 2*Pi*k/kk);
    if(val < mn, mn = val);
  );
  /* triangle-inequality majorant: 1 + 6 sum n (r^n/(1-r^n) + 5 r^{5n}/(1-r^{5n})) */
  maj = 1 + 6*sum(n=1, lim, n*( r^n/(1-r^n) + 5*r^(5*n)/(1-r^(5*n)) ));
  [mx, mn, amx, maj, lim];
}

rs = [0.60, 0.70, 0.75, 0.786, 0.80, 0.85, 0.90];
KK = 4000;
print("");
print("r        max|q dlogx/dq|        min|.|                 argmax   truncation  majorant(abs inside)");
{
for(i = 1, #rs,
  my(res);
  res = circ(rs[i], KK);
  print(rs[i], "  ", res[1], "  ", res[2], "  ", res[3], "  ", res[5], "  ", res[4]);
);
}
quit;
