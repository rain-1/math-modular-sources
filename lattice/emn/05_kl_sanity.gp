/* Independent sanity check of the Kubota-Leopoldt implementation used in 04_padic3.gp:
   verify L_3(1-n, chi_12) = -(1 - (chi teich^{-n})(3) 3^{n-1}) B_{n, chi teich^{-n}} / n. */
default(parisizemax, 4*10^9);
PREC = 60;
chr(a) = my(r = a % 12); if(r==1 || r==11, 1, if(r==5 || r==7, -1, 0));
chm4(a) = my(r = a % 4); if(r==1, 1, if(r==3, -1, 0));
teich(a) = if(a % 3 == 1, 1, if(a % 3 == 2, -1, 0));

lkl(sv, FF, JMAX) = my(tot = 0); for(a = 1, FF, if(a % 3 != 0 && chr(a) != 0, my(br = a/teich(a)); my(inner = sum(j = 0, JMAX, binomial(1 - sv, j) * (FF/a)^j * bernfrac(j))); tot += chr(a) * (br + O(3^PREC))^(1 - sv) * inner)); tot / (FF * (sv - 1));

/* generalized Bernoulli B_{n,chi} = f^{n-1} sum_{a=1}^{f} chi(a) B_n(a/f) */
bchi(n, ff, which) = my(t = 0); for(a = 1, ff, my(c = if(which == 4, chm4(a), chr(a))); if(c != 0, t += c * subst(bernpol(n), 'x, a/ff))); ff^(n-1) * t;

print("=== KL sanity: L_3(1-n, chi_12) vs generalized Bernoulli ===");
for(n = 1, 6, \
  my(lhs = lkl(1 - n, 12, 100)); \
  my(rhs); \
  if(n % 2 == 1, rhs = -(1 - chm4(3)*3^(n-1)) * bchi(n, 4, 4) / n, \
                 rhs = -(1) * bchi(n, 12, 12) / n); \
  print("  n=", n, ": L_3(", 1-n, ") = ", lift(lhs + O(3^12)), "  predicted = ", lift(rhs + O(3^12)), "  v_3(diff) = ", valuation(lhs - rhs, 3), "  [rhs exact = ", rhs, "]") \
);
quit;
