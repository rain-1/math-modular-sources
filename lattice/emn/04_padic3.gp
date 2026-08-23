/* (e) 3-adic claim:  sum a_n (3/2)^n  converges 3-adically; is it L_3(2,chi_12)?
   Independent Kubota-Leopoldt computation via Washington, Thm 5.11.
   Renamed: chi -> chr, omega -> teich, L -> lkl. */
default(parisizemax, 6*10^9);

PREC = 90;                                   /* 3-adic digits */
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1)) / (2^(n+1)*(n+1));

print("=== (e0) 3-adic convergence of a_n (3/2)^n ===");
{
for(n = 1, 12, print("  n=",n,"  v_3(a_n)=", valuation(amom(n),3), "  v_3(a_n (3/2)^n)=", valuation(amom(n)*(3/2)^n,3)));
my(mn = 10^6);
for(n = 1, 200, my(v = valuation(amom(n),3)); if(v < mn, mn = v));
print("  min v_3(a_n), n<=200: ", mn, "   (so v_3(a_n(3/2)^n) >= n - O(log n))");
}

print();
print("=== (e1) the EMN 3-adic sum ===");
{
my(NT = 260, s = 0);
for(n = 0, NT, s += amom(n)*(3/2)^n);
EMN3 = s + O(3^PREC);
print("  Phi_3(3/2) = sum_{n<=", NT, "} a_n (3/2)^n  (3-adically):");
print("  ", EMN3);
/* stability check */
my(s2 = 0); for(n = 0, 200, s2 += amom(n)*(3/2)^n);
print("  v_3(partial(260) - partial(200)) = ", valuation(s - s2, 3));
}

print();
print("=== (e2) Kubota-Leopoldt L_3(2, chi_12) via Washington Thm 5.11 ===");
/* chi_12 = chi_{-4} * chi_{-3}, the even quadratic character of conductor 12.
   chr(1)=1, chr(5)=-1, chr(7)=-1, chr(11)=1, else 0.  */
chr(a) = my(r = a % 12); if(r==1 || r==11, 1, if(r==5 || r==7, -1, 0));
/* Teichmueller at p=3: teich(a) = +1 if a=1 mod 3, -1 if a=2 mod 3 */
teich(a) = if(a % 3 == 1, 1, if(a % 3 == 2, -1, 0));

/* L_p(s,chi) = 1/(F(s-1)) * sum_{a=1,(a,p)=1}^F chr(a) <a>^{1-s} sum_{j>=0} binom(1-s,j) (F/a)^j B_j
   with F = lcm(f_chi, p), <a> = a/teich(a).  */
lkl(sv, FF, JMAX) = {
  my(tot = 0);
  for(a = 1, FF,
    if(a % 3 == 0, next);
    if(chr(a) == 0, next);
    my(br = a/teich(a));                     /* <a> in 1 + 3Z_3 */
    my(inner = 0);
    for(j = 0, JMAX,
      inner += binomial(1 - sv, j) * (FF/a)^j * bernfrac(j)
    );
    tot += chr(a) * (br + O(3^PREC))^(1 - sv) * inner
  );
  tot / (FF * (sv - 1));
};

L3chi12 = lkl(2, 12, 140);
print("  L_3(2, chi_12) = ", L3chi12);

print();
print("=== (e3) sanity check of the KL implementation at s = 1-n ===");
print("    Washington: L_p(1-n, chi) = -(1 - chi*teich^{-n}(p) p^{n-1}) B_{n, chi teich^{-n}} / n");
{
/* chi teich^{-n}: teich = chi_{-3} as a Dirichlet char mod 3.
   chi_12 * chi_{-3}^{-n} = chi_{-4} * chi_{-3}^{1-n}. For n odd -> chi_{-4}*chi_{-3}^{even}=chi_{-4} (cond 4)
   For n even -> chi_{-4}*chi_{-3} = chi_12 (cond 12). */
bernchi(n, ff, chf) = {  /* generalized Bernoulli B_{n,chi} */
  my(t = 0);
  for(a = 1, ff,
    if(chf(a) == 0, next);
    /* B_{n,chi} = f^{n-1} sum_{a=1}^{f} chi(a) B_n(a/f) */
    t += chf(a) * subst(bernpol(n), 'x, a/ff)
  );
  ff^(n-1) * t;
};
chm4(a) = my(r = a % 4); if(r==1, 1, if(r==3, -1, 0));
for(n = 1, 4,
  my(lhs = lkl(1 - n, 12, 140));
  my(rhs);
  if(n % 2 == 1,
    /* chi teich^{-n} = chi_{-4}, conductor 4; chi(3) = chm4(3) = -1 */
    rhs = -(1 - chm4(3)*3^(n-1)) * bernchi(n, 4, chm4) / n,
    /* n even: chi teich^{-n} = chi_12, conductor 12; chi_12(3)=0 */
    rhs = -(1 - 0) * bernchi(n, 12, chr) / n
  );
  print("  n=", n, ":  L_3(", 1-n, ",chi_12) = ", lhs);
  print("        predicted             = ", rhs + O(3^PREC));
  print("        v_3(difference)       = ", valuation(lhs - rhs, 3))
);
}

print();
print("=== (e4) COMPARISON ===");
print("  Phi_3(3/2)      = ", EMN3);
print("  L_3(2,chi_12)   = ", L3chi12);
print("  v_3(difference) = ", valuation(EMN3 - L3chi12, 3));
print("  ratio           = ", EMN3/L3chi12);
print();
print("  Also test simple rational multiples:");
foreach([1, 10/9, 9/10, 2, 1/2, 3/2, 2/3, -1, 4/3, 3/4, 5/3], c,
  my(d = EMN3 - c*L3chi12); print("    c=", c, "  v_3(Phi_3 - c*L_3) = ", valuation(d,3)));
quit;
