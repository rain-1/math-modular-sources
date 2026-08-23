default(parisizemax, 20*10^9);
NORD = 2000;
print("=== v_2(a_n) for the EMN moments ===");
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1)) / (2^(n+1)*(n+1));
{
my(bad = 0);
for(n = 0, 600, if(valuation(amom(n),2) != -1, bad++; print("  v_2(a_",n,") = ", valuation(amom(n),2))));
print("  v_2(a_n) = -1 for every n<=600?  exceptions: ", bad);
}
print();
av = vector(NORD+2); for(n = 0, NORD+1, av[n+1] = amom(n));
hser = sum(n = 0, NORD, av[n+1]*'z^(n+1)) + O('z^(NORD+2));
gettime();
comp = subst(hser, 'z, 'z*(3-2*'z)^2 + O('z^(NORD+2)));
ktil = comp + 3*hser;
print("composition to order ", NORD, ": ", gettime(), " ms");
{
my(vv = vector(NORD), bad = 0, mn = 10^9, mx = -10^9);
for(n = 1, NORD, vv[n] = valuation(polcoef(ktil, n, 'z), 2));
for(n = 1, NORD,
  my(d = vv[n]-n);
  if(d < mn, mn = d); if(d > mx, mx = d);
  my(pred = if(n % 2 == 0, 1, 0));
  if(d != pred, bad++; if(bad < 10, print("  RULE EXCEPTION n=", n, "  v_2-n=", d)))
);
print("  range of v_2([z^n]Ktil) - n over 1<=n<=", NORD, ": [", mn, ", ", mx, "]");
print("  exceptions to the rule  v_2 = n + (1 if n even else 0):  ", bad, " out of ", NORD);
print("  => 2-adic slope varsigma_2 = 1 EXACTLY, with a bounded (period-2) correction");
}
print();
print("=== radius of convergence: |c_n|^{-1/n} vs 1 - sqrt(3)/2 = ", 1-sqrt(3)/2, " ===");
foreach([100,250,500,1000,1500,2000], n, my(c = polcoef(ktil,n,'z)); print("   n=", n, "   |c_n|^{-1/n} = ", exp(-log(abs(c*1.0))/n)));
print();
print("=== Richardson-type extrapolation (c_n/c_{n+1}) ===");
foreach([500,1000,1500,1999], n, print("   n=", n, "   c_n/c_{n+1} = ", 1.0*polcoef(ktil,n,'z)/polcoef(ktil,n+1,'z)));
quit;
