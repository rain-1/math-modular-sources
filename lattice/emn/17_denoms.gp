/* denominator type of the EMN moments a_n (input to tau) */
default(parisizemax, 8*10^9);
NN = 400;
amom(n) = sum(k = 0, n, binomial(n,k)/(2*k+1)) / (2^(n+1)*(n+1));
av = vector(NN+1); for(n = 0, NN, av[n+1] = amom(n));
dlcm(m) = my(r = 1); for(k = 1, m, r = lcm(r,k)); r;

print("=== v_2(a_n) ===");
{ my(bad = 0); for(n = 0, NN, if(valuation(av[n+1],2) != -1, bad++)); print("  v_2(a_n) = -1 for all n<=", NN, "?  exceptions: ", bad); }

print();
print("=== denominator type tests ===");
{
my(t1 = 0, t2 = 0, t3 = 0, t4 = 0, t5 = 0);
for(n = 1, NN,
  my(d = denominator(av[n+1]));
  if((2*dlcm(2*n+1)) % d != 0, t1++);
  if((2*(n+1)*dlcm(2*n+1)) % d != 0, t2++);
  if((2*dlcm(n+1)*dlcm(2*n+1)) % d != 0, t3++);
  if((2*dlcm(2*n)) % d != 0, t4++);
  if((2*dlcm(n)) % d != 0, t5++);
);
print("  den(a_n) | 2*[1..2n+1]              : failures ", t1, " / ", NN);
print("  den(a_n) | 2*(n+1)*[1..2n+1]        : failures ", t2, " / ", NN);
print("  den(a_n) | 2*[1..n+1]*[1..2n+1]     : failures ", t3, " / ", NN);
print("  den(a_n) | 2*[1..2n]                : failures ", t4, " / ", NN);
print("  den(a_n) | 2*[1..n]                 : failures ", t5, " / ", NN);
}
print();
print("=== is den(a_n) essentially 2*[1..2n+1]^odd ? per-prime check ===");
{
for(nn = [50, 120, 300], 0);
foreach([50,120,300], n,
  my(d = denominator(av[n+1]), f = factor(d), bad = 0, mx = 0);
  for(i = 1, matsize(f)[1],
    my(p = f[i,1], e = f[i,2]);
    if(p == 2, next);
    my(pred = floor(log(2*n+1)/log(p)));
    if(e > pred, bad++);
    if(e > mx, mx = e);
  );
  print("  n=", n, ": distinct odd primes ", matsize(f)[1]-1, ", max exponent ", mx,
        ", primes exceeding floor(log_p(2n+1)): ", bad,
        ", largest prime ", f[matsize(f)[1],1], " vs 2n+1=", 2*n+1);
);
}
print();
print("=== rates ===");
print("  log den(a_n)/n and log(2*[1..2n+1])/n:");
foreach([50,100,200,300,400], n, print("   n=", n, "   log den/n = ", log(1.0*denominator(av[n+1]))/n, "   log(2 d_{2n+1})/n = ", log(2.0*dlcm(2*n+1))/n));
print();
print("=== the sigma-invariant: den(a_n) as a CDT type [1..b n] ===");
print("  b = lim log(den)/n / 1  (since log d_{bn} ~ bn):");
foreach([100,200,300,400], n, print("   n=", n, "   b_est = ", log(1.0*denominator(av[n+1]))/n));
print();
print("=== numerator side: is 2 d_{2n+1} a_n primitive? (gcd with the numerator) ===");
foreach([20,50,100], n, my(c = 2*dlcm(2*n+1)*av[n+1]); print("   n=", n, "   2 d_{2n+1} a_n in Z? ", denominator(c)==1, "   value bits ", #binary(numerator(c))));
quit;
