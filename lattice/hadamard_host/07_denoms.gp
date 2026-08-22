\\ 07_denoms.gp -- exact denominator types, v_p profiles (p<=50), 2-adic slopes,
\\ archimedean growth, for every function on the Hadamard host.
default(parisizemax, 12000000000);
\p 6000
outdir = "/home/ubuntu/code/math-modular-sources/lattice/hadamard_host/out/";
read(concat(outdir, "rows_raw.gp"));
read(concat(outdir, "rows_full.gp"));
read(concat(outdir, "had.gp"));
N = #HZZ;

\\ conditional function: F = b*(A_Z (*) B_N) - a*(A_Z (*) A_N), a/b a test rational for G
AB = bestappr(Catalan, 10^8); aa = numerator(AB); bb = denominator(AB);
print("test rational a/b = ", aa, "/", bb, "   |G - a/b| = ", 1.0*abs(Catalan - AB));
CND = vector(N, n, bb*HZB[n] - aa*HZZ[n]);
\\ pure host-side comparison: the plain conditional row form on the Nesterenko side
write(concat(outdir,"cond.gp"), "CND = ", CND, ";");

NAMES = ["A_Z(*)A_N","A_Z(*)B_N","A_N(*)B_Z","B_Z(*)B_N","W","COND"];
SEQS  = [HZZ, HZB, HNB, HBB, WW, CND];

DT = vector(6*N+1); DT[1] = 1;
for(k = 1, 6*N, DT[k+1] = lcm(DT[k], k));
Dn(n) = DT[6*n+1];

print("");
print("### 1. archimedean growth log|c_n|/n and 2-adic slope v_2(c_n)/n ###");
print("           n=100            n=200            n=300            n=400");
{ for(s = 1, 6,
    printf("%-10s", NAMES[s]);
    for(k = 1, 4, my(n = 100*k, c = SEQS[s][n]);
      printf("  %8.5f/%7.4f", log(abs(c))/n,
        (valuation(numerator(c),2)-valuation(denominator(c),2))/n));
    print("")); }

print("");
print("### 2. denominator structure: den(c_n) = 2^A * (odd part) ###");
print("  columns: n | v2(den) | v2(den)/n | largest odd prime P in den | P/n | max_{p>3n} v_p(den) | v_p exponent profile vs D_{6n}^k");
{ for(s = 1, 6,
    print("  -- ", NAMES[s], " --");
    for(k = 1, 4, my(n = 100*k, c = SEQS[s][n], dd = denominator(c), v2 = valuation(dd,2),
                     od, fa, big = 0, kmax = 0, kk);
      od = dd >> v2;
      fa = factor(od);
      if(matsize(fa)[1] > 0,
        big = fa[matsize(fa)[1],1];
        for(t = 1, matsize(fa)[1],
          if(fa[t,1] > 3*n, kk = fa[t,2]; if(kk > kmax, kmax = kk))));
      printf("   n=%3d  v2=%7d (%7.4f/n)  Pmax=%6d (%.4f n)  k_top=%d\n",
        n, v2, 1.0*v2/n, big, 1.0*big/n, kmax))); }

print("");
print("### 3. sharp LCM exponent: smallest k with D_{6n}^k * c_n in Z[1/2] ###");
{ for(s = 1, 6,
    my(row = "");
    for(k = 1, 8, my(n = 50*k, c = SEQS[s][n], dd = denominator(c), od, kk = 0);
      od = dd >> valuation(dd,2);
      while(od % 1 == 0 && (Dn(n)^kk) % od != 0 && kk < 12, kk++);
      row = concat(row, concat(Str(kk), " ")));
    print("  ", NAMES[s], ":  k(n) for n=50,100,...,400 : ", row)); }

print("");
print("### 4. v_p profile for p <= 50 at n = 200, 400 ###");
{ my(PL = primes(15));
  for(s = 1, 6,
    print("  -- ", NAMES[s], " --");
    for(k = 2, 2, my(n = 200*k, c = SEQS[s][n], dd = denominator(c), nu = numerator(c), row = "");
      for(t = 1, #PL, my(p = PL[t], e = valuation(dd,p) - valuation(nu,p),
                          fl = logint(6*n, p));
        row = concat(row, Str(p)); row = concat(row, ":");
        row = concat(row, Str(e)); row = concat(row, "/");
        row = concat(row, Str(fl)); row = concat(row, "  "));
      print("   n=", n, "  (v_p(c_n) as -e; e/floor(log_p 6n)) : ", row))); }
\q
