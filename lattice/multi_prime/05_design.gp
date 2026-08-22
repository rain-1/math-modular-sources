\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/05_design.gp
\\ Master formula (paper/sections/05_two_row.tex eq. master) for the conductor-6
\\ decayer against the modular engines C, B, F, comparing
\\    (i) no bridge          G = 0
\\   (ii) single prime p=3   G = min(sig3eng r, sig3dec) log 3
\\  (iii) single prime p=2   G = min(sig2eng r, sig2dec) log 2
\\   (iv) both primes        G = sum of the two
\\ Decayer sampled at index b (gamma = 1), engine at r*b.
\\ K(r) = log lcm(d_(rb)^2, lcm(den Q_b, den P_b)) / b   measured exactly.
default(parisizemax, 8000000000);
default(realprecision, 4000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

Lc = (zetahurwitz(2,1/3) - zetahurwitz(2,2/3))/9;
dn(m) = if(m < 1, 1, lcm(vector(m, i, i)));

\\ engines: [name, k, sig2, sig3, log rho2, r_inf]
{ENG = [["C", 2, 0, 2, 0.0, 1/2], ["B", 2, 0, 3, 1.5*log(3), 1/2], ["F", 2, 3, 2, log(8), 5/8]];}

\\ measured decayer profile at index b
prof6(pp, qq, kk) = {
  my(bb = qq*kk, r = chi6row(pp*kk, qq*kk), qv = r[1], pv = r[2], den, dp);
  den = denominator(qv); dp = denominator(pv);
  [bb, log(abs(1.0*qv))/bb, log(abs(1.0*(qv*Lc - pv)))/bb, lcm(den,dp)];
};

\\ delta from the master formula
delt(kk, rr, lrho, llam, lLam, gg) = {
  my(ff = 0.5*(kk + rr*lrho + llam - gg), hh);
  hh = ff + lLam - llam;
  [ff, hh, 1 - ff/hh];
};

sig2dec = 7.95;   \\ measured, per unit b (uniform over the family)
sig3dec = 2.90;

scanmember(pp, qq, kk) = {
  my(pr = prof6(pp, qq, kk), bb = pr[1], lLam = pr[2], llam = pr[3], Den = pr[4]);
  print("");
  printf("== conductor-6 decayer  alpha = %d/%d  (b = %d):  logLam = %.4f  loglam = %.4f  eta = %.4f\n",
         pp, qq, bb, lLam, llam, log(1.0*Den)/bb);
  for(e = 1, #ENG,
    my(nm = ENG[e][1], keng = ENG[e][2], s2e = ENG[e][3], s3e = ENG[e][4], lrho = ENG[e][5],
       best = [0,0,0,0,0], bestlab = "");
    for(j = 1, 400,
      my(rr = j*0.02, kk2, g0, g2, g3, g23, res);
      kk2 = log(1.0*lcm(dn(floor(rr*bb))^keng, Den))/bb;
      g3 = min(s3e*rr, sig3dec)*log(3);
      g2 = min(s2e*rr, sig2dec)*log(2);
      g23 = g2 + g3;
      for(w = 1, 4,
        my(gg = [0.0, g3, g2, g23][w], lab = ["none","p=3","p=2","p=2&3"][w]);
        res = delt(kk2, rr, lrho, llam, lLam, gg);
        if(res[3] > best[w], best[w] = res[3]; if(w == 4, bestlab = Str("r=", rr)))));
    printf("   engine %s :  delta(no bridge) = %.4f   delta(p=3) = %.4f   delta(p=2) = %.4f   delta(BOTH) = %.4f   %s\n",
           nm, best[1], best[2], best[3], best[4], bestlab));
};

print("=== master-formula scan, conductor-6 decayer x modular engines ===");
print("decayer slopes used: sig2 = ", sig2dec, "  sig3 = ", sig3dec, "  (per unit of b)");
scanmember(2,1,60);
scanmember(3,2,40);
scanmember(5,3,32);
scanmember(4,3,32);
scanmember(5,4,28);
scanmember(1,1,40);
scanmember(2,3,32);
scanmember(1,2,32);
scanmember(1,3,24);
scanmember(1,4,20);
scanmember(1,6,16);
quit
