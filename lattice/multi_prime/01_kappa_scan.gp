\\ /home/ubuntu/code/math-modular-sources/lattice/multi_prime/01_kappa_scan.gp
\\ Step 1a: for every candidate DECAYER, find ALL primes carrying p-power
\\ denominators at linear rate (kappa_p > 0), by factoring den(Q_n) exactly.
\\ Run: timeout 3000 gp -q 01_kappa_scan.gp
default(parisizemax, 8000000000);
read("/home/ubuntu/code/math-modular-sources/lattice/multi_prime/lib.gp");

primesin(x) = if(x == 0, [], Vec(factor(denominator(x))[,1]~));

report(nm, qv, idxs) = {
  my(allp = [], hi = idxs[#idxs]);
  print("");
  print("### ", nm);
  for(k = 1, #idxs, allp = setunion(allp, Set(primesin(qv[idxs[k]+1]))));
  allp = vecsort(Vec(allp));
  print("  primes in den(Q) : ", allp);
  for(i = 1, #allp,
    my(p = allp[i], vr = vector(#idxs, k, valuation(denominator(qv[idxs[k]+1]), p)));
    print("    p=", p, "  v_p(den Q) at ", idxs, " = ", vr, "   kappa_p ~ ", 1.0*vr[#vr]/hi));
  for(i = 1, #allp,
    my(p = allp[i], vr = vector(#idxs, k, valuation(numerator(qv[idxs[k]+1]), p)));
    print("    p=", p, "  v_p(num Q) = ", vr));
  \\ small primes not in den: numerator content (negative kappa)
  for(i = 1, 6,
    my(p = prime(i), vr = vector(#idxs, k, valuation(qv[idxs[k]+1], p)));
    print("    p=", p, "  v_p(Q)     = ", vr, "   rate ", 1.0*vr[#vr]/hi));
};

zudpart() = {
  my(MZ = 200, zz, zq, zp, ok, dev);
  zz = zudrow(MZ); zq = zz[1]; zp = zz[2];
  report("Zudilin Catalan row Q_m", zq, [20,40,80,120,160,200]);
  print("  den(P_m) primes at m=40 : ", primesin(zp[41]));
  dev = vecsort(Set(vector(MZ, m, valuation(denominator(zq[m+1]),2) - (4*m - 2*hammingweight(m)))));
  print("  v_2(den Q_m) - (4m - 2 s_2(m)), m<=200, distinct values : ", dev);
  ok = sum(m = 1, MZ, if(denominator(zq[m+1]) >> valuation(denominator(zq[m+1]),2) != 1, 1, 0));
  print("  number of m<=200 with odd part of den(Q_m) not 1 : ", ok);
};

chi3part() = {
  my(alphas = [[2,1],[9,5],[7,4],[5,3],[3,2],[4,3],[5,4],[1,1],[3,4],[1,2]], KM = 18, idxs = [6,12,18]);
  for(j = 1, #alphas,
    my(pp = alphas[j][1], qq = alphas[j][2], qv = vector(KM+1));
    qv[1] = 27;
    for(k = 1, KM, qv[k+1] = chi3row(pp*k, qq*k)[1]);
    report(Str("chi_-3 family alpha = ", pp, "/", qq, " ; a=", pp, "k b=", qq, "k"), qv, idxs));
};

zudpart();
chi3part();
quit
