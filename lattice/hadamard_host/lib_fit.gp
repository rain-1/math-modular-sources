\\ lib_fit.gp -- minimal-recurrence fitting over Q for integer/rational sequences,
\\ adapted from lattice/two_prime_holonomy/02_fit_recur.gp.
\\ Sequences are 1-based vectors: S[n] = u_n, n = 1..#S.
\\ Recurrence sought:  sum_{i=0}^{r} c_i(n) u_{n+i} = 0,  deg c_i <= d.

\\ coefficient matrix mod pp (pp = 0 -> over Q)
hh_mat(seqs, r, d, nlist, pp) =
{ my(uu = (d+1)*(r+1), nr = #nlist*#seqs, mm = matrix(nr, uu), ro = 0, val, u, n);
  for(s = 1, #seqs,
    u = seqs[s];
    for(k = 1, #nlist,
      n = nlist[k]; ro++;
      for(i = 0, r,
        val = u[n+i];
        if(pp, val = Mod(numerator(val), pp) / Mod(denominator(val), pp));
        for(j = 0, d, mm[ro, i*(d+1)+j+1] = val * n^j))));
  mm;
}

\\ sweep for the smallest (r,d) with a nontrivial kernel
hh_sweep(seqs, rmax, dmax, nlist, pp) =
{ my(res = List(), mm, kk, nu);
  for(r = 1, rmax,
    for(d = 0, dmax,
      if((r+1)*(d+1) > #nlist*#seqs - 6, print("   r=", r, ": out of data at d=", d); break);
      mm = hh_mat(seqs, r, d, nlist, pp); kk = matker(mm); nu = #kk;
      if(nu > 0,
        print("   r=", r, " d=", d, "  unknowns=", (r+1)*(d+1),
              "  rows=", #nlist*#seqs, "  nullity=", nu);
        listput(res, [r, d, nu]);
        break)));
  Vec(res);
}

\\ exact kernel by multimodular CRT (nullity 1 assumed)
hh_ker(seqs, r, d, nlist0, nprimes) =
{ my(uu = (d+1)*(r+1), p0 = 2^62, vv, kk, pos = 0, acc, modu = 1, cand, prev = 0,
     nl, nlist);
  nl = ceil(uu/#seqs) + 20;
  nlist = if(#nlist0 > nl, vector(nl, k, nlist0[k]), nlist0);
  p0 = nextprime(p0);
  kk = matker(hh_mat(seqs, r, d, nlist, p0));
  if(#kk != 1, error("hh_ker: nullity ", #kk, " != 1"));
  vv = kk[,1];
  for(t = 1, uu, if(lift(vv[t]) != 0, pos = t; break));
  acc = vector(uu);
  for(c = 1, nprimes,
    p0 = nextprime(p0+1);
    kk = matker(hh_mat(seqs, r, d, nlist, p0));
    if(#kk != 1, next);
    vv = kk[,1]; if(lift(vv[pos]) == 0, next);
    vv = vv / vv[pos];
    for(t = 1, uu, acc[t] = lift(chinese(Mod(acc[t], modu), vv[t])));
    modu = modu*p0;
    cand = vector(uu, t, bestappr(Mod(acc[t], modu)));
    if(type(cand[1]) != "t_INT" && type(cand[1]) != "t_FRAC", next);
    if(cand == prev, return([cand, modu, c]));
    prev = cand);
  print("   WARNING: reconstruction did not stabilise after ", nprimes, " primes");
  [prev, modu, nprimes];
}

hh_tocoefs(v, r, d) = vector(r+1, i, sum(j = 0, d, v[(i-1)*(d+1)+j+1] * 'nn^j));

\\ clear denominators / content
hh_clean(cs) =
{ my(r1 = #cs, den = 1, gg = 0, out);
  for(i = 1, r1, den = lcm(den, denominator(content(cs[i]))));
  out = vector(r1, i, cs[i]*den);
  for(i = 1, r1, gg = gcd(gg, content(out[i])));
  if(gg != 0, out = vector(r1, i, out[i]/gg));
  out;
}

\\ verify over Q on n = nlo..nhi
hh_verify(seqs, cs, r, nlo, nhi) =
{ my(bad = 0, tt);
  for(s = 1, #seqs,
    for(n = nlo, nhi,
      tt = sum(i = 0, r, subst(cs[i+1], 'nn, n) * seqs[s][n+i]);
      if(tt != 0, bad++; if(bad <= 3, print("   FAIL seq=", s, " n=", n)))));
  bad;
}

\\ singular structure of the generating function.
\\ For f(x) = sum u_n x^n and sum_i c_i(n) u_{n+i} = 0, the ODE in theta = x d/dx is
\\   sum_i x^(r-i) c_i(theta - i) f = 0 ; the symbol in theta^D is
\\   Pgf(x) = sum_{i : deg c_i = D} lead(c_i) x^(r-i).
hh_sing(cs, r) =
{ my(dg, dd, pgf = 0, lead, rts);
  dg = vector(r+1, i, poldegree(cs[i], 'nn));
  dd = vecmax(dg);
  lead = vector(r+1, i, polcoef(cs[i], dg[i], 'nn));
  for(i = 1, r+1, if(dg[i] == dd, pgf += lead[i] * 'x^(r-i+1)));
  print("   deg c_i     = ", dg);
  print("   lead c_i    = ", lead);
  print("   Pgf(x)      = ", pgf, "   (D = ", dd, ")");
  if(pgf != 0 && poldegree(pgf,'x) > 0,
    print("   Pgf factored= ", factor(pgf));
    rts = polroots(pgf);
    rts = select(z -> abs(z) > 1e-30, Vec(rts));
    print("   sing pts    = ", rts);
    print("   |sing|      = ", vector(#rts, t, abs(rts[t])));
    print("   min |sing|  = ", if(#rts, vecmin(vector(#rts,t,abs(rts[t]))), "-")),
    print("   Pgf constant: no finite nonzero singular point"));
  [dg, pgf];
}

\\ full driver: fit, verify, report.  returns [r, d, cs]
hh_fit(nm, seqs, rmax, dmax, nfit, nall) =
{ my(nlist, res, r = 0, d, ek, cs);
  print("");
  print("=== fit: ", nm, "   (#seqs=", #seqs, ", fit n<=", nfit, ", data n<=", nall, ") ===");
  nlist = vector(nfit - rmax, k, k);
  res = hh_sweep(seqs, rmax, dmax, nlist, nextprime(2^61));
  if(#res == 0, print("   NO recurrence in range r<=", rmax, " d<=", dmax); return(0));
  for(t = 1, #res, if(res[t][3] == 1, r = res[t][1]; d = res[t][2]; break));
  if(r == 0, r = res[1][1]; d = res[1][2]; print("   (no nullity-1 hit; using first)"));
  nlist = vector(nfit - r, k, k);
  ek = hh_ker(seqs, r, d, nlist, 60);
  cs = hh_clean(hh_tocoefs(ek[1], r, d));
  for(i = 0, r, print("   c_", i, "(n) = ", cs[i+1]));
  print("   verify on n=1..", nall-r, " (fit used n<=", nfit-r, "): failures = ",
        hh_verify(seqs, cs, r, 1, nall-r));
  hh_sing(cs, r);
  [r, d, cs];
}

\\ extend a sequence using a fitted recurrence (leading coefficient must not vanish)
hh_extend(S, cs, r, N) =
{ my(u = vector(N), lc);
  for(k = 1, #S, u[k] = S[k]);
  for(n = #S - r + 1, N - r,
    lc = subst(cs[r+1], 'nn, n);
    if(lc == 0, error("hh_extend: leading coefficient vanishes at n=", n));
    u[n+r] = -sum(i = 0, r-1, subst(cs[i+1], 'nn, n)*u[n+i])/lc);
  u;
}

\\ binary search on d for fixed r: smallest d with nontrivial kernel in [dlo,dhi]
hh_mind(seqs, r, dlo, dhi, nlist, pp) =
{ my(lo = dlo, hi = dhi, mid, nu, best = -1, bestnu = 0);
  \\ first check dhi
  if((r+1)*(dhi+1) > #nlist*#seqs - 4, print("   r=", r, ": dhi=", dhi, " exceeds data"); return([-1,0]));
  nu = #matker(hh_mat(seqs, r, dhi, nlist, pp));
  print("   r=", r, " d=", dhi, " nullity=", nu);
  if(nu == 0, return([-1,0]));
  best = dhi; bestnu = nu;
  while(lo < hi,
    mid = (lo+hi)\2;
    nu = #matker(hh_mat(seqs, r, mid, nlist, pp));
    print("   r=", r, " d=", mid, " nullity=", nu);
    if(nu > 0, hi = mid; best = mid; bestnu = nu, lo = mid+1));
  [best, bestnu];
}
