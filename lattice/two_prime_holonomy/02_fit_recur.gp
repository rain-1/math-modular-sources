default(parisizemax, 8000000000);
\\ 02_fit_recur.gp -- fit a JOINT polynomial-coefficient recurrence
\\   sum_{i=0}^{r} c_i(n) u_{n+i} = 0,   deg c_i <= d
\\ satisfied by BOTH (Q_b) and (P_b) of a well-poised family.
\\ Stage 1: sweep (r,d) modulo a large prime, record kernel dimension.
\\ Stage 2: at the minimal (r,d) with nullity 1, reconstruct the kernel exactly
\\          over Q by multi-modular CRT + rational reconstruction.
\\ Stage 3: verify the exact recurrence on ALL indices (including held-out ones).
\\ Stage 4: characteristic polynomial, its roots, generating-function singularities.

outdir = "/home/ubuntu/code/math-modular-sources/lattice/two_prime_holonomy/out/";

NFIT  = 0;   \\ set per run: highest index used in the FIT
NALL  = 0;   \\ highest index available (verification uses all of them)

\\ ---- build the coefficient matrix mod pp (or over Q if pp==0)
\\ rows: for each n in nlist and each sequence in seqs, entries
\\       col (i*(d+1)+j+1) = n^j * u_{n+i}
buildmat(seqs, r, d, nlist, pp) =
{ my(uu = (d+1)*(r+1), nr = #nlist*#seqs, mm = matrix(nr, uu), ro = 0, val);
  for(s = 1, #seqs,
    my(u = seqs[s]);
    for(k = 1, #nlist,
      my(n = nlist[k]); ro++;
      for(i = 0, r,
        val = u[n+i+1];
        if(pp, val = Mod(numerator(val), pp) / Mod(denominator(val), pp));
        for(j = 0, d, mm[ro, i*(d+1)+j+1] = val * n^j))));
  mm;
};

\\ ---- stage 1 sweep
sweep(seqs, rmax, dmax, nlist, pp) =
{ my(res = List());
  for(r = 1, rmax,
    for(d = 2, dmax,
      if((r+1)*(d+1) > #nlist*#seqs - 10, print("   r=", r, ": out of data at d=", d); break);
      my(mm = buildmat(seqs, r, d, nlist, pp), kk = matker(mm), nu = #kk);
      if(nu > 0,
        print("   r=", r, " d=", d, "  unknowns=", (r+1)*(d+1),
              "  rows=", #nlist*#seqs, "  nullity=", nu);
        listput(res, [r, d, nu]);
        break)));
  res;
};

\\ ---- exact kernel vector by multi-modular reconstruction (assumes nullity 1)
exactker(seqs, r, d, nlist0, nprimes) =
{ my(uu = (d+1)*(r+1), p0 = 2^62, vv, kk, pos = 0, acc, modu = 1, cand, prev = 0,
     nl = ceil(uu/#seqs) + 25, nlist);
  nlist = if(#nlist0 > nl, vector(nl, k, nlist0[k]), nlist0);
  \\ first pass to locate a nonzero coordinate
  p0 = nextprime(p0);
  kk = matker(buildmat(seqs, r, d, nlist, p0));
  if(#kk != 1, error("exactker: nullity ", #kk, " != 1"));
  vv = kk[,1];
  for(t = 1, uu, if(lift(vv[t]) != 0, pos = t; break));
  acc = vector(uu);
  for(c = 1, nprimes,
    p0 = nextprime(p0+1);
    kk = matker(buildmat(seqs, r, d, nlist, p0));
    if(#kk != 1, next);
    vv = kk[,1]; if(lift(vv[pos]) == 0, next);
    vv = vv / vv[pos];
    my(nm = modu*p0);
    for(t = 1, uu, acc[t] = lift(chinese(Mod(acc[t], modu), vv[t])));
    modu = nm;
    cand = vector(uu, t, bestappr(Mod(acc[t], modu)));
    if(type(cand[1]) != "t_INT" && type(cand[1]) != "t_FRAC", next);
    if(cand == prev, return([cand, modu, c]));
    prev = cand);
  print("   WARNING: reconstruction not stabilised after ", nprimes, " primes");
  [prev, modu, nprimes];
};

\\ ---- turn a flat coefficient vector into the polynomials c_i(n)
tocoefs(v, r, d) = vector(r+1, i, sum(j = 0, d, v[(i-1)*(d+1)+j+1] * 'nn^j));

\\ ---- exact verification over Q
verify(seqs, cs, r, nlo, nhi) =
{ my(bad = 0);
  for(s = 1, #seqs,
    my(u = seqs[s]);
    for(n = nlo, nhi,
      my(tt = sum(i = 0, r, subst(cs[i+1], 'nn, n) * u[n+i+1]));
      if(tt != 0, bad++; if(bad <= 3, print("   FAIL seq=", s, " n=", n)))));
  bad;
};

\\ ---- driver
run(nm, rmax, dmax, sel) =
{ my(fn = concat(concat(outdir, nm), ".seq"), seqs, nlist, res, r, d, ek, cs,
     lead, chi, rts, sing, den, gg, tag);
  read(fn);
  if(sel == 0, seqs = [QV, PV]; tag = "JOINT (Q,P)",
     sel == 1, seqs = [QV];     tag = "Q only",
                seqs = [PV];    tag = "P only");
  NALL = #QV - 1;
  NFIT = NALL - 8;                      \\ hold out the last 8 indices
  print("");
  print("=== ", nm, " / ", tag, "  (NALL=", NALL, ", fit on n<=", NFIT,
        ", hold out ", NALL-NFIT, ") ===");
  print("   Q_0..Q_3 = ", vector(4, t, QV[t]));
  print("   P_0..P_3 = ", vector(4, t, PV[t]));
  nlist = vector(NFIT - rmax, k, k);    \\ n = 1 .. NFIT-rmax
  print("  -- stage 1: mod-p sweep, r=1..", rmax, " d=2..", dmax);
  res = sweep(seqs, rmax, dmax, nlist, nextprime(2^61));
  if(#res == 0, print("   NO recurrence found in the sweep range"); return(0));
  \\ pick the entry with nullity 1 and smallest r
  r = 0;
  for(t = 1, #res, if(res[t][3] == 1, r = res[t][1]; d = res[t][2]; break));
  if(r == 0, print("   no nullity-1 fit found; using the first hit"); r = res[1][1]; d = res[1][2]);
  print("  -- stage 2: exact reconstruction at r=", r, " d=", d);
  \\ trim nlist to this r
  nlist = vector(NFIT - r, k, k);
  ek = exactker(seqs, r, d, nlist, 40);
  print("   primes used: ", ek[3]);
  cs = tocoefs(ek[1], r, d);
  den = lcm(vector(#ek[1], t, denominator(ek[1][t])));
  cs = vector(r+1, i, cs[i]*den);
  gg = 0; for(i = 1, r+1, gg = gcd(gg, content(cs[i])));
  if(gg != 0, cs = vector(r+1, i, cs[i]/gg));
  for(i = 0, r, print("   c_", i, "(n) = ", cs[i+1]));
  for(i = 0, r, print("   factor c_", i, " = ", factor(cs[i+1])));
  print("  -- stage 3: exact verification over Q on n = 1 .. ", NALL - r,
        "   (fit used n <= ", NFIT - r, ")");
  print("   failures = ", verify(seqs, cs, r, 1, NALL - r));
  print("   failures on HELD-OUT n = ", NFIT-r+1, "..", NALL-r, " : ",
        verify(seqs, cs, r, NFIT-r+1, NALL-r));
  print("  -- stage 4: singularity structure");
  \\ d_i = actual degree of c_i ; D = max ; the ODE for f(x)=sum u_n x^n is
  \\   sum_i x^(r-i) c_i(th - i) f = 0,  th = x d/dx.
  \\ Coefficient of th^D is  Pgf(x) = sum_{i : deg c_i = D} lead_i x^(r-i);
  \\ in d/dx form the leading coefficient is x^D * Pgf(x), so the finite singular
  \\ points of f are x = 0 and the roots of Pgf.
  my(dg = vector(r+1, i, poldegree(cs[i], 'nn)), dd = vecmax(dg), pgf = 0, npts);
  print("   deg c_i = ", dg);
  lead = vector(r+1, i, polcoef(cs[i], dg[i], 'nn));
  print("   leading coeff of each c_i = ", lead);
  for(i = 1, r+1, if(dg[i] == dd, pgf += lead[i] * 'x^(r-i+1)));
  print("   Pgf(x) = ", pgf, "   (D = ", dd, ")");
  \\ classical characteristic polynomial (only meaningful when all deg equal)
  if(vecmin(dg) == dd,
     chi = sum(i = 1, r+1, lead[i] * 'x^(i-1));
     print("   chi(x) = sum lead_i x^i = ", chi);
     print("   chi factored = ", factor(chi));
     rts = polroots(chi);
     print("   characteristic roots lambda = ", rts);
     print("   |lambda| = ", vector(#rts, t, abs(rts[t])));
     print("   log|lambda| = ", vector(#rts, t, log(abs(rts[t])))),
     print("   *** deg c_i NOT all equal: no classical characteristic polynomial.");
     print("   *** Newton polygon of (i, deg c_i):");
     for(i = 1, r, print("       edge ", i-1, "->", i, "  slope = ", dg[i+1]-dg[i],
                         "   => u_(n+1)/u_n ~ rho * n^", dg[i]-dg[i+1]));
     print("   *** balance polynomial sum_i lead_i rho^i (all i, formal) = ",
           sum(i = 1, r+1, lead[i]*'x^(i-1)));
     print("   *** its roots (rho) = ", polroots(sum(i = 1, r+1, lead[i]*'x^(i-1))));
     print("   *** |rho| = ", vector(r, t, abs(polroots(sum(i = 1, r+1, lead[i]*'x^(i-1)))[t]))));
  print("  -- generating-function singular points");
  if(poldegree(pgf, 'x) <= 0 || pgf == 0,
     print("   Pgf is constant: NO finite nonzero singular point."),
     sing = polroots(pgf);
     sing = select(z -> abs(z) > 1e-30, Vec(sing));
     print("   FINITE NONZERO SINGULAR POINTS = ", sing);
     print("   |sing| = ", vector(#sing, t, abs(sing[t])));
     print("   exact factorisation of Pgf = ", factor(pgf));
     print("   MIN modulus of nonzero finite singularity = ",
           vecmin(vector(#sing, t, abs(sing[t])))));
  print("   (x = 0 and x = infinity are always singular for these operators.)");
  1;
};
