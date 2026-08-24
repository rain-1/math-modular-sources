/* ==================================================================
 * 07_catalan.gp -- Catalan-focused identification battery for Apery
 *                  limits xi of four-term rows.
 *
 * INPUT FILE FORMAT: whitespace separated "label value", one per line,
 *                    value a real decimal with >= 60 correct digits.
 *
 * ------------------------------------------------------------------
 * INVOCATION (batch; default input lattice/four_term/out/xi.txt):
 *
 *   gp -q -s 2000000000 \
 *      /home/ubuntu/code/math-modular-sources/lattice/four_term_deep/07_catalan.gp
 *
 * INVOCATION on any other file of limits:
 *
 *   echo 'XIFILE="/abs/path/limits.txt"; read("/home/ubuntu/code/math-modular-sources/lattice/four_term_deep/07_catalan.gp")' \
 *     | gp -q -s 2000000000
 *
 * LIBRARY MODE (define everything, run nothing):
 *
 *   echo 'CATNORUN=1; read(".../07_catalan.gp"); print(catbattery(Catalan+zeta(2)))' \
 *     | gp -q -s 2000000000
 *
 * EXPORTED:
 *   catbattery(xv)      -> List() of verified hit-strings for one real xv
 *                          (T1 singles, T2 pairs, T4 algdep).
 *   catreport(lab, xv)  -> full printed report T1..T4 for one value.
 *   catrunfile(fname)   -> report on every line of fname.
 *   BATNAMES, BATVALS   -> the battery (labels, 250-digit values).
 *   PAIRIDX             -> greedily independent sublist used by T2.
 *
 * CONVENTIONS / GUARDS:
 *   * every candidate relation is VERIFIED by reconstructing xi from it
 *     and reporting the number of agreeing significant digits; only
 *     agreement > 55 digits at 60-digit input is reported as a hit.
 *   * every reported relation has a NONZERO coefficient on xi.
 *   * the T2 pair basis is reduced to a greedily independent subset at
 *     220 digits (height cap 10^6) BEFORE use -- see FOUR_TERM_SCAN.md
 *     section 6: a dependent basis makes lindep return relations among
 *     the basis elements instead of anything about xi.
 *   * PARI builtin names (psi, M, Phi, S, cmp, I, O, Pi, Catalan, ...)
 *     are never used as variables.
 * ================================================================== */

default(realprecision, 250);
PREC   = 60;      /* digits at which xi is trusted / lindep is run       */
PBIG   = 220;     /* digits for the basis-independence reduction         */
HCAP1  = 10^8;    /* T1 coefficient-height cap                           */
HCAP2  = 10^5;    /* T2 coefficient-height cap                           */
HCAPD  = 10^6;    /* independence-reduction height cap                   */
HCAPA  = 10^6;    /* T4 algdep height cap                                */
MINDIG = 55;      /* required verified agreement, in significant digits  */

{
if(type(XIFILE) != "t_STR",
   XIFILE = "/home/ubuntu/code/math-modular-sources/lattice/four_term/out/xi.txt");
}

/* ---------------- battery construction ---------------- */

BATNAMES = List(); BATVALS = List();

addc(nm, vv) =
{ my(rv = vv*1.0);
  for(i = 1, #BATVALS,
      if(abs(BATVALS[i] - rv) < 1.0e-100*(1+abs(rv)), return(0)));
  listput(BATNAMES, nm); listput(BATVALS, rv); 1;
}

GCAT = Catalan*1.0;
Z2   = zeta(2);
PI2  = Pi^2;
/* period of the elliptic-K3 row, consolidation/K3_ROW_PERIOD.md */
LG2  = (Pi/32)*gamma(1/8)*gamma(3/8)/(gamma(5/8)*gamma(7/8));

ti2(xx) = imag(polylog(2, xx*I));

{
  /* --- Catalan and its weight-shifted companions --- */
  addc("G", GCAT);
  addc("Pi*G",      Pi*GCAT);
  addc("G*sqrt(2)", GCAT*sqrt(2));
  addc("G*sqrt(3)", GCAT*sqrt(3));
  addc("G/sqrt(2)", GCAT/sqrt(2));
  addc("G/Pi",      GCAT/Pi);
  addc("G/Pi^2",    GCAT/PI2);

  /* --- zeta(2), Pi^2 and its real quadratic twists --- */
  addc("zeta(2)", Z2);
  addc("Pi^2",    PI2);
  foreach([1,2,3,4,6,8,12,24], dd,
    addc(Str("Pi^2*sqrt(", dd, ")"), PI2*sqrt(dd));
    addc(Str("Pi^2/sqrt(", dd, ")"), PI2/sqrt(dd)));

  /* --- L(2, chi_D), all fundamental discriminants |D| <= 24 --- */
  foreach([-3,-4,-7,-8,-11,-15,-19,-20,-23,-24,5,8,12,13,17,21,24], dd,
    addc(Str("L(2,chi_", dd, ")"), lfun(dd, 2)));

  /* --- weight three --- */
  foreach([-3,-4,-7,-8,5,8], dd,
    addc(Str("L(3,chi_", dd, ")"), lfun(dd, 3)));
  addc("zeta(3)",       zeta(3));
  addc("zeta(3)/Pi",    zeta(3)/Pi);
  addc("zeta(3)/Pi^2",  zeta(3)/PI2);

  /* --- inverse tangent integral Ti_2(x) = Im Li_2(i x) --- */
  foreach([1, 2, 1/2, 3, 1/3, 4, 1/4, 5, 1/5], xx,
    addc(Str("Ti2(", xx, ")"), ti2(xx)));
  addc("Ti2(2-sqrt(3))", ti2(2-sqrt(3)));
  addc("Ti2(sqrt(2)-1)", ti2(sqrt(2)-1));
  addc("G+(Pi/2)*log(2)  [=Ti2(2)]", GCAT + (Pi/2)*log(2));

  /* --- the elliptic-K3 row period --- */
  addc("L(g,2)",          LG2);
  addc("L(g,2)*sqrt(2)",  LG2*sqrt(2));
  addc("L(g,2)/sqrt(2)",  LG2/sqrt(2));
  addc("L(g,2)*Pi",       LG2*Pi);
  addc("L(g,2)/Pi",       LG2/Pi);

  /* --- CM period / Gamma products --- */
  addc("Gamma(1/4)^4/Pi",   gamma(1/4)^4/Pi);
  addc("Gamma(1/3)^6/Pi^2", gamma(1/3)^6/PI2);
  addc("G(1/8)G(3/8)/(G(5/8)G(7/8))",
       gamma(1/8)*gamma(3/8)/(gamma(5/8)*gamma(7/8)));
  addc("G(1/8)^2G(3/8)^2/Pi^3",
       gamma(1/8)^2*gamma(3/8)^2/Pi^3);
  addc("G(1/6)G(1/3)/(G(2/3)G(5/6))",
       gamma(1/6)*gamma(1/3)/(gamma(2/3)*gamma(5/6)));
  addc("G(1/12)G(5/12)/(G(7/12)G(11/12))",
       gamma(1/12)*gamma(5/12)/(gamma(7/12)*gamma(11/12)));
  addc("Gamma(1/6)^6/Pi^2", gamma(1/6)^6/PI2);

  /* --- logarithmic weight two --- */
  foreach([2,3,5,6,8,9,10,12,16,18,24,27,32], mm,
    addc(Str("log(", mm, ")^2"), log(mm)^2));
  foreach([2,3,5,6,8,9,10,12,16,18,24,27,32], mm,
    addc(Str("Pi*log(", mm, ")"), Pi*log(mm)));
  foreach([2,3,5,6,7,8,10,12], m1, foreach([2,3,5,6,7,8,10,12], m2,
    if(m1 < m2, addc(Str("log(", m1, ")*log(", m2, ")"), log(m1)*log(m2)))));

  /* --- real quadratic units --- */
  foreach([2,3,5,8,12], mm,
    my(uu = (sqrt(mm)+1)/abs(sqrt(mm)-1));
    addc(Str("Pi*log((sqrt(", mm, ")+1)/(sqrt(", mm, ")-1))"), Pi*log(uu));
    addc(Str("log((sqrt(", mm, ")+1)/(sqrt(", mm, ")-1))^2"), log(uu)^2));

  /* --- arctangent weight two --- */
  foreach([1,2,3,4,5], mm, addc(Str("Pi*atan(1/", mm, ")"), Pi*atan(1/mm)));
  foreach([1,2,3,4,5], mm, addc(Str("atan(1/", mm, ")^2"), atan(1/mm)^2));
}

NBAT = #BATVALS;

/* ---------------- helpers ---------------- */

PR(vv, dg) = apply(zz -> precision(zz*1.0, dg), vv);

digagree(aa, bb) =
{ my(dd = abs(aa-bb), sc = abs(bb));
  if(sc == 0, sc = 1);
  if(dd == 0, return(1000));
  -log(dd/sc)/log(10);
}

/* T1: 1, xi, X */
try1(xv, ii) =
{ my(rel, hh, rec, dg);
  rel = lindep(PR([1, xv, BATVALS[ii]], PREC));
  if(type(rel) != "t_COL" && type(rel) != "t_VEC", return(0));
  if(#rel != 3, return(0));
  if(rel[2] == 0 || rel[3] == 0, return(0));
  hh = vecmax(apply(abs, rel));
  if(hh > HCAP1, return(0));
  rec = -(rel[1] + rel[3]*BATVALS[ii])/rel[2];
  dg  = digagree(rec, xv);
  if(dg <= MINDIG, return(0));
  [Str("xi = (", -rel[1], " + ", -rel[3], "*[", BATNAMES[ii], "])/", rel[2]),
   rel, dg];
}

/* T2: 1, xi, X, Y */
try2(xv, ii, jj) =
{ my(rel, hh, rec, dg);
  rel = lindep(PR([1, xv, BATVALS[ii], BATVALS[jj]], PREC));
  if(#rel != 4, return(0));
  if(rel[2] == 0, return(0));
  if(rel[3] == 0 || rel[4] == 0, return(0));   /* genuine pairs only; the
        single-constant cases are already covered by T1, which uses the whole
        145-constant battery rather than the 22-element pair basis */
  hh = vecmax(apply(abs, rel));
  if(hh > HCAP2, return(0));
  rec = -(rel[1] + rel[3]*BATVALS[ii] + rel[4]*BATVALS[jj])/rel[2];
  dg  = digagree(rec, xv);
  if(dg <= MINDIG, return(0));
  [Str("xi = (", -rel[1], " + ", -rel[3], "*[", BATNAMES[ii], "] + ",
       -rel[4], "*[", BATNAMES[jj], "])/", rel[2]), rel, dg];
}

/* ---------------- T2 basis: greedy independent reduction ---------------- */

PAIRCAND = List();
{
  listput(PAIRCAND, "G");
  listput(PAIRCAND, "Pi*G");
  listput(PAIRCAND, "zeta(2)");
  listput(PAIRCAND, "Pi^2");
  foreach([-3,-4,-7,-8,-11,-15,-19,-20,-23,-24,5,8,12,13,17,21,24], dd,
          listput(PAIRCAND, Str("L(2,chi_", dd, ")")));
  listput(PAIRCAND, "L(g,2)");
  listput(PAIRCAND, "log(2)^2");
  listput(PAIRCAND, "Pi*log(2)");
}

nameidx(nm) =
{ for(i = 1, NBAT, if(BATNAMES[i] == nm, return(i))); 0; }

PAIRIDX = List(); PAIRDROP = List();
{
  my(basis = [], rel, hh, res, ok);
  foreach(Vec(PAIRCAND), nm,
    my(ii = nameidx(nm));
    if(ii == 0, listput(PAIRDROP, Str(nm, "  [duplicate value, deduped]")),
      ok = 1;
      rel = lindep(PR(concat(concat([1], basis), [BATVALS[ii]]), PBIG));
      if(type(rel) == "t_COL" && #rel == #basis + 2 && rel[#rel] != 0,
         hh = vecmax(apply(abs, rel));
         res = abs(rel[1] + sum(k = 1, #basis, rel[k+1]*basis[k])
                          + rel[#rel]*BATVALS[ii]);
         if(hh <= HCAPD && res < 1.0e-160, ok = 0;
            listput(PAIRDROP, Str(nm, "  [dependent, rel = ", rel, "]"))));
      if(ok, basis = concat(basis, [BATVALS[ii]]); listput(PAIRIDX, ii)));
  );
}

/* ---------------- the four tests ---------------- */

catbattery(xv) =
{ my(hits = List(), r, pv, pl);
  for(i = 1, NBAT,
    r = try1(xv, i);
    if(r != 0, listput(hits, Str("T1  ", r[1], "   [", r[2], "]  agree ",
                                 floor(r[3]), " digits"))));
  pv = Vec(PAIRIDX);
  for(a = 1, #pv, for(b = a+1, #pv,
    r = try2(xv, pv[a], pv[b]);
    if(r != 0, listput(hits, Str("T2  ", r[1], "   [", r[2], "]  agree ",
                                 floor(r[3]), " digits")))));
  for(k = 1, 6,
    pl = algdep(precision(xv*1.0, PREC), k);
    if(poldegree(pl) >= 1,
      my(hh = vecmax(apply(abs, Vec(pl))), rts = polroots(pl), bd = 0, dg);
      for(m = 1, #rts,
          dg = if(abs(imag(rts[m])) > 1.0e-40, 0, digagree(real(rts[m]), xv));
          if(dg > bd, bd = dg));
      if(hh <= HCAPA && bd > MINDIG,
         listput(hits, Str("T4  algdep deg<=", k, ":  ", pl,
                           "   height ", hh, "  agree ", floor(bd), " digits")))));
  hits;
}

/* T3: the explicit Catalan shapes, always printed */
cat_t3(lab, xv) =
{ my(r1, r2);
  print("  T3 ratios:  xi/G       = ", precision(xv/GCAT, 40));
  print("              xi/(Pi*G)  = ", precision(xv/(Pi*GCAT), 40));
  print("              xi*G       = ", precision(xv*GCAT, 40));
  print("              xi/zeta(2) = ", precision(xv/Z2, 40));
  print("              xi/Pi^2    = ", precision(xv/PI2, 40));
  print("              xi/L(g,2)  = ", precision(xv/LG2, 40));
  r1 = lindep(PR([1, xv, GCAT, Z2], PREC));
  r2 = lindep(PR([xv, GCAT, Z2], PREC));
  print("  T3 lindep([1,xi,G,z2]) = ", r1,
        "    residual ", precision(abs(r1[1] + r1[2]*xv + r1[3]*GCAT + r1[4]*Z2), 5));
  print("  T3 lindep([xi,G,z2])   = ", r2,
        "    residual ", precision(abs(r2[1]*xv + r2[2]*GCAT + r2[3]*Z2), 5));
}

/* T4 raw dump (unverified), always printed */
cat_t4(lab, xv) =
{ my(pl, rt);
  for(k = 1, 6,
    pl = algdep(precision(xv*1.0, PREC), k);
    print("  T4 algdep(xi,", k, ") = ", pl,
          "   min|root-xi| = ",
          if(poldegree(pl) < 1, "n/a",
             precision(vecmin(apply(zz -> abs(zz - xv), polroots(pl))), 5))));
  foreach([["xi/G", xv/GCAT], ["xi/zeta(2)", xv/Z2], ["xi/Pi^2", xv/PI2],
           ["xi/L(g,2)", xv/LG2]], pr,
    pl = algdep(precision(pr[2], PREC), 4);
    print("  T4 algdep(", pr[1], ",4) = ", pl,
          "   min|root-r| = ",
          if(poldegree(pl) < 1, "n/a",
             precision(vecmin(apply(zz -> abs(zz - pr[2]), polroots(pl))), 5))));
}

catreport(lab, xv) =
{ my(hits, rr);
  print("");
  print("=== ", lab);
  print("  xi = ", precision(xv, 70));
  rr = lindep(PR([1, xv], PREC));
  if(type(rr) == "t_COL" && #rr == 2 && rr[2] != 0 && vecmax(apply(abs, rr)) < 10^12
     && digagree(-rr[1]/rr[2], xv) > MINDIG,
     print("  RATIONAL: xi = ", -rr[1], "/", rr[2]));
  hits = catbattery(xv);
  if(#hits == 0,
     print("  T1/T2/T4: NO HIT (unidentified)"),
     foreach(Vec(hits), h, print("  HIT  ", h)));
  cat_t3(lab, xv);
  cat_t4(lab, xv);
  hits;
}

catrunfile(fname) =
{ my(ln = readstr(fname), nh = 0, labs = List());
  for(i = 1, #ln,
    my(vv = select(zz -> zz != "", strsplit(ln[i], " ")));
    if(#vv >= 2,
       my(hh = catreport(vv[1], eval(Str(vv[2], "*1.0"))));
       if(#hh > 0, nh = nh + #hh; listput(labs, vv[1]))));
  print("");
  print("################ SUMMARY ################");
  print("battery size          : ", NBAT, " constants");
  print("T2 independent basis  : ", #PAIRIDX, " constants, ",
        binomial(#PAIRIDX,2), " pairs per xi");
  print("input file            : ", fname);
  print("lines processed       : ", #ln);
  print("total verified hits   : ", nh);
  if(nh == 0, print("=> ALL LIMITS UNIDENTIFIED against the Catalan battery."),
              print("=> rows with hits: ", Vec(labs)));
  nh;
}

/* ---------------- main ---------------- */

{
if(type(CATNORUN) != "t_INT",
  print("07_catalan.gp  -- Catalan-focused identification battery");
  print("realprecision = ", default(realprecision), ",  lindep at ", PREC,
        " digits,  verification threshold ", MINDIG, " digits");
  print("battery: ", NBAT, " constants");
  print("T2 pair candidates: ", #PAIRCAND, "; dropped as dependent/duplicate: ",
        #PAIRDROP);
  foreach(Vec(PAIRDROP), dnm, print("   drop  ", dnm));
  print("T2 independent basis (", #PAIRIDX, "): ",
        Vec(apply(zz -> BATNAMES[zz], Vec(PAIRIDX))));
  print("");
  print("---- battery listing ----");
  for(i = 1, NBAT, print("  [", i, "] ", BATNAMES[i], " = ",
                         precision(BATVALS[i], 30)));
  catrunfile(XIFILE);
  quit;
);
}
