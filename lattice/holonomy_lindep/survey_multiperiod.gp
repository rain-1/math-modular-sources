\\ ===================================================================
\\ survey_multiperiod.gp
\\ Census of multi-period rank-4 CY rows (degz >= 3) from mum_survey/ops.gp
\\   Step 1  filter degz >= 3
\\   Step 2  characteristic polynomial, lambda_1, lambda_2, simplicity
\\   Step 3  companion limits xi_j = lim X^(j)_n / A_n, certified digits,
\\           lindep -> period rank prk, identification vs. standard constants
\\   Step 4  sharp LCM denominator exponent k for each companion
\\ ===================================================================
default(parisizemax, 8G);
read("/home/ubuntu/code/math-modular-sources/lattice/mum_survey/ops.gp");
read("/home/ubuntu/code/math-modular-sources/lattice/multislope/sc_rows.gp");

\\ ---- basic constructors -------------------------------------------
\\ recurrence  sum_{i=0}^{degz} P_i(n-i) u_{n-i} = 0  ->  Cf[i+1](n) = P_i(n-i)
mkcf(op) = { my(dz=op[3], Ps=op[4]); vector(dz+1, i, subst(Ps[i], 'X, 'n-(i-1))); };
findop(cl) = { for(t=1,#OPS, if(OPS[t][2]==cl, return(OPS[t]))); 0; };

\\ chi(x) = sum_i (coeff of X^4 in P_i) x^{degz-i}
chipoly(op) = {
  my(dz = op[3], Ps = op[4], dmax = 0);
  for(i=1, #Ps, dmax = max(dmax, poldegree(Ps[i], 'X)));
  sum(i=0, dz, polcoeff(Ps[i+1], dmax, 'X) * 'x^(dz-i));
};

\\ ---- lindep with genuineness test ---------------------------------
prim(v) = { my(g=0); for(i=1,#v, g=gcd(g,v[i])); if(g==0, return(v)); v=v/g;
            my(i=1); while(i<=#v && v[i]==0, i++); if(i<=#v && v[i]<0, v=-v); v; };
\\ returns the relation if the SAME primitive vector is produced at two
\\ well-separated precisions, else 0 (spurious / no relation)
lidep(vals, D) = {
  my(dl = max(20, floor(0.45*D)), dh = max(30, floor(0.85*D)), pv, ph, sv);
  sv = default(realprecision);
  default(realprecision, dl); pv = prim(Vec(lindep(vector(#vals,i,vals[i]*1.0))));
  default(realprecision, dh); ph = prim(Vec(lindep(vector(#vals,i,vals[i]*1.0))));
  default(realprecision, sv);
  if(pv == ph, ph, 0);
};

\\ ---- constants for identification ----------------------------------
IDPREC = 300;
default(realprecision, IDPREC + 20);
CNAM = ["zeta(3)","zeta(5)","Pi^2","Pi^4","Catalan","L(chi-3,2)","L(chi-3,3)","L(chi5,3)","L(chi8,2)"];
CVAL = [zeta(3), zeta(5), Pi^2, Pi^4, lfun(-4,2), lfun(-3,2), lfun(-3,3), lfun(5,3), lfun(8,2)];

identify(xi, D) = {
  my(Did = min(D, IDPREC), out = List());
  my(r = lidep([xi, 1], Did));
  if(!(r === 0), listput(out, Str("RATIONAL: xi = ", -r[2], "/", r[1])));
  for(c=1, #CVAL,
    my(v = lidep([xi, 1, CVAL[c]], Did));
    if(!(v === 0) && v[3] != 0, listput(out, Str(v[1],"*xi + ",v[2]," + ",v[3],"*",CNAM[c]," = 0")));
  );
  for(c1=1, #CVAL, for(c2=c1+1, #CVAL,
    my(v = lidep([xi, 1, CVAL[c1], CVAL[c2]], Did));
    if(!(v === 0) && (v[3] != 0 || v[4] != 0),
       listput(out, Str(v[1],"*xi + ",v[2]," + ",v[3],"*",CNAM[c1]," + ",v[4],"*",CNAM[c2]," = 0")));
  ));
  Vec(out);
};

\\ ---- main analysis of one row ---------------------------------------
NTERM = 800; NBACK = 750; NDEN = 200; WPREC = 1600;
analyse(cl) = {
  my(op = findop(cl));
  if(op == 0, print("MISSING ", cl); return(0));
  my(dz = op[3], cf = mkcf(op), cp = chipoly(op));
  default(realprecision, 40);
  my(rr = polroots(cp), m = vector(#rr,i,abs(rr[i])), idx = vecsort(m,,5));
  my(l1 = rr[idx[1]], l2 = rr[idx[2]], a1 = m[idx[1]], a2 = m[idx[2]]);
  my(rat = a2/a1, l1simple = (abs(a2-a1) > 1e-20*a1));
  my(l1real = (abs(imag(l1)) < 1e-20*max(1,a1)), l2real = (abs(imag(l2)) < 1e-20*max(1,a2)));
  \\ sequences
  my(a = genseq(cf, dz, [1], NTERM));
  my(xs = vector(dz-1, j, compan(cf, dz, j, NTERM)));
  my(ks = vector(dz-1, j, denexp(xs[j], NDEN)));
  default(realprecision, WPREC);
  my(xi = vector(dz-1), cert = vector(dz-1));
  for(j=1, dz-1,
    my(rN = (xs[j][NTERM+1]*1.0)/a[NTERM+1], rB = (xs[j][NBACK+1]*1.0)/a[NBACK+1]);
    xi[j] = rN;
    my(d = abs(rN-rB));
    cert[j] = if(d == 0, WPREC-100, min(WPREC-100, floor(-log(d)/log(10))));
  );
  my(D = vecmin(cert), kmax = vecmax(ks));
  print("### AESZ ", op[1], "  cluster ", cl, "  degz=", dz);
  default(realprecision, 25);
  print("  chi(x) = ", cp);
  print("  lambda_1 = ", l1*1.0, "   lambda_2 = ", l2*1.0);
  print("  |l2/l1| = ", rat*1.0, "   l1 simple&real: ", l1simple, "/", l1real, "   l2 real: ", l2real);
  print("  sharp k per companion = ", ks, "   kmax = ", kmax);
  print("  -log|l2/l1| = ", -log(rat)*1.0, "   SCORE = ", (-log(rat) - kmax)*1.0);
  print("  certified digits = ", cert, "  (D = ", D, ")");
  default(realprecision, min(D, 60));
  for(j=1, dz-1, print("  xi_", j, " = ", xi[j]*1.0));
  \\ period rank: full relation lattice by iterated elimination
  my(vals = concat([1], vector(dz-1,j,xi[j])));
  my(labs = concat(["1"], vector(dz-1,j,Str("xi",j))));
  my(nrel = 0, H0 = 0);
  while(#vals > 1,
    my(r = lidep(vals, D));
    if(r === 0,
      default(realprecision, min(D, IDPREC));
      my(raw = Vec(lindep(vector(#vals,i,vals[i]*1.0))));
      H0 = vecmax(vector(#raw,i,abs(raw[i])));
      default(realprecision, WPREC);
      break);
    default(realprecision, min(D,IDPREC));
    my(res = abs(sum(i=1,#vals, r[i]*vals[i])));
    print("  GENUINE RELATION #", nrel+1, ": ", r, " on ", labs,
          "   height=", vecmax(vector(#r,i,abs(r[i]))),
          "  log10|res|=", if(res==0,-999,round(log(res)/log(10))));
    default(realprecision, WPREC);
    my(pp=0); for(i=1,#r, if(r[i]!=0, pp=i; break));
    my(keep = setminus(Set(vector(#vals,i,i)), Set([pp])));
    vals = vector(#keep,i,vals[keep[i]]); labs = vector(#keep,i,labs[keep[i]]);
    nrel++;
  );
  print("  #independent Q-relations on [1,xi_1..xi_",dz-1,"] = ", nrel,
        "   => prk = ", (dz-1) - nrel);
  if(nrel == 0 && H0 > 0,
     print("  (no relation of height below 10^", floor(log(H0)/log(10)),
           "; LLL certificate height ", H0, " at ", min(D,IDPREC), " digits)"));
  default(realprecision, IDPREC+20);
  for(j=1, dz-1,
    my(ids = identify(xi[j], D));
    if(#ids > 0, print("  ID xi_",j,": ", ids), print("  ID xi_",j,": none")));
  print("");
  1;
};

{TARGETS = ["4.5.46","4.3.18","4.3.9","4.5.82","4.3.1","4.8.2","4.5.111","4.8.1",
            "4.8.6","4.3.10","4.4.37","4.8.29","4.5.109","4.5.112","4.5.3","4.5.73",
            "4.8.3","4.4.50","4.3.28","4.3.29","4.3.8","4.8.10","4.5.33","4.3.22",
            "4.3.7","4.4.46","4.4.41","4.4.40","4.5.16","4.5.101","4.5.89"];}
           "4.8.6","4.4.46","4.3.7","4.8.1","4.8.2","4.4.41","4.3.10","4.4.40",
           "4.5.16","4.5.73","4.5.101","4.5.89"];}
{
foreach(TARGETS, cl, analyse(cl));
}
quit;
