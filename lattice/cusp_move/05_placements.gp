/* 05_placements.gp -- the three placements of every corpus row, over Q or over the
   quadratic field K = Q(sqrt(a^2-4d)); the gap picture; the formal best score;
   and the rigidity theorem  (rational roots  =>  |lambda_2| >= 1 at every placement). */
read("lib.gp"); read("rows.gp");
default(realprecision, 40);

{ show(rw) =
  my(aa = polcoef(rw[1],2,nv), dd = polcoef(rw[2],2,nv), ds, lam, mu, g, sc, gaps, best);
  ds = aa^2 - 4*dd;
  print(""); print("=== ", rw[3], "   ", rw[4]);
  print("   a=", aa, "  d=", dd, "  disc = a^2-4d = ", ds,
        if(ds < 0, "  (complex conjugate roots: no archimedean limit)",
           if(issquare(ds), "  (SQUARE: roots in Z, rational orbit)",
              Str("  (not a square: roots in Q(sqrt(", core(ds), ")))"))));
  if(ds < 0,
     lam = (aa + sqrt(ds*1.0))/2; mu = conj(lam);
     print("   |lambda| = |mu| = ", abs(lam)),
     lam = (aa + sqrt(ds*1.0))/2; mu = (aa - sqrt(ds*1.0))/2;
     gaps = [abs(lam), abs(mu), abs(lam-mu)];
     print("   config {0, lam, mu} = {0, ", lam, ", ", mu, "}");
     print("   the three gaps |lam|, |mu|, |lam-mu| = ", gaps);
     g = vecmin(gaps);
     print("   placement at v=0    : roots (", lam, ", ", mu, ")   |lam2| = ",
           min(abs(lam),abs(mu)),  "   score(k=2) = ", log(1/min(abs(lam),abs(mu)))-2);
     print("   placement at v=lam  : roots (", -lam, ", ", mu-lam, ")   |lam2| = ",
           min(abs(lam),abs(mu-lam)), "   score(k=2) = ", log(1/min(abs(lam),abs(mu-lam)))-2);
     print("   placement at v=mu   : roots (", -mu, ", ", lam-mu, ")   |lam2| = ",
           min(abs(mu),abs(lam-mu)), "   score(k=2) = ", log(1/min(abs(mu),abs(lam-mu)))-2);
     print("   best over the orbit : |lam2|_min = g = ", g, "   score = ", log(1/g)-2,
           if(issquare(ds), "   [rational: g is a positive integer, so score <= -2]",
              "   [attained only over K]"));
  );
}
for(i=1,#corpus, show(corpus[i]));

print(""); print("=== rigidity check: for every rational-root row in the corpus,");
print("    are all three gaps integers >= 1 ?");
{ for(i=1,#corpus,
  my(rw=corpus[i], aa=polcoef(rw[1],2,nv), dd=polcoef(rw[2],2,nv), ds, sq, lam, mu);
  ds = aa^2-4*dd;
  if(ds<=0 || !issquare(ds,&sq), next);
  lam = (aa+sq)/2; mu = (aa-sq)/2;
  print("   ", rw[3], ": lam=", lam, " mu=", mu, "  gaps=", [abs(lam),abs(mu),abs(lam-mu)],
        "  all integers>=1 : ", (denominator(lam)==1) && (denominator(mu)==1)
          && vecmin([abs(lam),abs(mu),abs(lam-mu)])>=1));
}
quit;
